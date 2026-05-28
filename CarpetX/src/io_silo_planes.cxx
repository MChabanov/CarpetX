#include "io_silo_planes.hxx"

#include "driver.hxx"
#include "io_meta.hxx"
#include "io_planes.hxx"
#include "mpi_types.hxx"
#include "timer.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef HAVE_CAPABILITY_Silo

#include <AMReX.H>
#include <AMReX_IntVect.H>

#include <mpi.h>

#include <silo.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <memory>
#include <mutex>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace CarpetX {

namespace {

constexpr int planes_ndims = 2;

template <typename T> struct db_datatype;
template <> struct db_datatype<char> : std::integral_constant<int, DB_CHAR> {};
template <> struct db_datatype<int> : std::integral_constant<int, DB_INT> {};
template <> struct db_datatype<long> : std::integral_constant<int, DB_LONG> {};
template <>
struct db_datatype<long long> : std::integral_constant<int, DB_LONG_LONG> {};
template <>
struct db_datatype<float> : std::integral_constant<int, DB_FLOAT> {};
template <>
struct db_datatype<double> : std::integral_constant<int, DB_DOUBLE> {};
template <typename T> constexpr int db_datatype_v = db_datatype<T>::value;

inline std::array<int, 2> in_plane_axes(int normal_axis) {
  switch (normal_axis) {
  case 0:
    return {1, 2};
  case 1:
    return {0, 2};
  case 2:
    return {0, 1};
  }
  assert(0);
  return {0, 1};
}

struct plane_mesh_props_t {
  std::array<int, 2> nghosts;
  friend bool operator<(const plane_mesh_props_t &p,
                        const plane_mesh_props_t &q) {
    return p.nghosts < q.nghosts;
  }
};

std::string make_plane_subdirname(const std::string &file_name,
                                  const std::string &plane_tag, int iteration) {
  std::ostringstream buf;
  buf << file_name << "." << plane_tag << ".it" << std::setw(8)
      << std::setfill('0') << iteration << ".silo_planes.dir";
  return buf.str();
}

std::string make_plane_filename(const std::string &file_name,
                                const std::string &plane_tag, int iteration,
                                int ioserver = -1) {
  std::ostringstream buf;
  buf << file_name << "." << plane_tag << ".it" << std::setw(8)
      << std::setfill('0') << iteration;
  if (ioserver >= 0)
    buf << ".p" << std::setw(6) << std::setfill('0') << ioserver;
  buf << ".silo";
  return buf.str();
}

std::string make_plane_meshname(const std::string &plane_tag,
                                const std::array<int, 2> &nghosts,
                                int patch = -1, int reflevel = -1,
                                int component = -1) {
  assert((patch == -1) == (reflevel == -1));
  assert((patch == -1) == (component == -1));
  std::ostringstream buf;
  buf << (patch >= 0 ? "box" : "gh") << "." << plane_tag << ".ghosts";
  for (int d = 0; d < 2; ++d) {
    if (d > 0)
      buf << "_";
    buf << std::setw(2) << std::setfill('0') << nghosts[d];
  }
  if (patch >= 0)
    buf << ".m" << std::setw(4) << std::setfill('0') << patch << ".rl"
        << std::setw(2) << std::setfill('0') << reflevel << ".c" << std::setw(8)
        << std::setfill('0') << component;
  return DB::legalize_name(buf.str());
}

std::string make_plane_varname(int gi, int vi, const std::string &plane_tag,
                               int patch = -1, int reflevel = -1,
                               int component = -1) {
  assert((patch == -1) == (reflevel == -1));
  assert((patch == -1) == (component == -1));
  std::string varname;
  {
    assert(vi >= 0);
    const int v0 = CCTK_FirstVarIndexI(gi);
    varname = CCTK_FullVarName(v0 + vi);
    varname = std::regex_replace(varname, std::regex("::"), "-");
    for (auto &ch : varname)
      ch = std::tolower(static_cast<unsigned char>(ch));
  }
  std::ostringstream buf;
  buf << varname << "." << plane_tag;
  if (reflevel >= 0)
    buf << ".m" << std::setw(4) << std::setfill('0') << patch << ".rl"
        << std::setw(2) << std::setfill('0') << reflevel << ".c" << std::setw(8)
        << std::setfill('0') << component;
  return DB::legalize_name(buf.str());
}

} // namespace

void OutputSiloPlanes(const cGH *const cctkGH,
                      const std::vector<bool> &output_group,
                      const std::vector<plane_spec_t> &planes,
                      const std::string &output_dir,
                      const std::string &output_file) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  static Timer timer("OutputSiloPlanes");
  Interval interval(timer);

  if (planes.empty())
    return;
  if (std::count(output_group.begin(), output_group.end(), true) == 0)
    return;

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(cctkGH);
  const int nprocs = CCTK_nProcs(cctkGH);

  const int ioproc_every = [&]() {
    if (CCTK_EQUALS(out_mode, "proc"))
      return 1;
    if (CCTK_EQUALS(out_mode, "np"))
      return out_proc_every;
    if (CCTK_EQUALS(out_mode, "onefile"))
      return nprocs;
    assert(0);
  }();
  assert(ioproc_every > 0);

  const bool write_file = myproc % ioproc_every == 0;
  const int metafile_ioproc = nprocs == 1 || ioproc_every == 1 ? 0 : 1;
  const bool write_metafile = myproc == metafile_ioproc;

  DBShowErrors(DB_ALL_AND_DRVR, nullptr);
  DBSetCompression(out_silo_compression_options);
  DBSetEnableChecksums(1);

  static std::set<int> warned_noncart_patches;
  static std::set<int> warned_skipped_groups;
  static std::set<std::string> warned_outside_domain;

  static std::once_flag create_out_dir;
  std::call_once(create_out_dir, [&]() {
    const int rc = CCTK_CreateDirectory(0755, output_dir.c_str());
    assert(rc >= 0);
  });

  for (const auto &plane : planes) {
    const auto in_axes = in_plane_axes(plane.normal_axis);
    const int axis_a = in_axes[0], axis_b = in_axes[1];

    const std::string subdirname =
        make_plane_subdirname(output_file, plane.tag, cctk_iteration);
    const std::string pathname = output_dir + "/" + subdirname;
    ierr = CCTK_CreateDirectory(0755, pathname.c_str());
    assert(ierr >= 0);

    DB::ptr<DBfile> file;
    if (write_file) {
      const std::string filename =
          pathname + "/" +
          make_plane_filename(output_file, plane.tag, cctk_iteration,
                              myproc / ioproc_every);
      file = DB::make(DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL,
                               output_file.c_str(), DB_HDF5));
      assert(file);
    }

    bool any_slab_emitted = false;

    for (const auto &patchdata : ghext->patchdata) {
      if (!patchdata.is_cartesian) {
        if (warned_noncart_patches.insert(patchdata.patch).second)
          CCTK_VWARN(CCTK_WARN_ALERT,
                     "OutputSiloPlanes: skipping non-Cartesian patch %d",
                     patchdata.patch);
        continue;
      }

      for (const auto &leveldata : patchdata.leveldata) {
        const amrex::Geometry &geom = patchdata.amrcore->Geom(leveldata.level);
        const amrex::Real *const x0 = geom.ProbLo();
        const amrex::Real *const dx = geom.CellSize();

        std::set<plane_mesh_props_t> have_meshes;
        for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
          if (!output_group.at(gi))
            continue;
          if (CCTK_GroupTypeI(gi) != CCTK_GF)
            continue;

          const auto &groupdata = *leveldata.groupdata.at(gi);
          if (groupdata.mfab.empty())
            continue;
          const int numvars = groupdata.numvars;
          const int tl = 0;
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::IndexType &indextype = mfab.ixType();
          const amrex::DistributionMapping &dm = mfab.DistributionMap();

          const int rank3 = int(indextype.cellCentered(0)) +
                            int(indextype.cellCentered(1)) +
                            int(indextype.cellCentered(2));
          if (rank3 != 0 && rank3 != 3) {
            if (warned_skipped_groups.insert(gi).second)
              CCTK_VWARN(CCTK_WARN_ALERT,
                         "OutputSiloPlanes: skipping staggered group %s "
                         "(rank-%d); supported in a follow-up commit",
                         CCTK_FullGroupName(gi), rank3);
            continue;
          }

          const int slice_idx = snap_to_grid_index(plane, geom, indextype);
          if (slice_idx < 0)
            continue;

          const bool in_plane_cc = (rank3 == 3);
          const int centering = in_plane_cc ? DB_ZONECENT : DB_NODECENT;

          const std::array<int, dim> &nghosts3 = groupdata.nghostzones;
          const plane_mesh_props_t mesh_props{
              {nghosts3[axis_a], nghosts3[axis_b]}};
          const bool have_mesh = have_meshes.count(mesh_props);

          const int ncomponents = dm.size();
          for (int component = 0; component < ncomponents; ++component) {
            const int proc = dm[component];
            const int ioproc = proc / ioproc_every * ioproc_every;
            const bool send_this_fab = proc == myproc;
            const bool write_this_fab = ioproc == myproc;
            if (!(send_this_fab || write_this_fab))
              continue;

            const amrex::Box &fabbox = mfab.fabbox(component);
            if (slice_idx < fabbox.smallEnd(plane.normal_axis) ||
                slice_idx > fabbox.bigEnd(plane.normal_axis))
              continue;

            const std::array<int, 2> dims = {fabbox.length(axis_a),
                                             fabbox.length(axis_b)};
            const std::ptrdiff_t zonecount = std::ptrdiff_t(dims[0]) * dims[1];
            assert(numvars * zonecount <= INT_MAX);

            std::vector<CCTK_REAL> local_buf;
            std::vector<CCTK_REAL> recv_buf;
            const CCTK_REAL *data = nullptr;

            if (send_this_fab) {
              const amrex::FArrayBox &fab = mfab[component];
              extract_slab(fab, plane.normal_axis, slice_idx, numvars,
                           local_buf);
              if (write_this_fab) {
                data = local_buf.data();
              } else {
                const int mpi_tag = 22910;
                MPI_Send(local_buf.data(), int(numvars * zonecount),
                         mpi_datatype_v<CCTK_REAL>, ioproc, mpi_tag, mpi_comm);
              }
            } else {
              assert(write_this_fab);
              const int mpi_tag = 22910;
              recv_buf.resize(std::size_t(numvars) * zonecount);
              MPI_Recv(recv_buf.data(), int(numvars * zonecount),
                       mpi_datatype_v<CCTK_REAL>, proc, mpi_tag, mpi_comm,
                       MPI_STATUS_IGNORE);
              data = recv_buf.data();
            }

            if (!write_this_fab)
              continue;

            any_slab_emitted = true;

            if (!have_mesh) {
              const std::array<int, 2> dims_vc = {
                  dims[0] + int(indextype.cellCentered(axis_a)),
                  dims[1] + int(indextype.cellCentered(axis_b))};

              std::array<std::vector<CCTK_REAL>, 2> coords;
              std::array<const void *, 2> coord_ptrs;
              for (int d = 0; d < 2; ++d) {
                const int ax = (d == 0) ? axis_a : axis_b;
                coords[d].resize(dims_vc[d]);
                for (int i = 0; i < dims_vc[d]; ++i)
                  coords[d][i] = x0[ax] + (fabbox.smallEnd(ax) + i) * dx[ax];
                coord_ptrs[d] = coords[d].data();
              }

              const std::string meshname = make_plane_meshname(
                  plane.tag, mesh_props.nghosts, patchdata.patch,
                  leveldata.level, component);

              const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
              assert(optlist);
              int cartesian = DB_CARTESIAN;
              ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
              assert(!ierr);
              int cycle = cctk_iteration;
              ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
              assert(!ierr);
              std::array<int, 2> min_index = {nghosts3[axis_a],
                                              nghosts3[axis_b]};
              std::array<int, 2> max_index = min_index;
              ierr =
                  DBAddOption(optlist.get(), DBOPT_LO_OFFSET, min_index.data());
              assert(!ierr);
              ierr =
                  DBAddOption(optlist.get(), DBOPT_HI_OFFSET, max_index.data());
              assert(!ierr);
              int column_major = 0;
              ierr =
                  DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
              assert(!ierr);
              double dtime = cctk_time;
              ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
              assert(!ierr);
              int hide_from_gui = 1;
              ierr = DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI,
                                 &hide_from_gui);
              assert(!ierr);

              ierr = DBPutQuadmesh(file.get(), meshname.c_str(), nullptr,
                                   coord_ptrs.data(), dims_vc.data(),
                                   planes_ndims, db_datatype_v<CCTK_REAL>,
                                   DB_COLLINEAR, optlist.get());
              assert(!ierr);
            }

            const std::string meshname = make_plane_meshname(
                plane.tag, mesh_props.nghosts, patchdata.patch, leveldata.level,
                component);

            const DB::ptr<DBoptlist> var_optlist = DB::make(DBMakeOptlist(10));
            assert(var_optlist);
            int cartesian = DB_CARTESIAN;
            ierr = DBAddOption(var_optlist.get(), DBOPT_COORDSYS, &cartesian);
            assert(!ierr);
            int cycle = cctk_iteration;
            ierr = DBAddOption(var_optlist.get(), DBOPT_CYCLE, &cycle);
            assert(!ierr);
            int column_major = 0;
            ierr =
                DBAddOption(var_optlist.get(), DBOPT_MAJORORDER, &column_major);
            assert(!ierr);
            double dtime = cctk_time;
            ierr = DBAddOption(var_optlist.get(), DBOPT_DTIME, &dtime);
            assert(!ierr);
            int hide_from_gui = 1;
            ierr = DBAddOption(var_optlist.get(), DBOPT_HIDE_FROM_GUI,
                               &hide_from_gui);
            assert(!ierr);

            for (int vi = 0; vi < numvars; ++vi) {
              const std::string varname =
                  make_plane_varname(gi, vi, plane.tag, patchdata.patch,
                                     leveldata.level, component);
              const void *const data_ptr = data + vi * zonecount;
              ierr = DBPutQuadvar1(
                  file.get(), varname.c_str(), meshname.c_str(), data_ptr,
                  dims.data(), planes_ndims, nullptr, 0,
                  db_datatype_v<CCTK_REAL>, centering, var_optlist.get());
              assert(!ierr);
            }
          }
          have_meshes.insert(mesh_props);
        }
      }
    }

    if (!any_slab_emitted) {
      const std::string key = plane.tag;
      if (warned_outside_domain.insert(key).second)
        CCTK_VWARN(CCTK_WARN_ALERT,
                   "OutputSiloPlanes: plane %s lies outside all Cartesian "
                   "(patch, level) extents on this iteration",
                   plane.tag.c_str());
    }

    if (write_metafile) {
      const std::string metafilename =
          output_dir + "/" +
          make_plane_filename(output_file, plane.tag, cctk_iteration);
      const DB::ptr<DBfile> metafile =
          DB::make(DBCreate(metafilename.c_str(), DB_CLOBBER, DB_LOCAL,
                            output_file.c_str(), DB_HDF5));
      assert(metafile);

      {
        const int mdims = 1;
        const int value = 1;
        ierr = DBWrite(metafile.get(), "MetadataIsTimeVarying", &value, &mdims,
                       1, DB_INT);
        assert(!ierr);
      }

      std::set<plane_mesh_props_t> meta_have_meshes;
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        if (!output_group.at(gi))
          continue;
        if (CCTK_GroupTypeI(gi) != CCTK_GF)
          continue;

        const auto &patchdata0 = ghext->patchdata.at(0);
        const auto &leveldata0 = patchdata0.leveldata.at(0);
        const auto &groupdata0 = *leveldata0.groupdata.at(gi);
        if (groupdata0.mfab.empty())
          continue;
        const int numvars = groupdata0.numvars;
        const int tl = 0;
        const amrex::MultiFab &mfab0 = *groupdata0.mfab[tl];
        const amrex::IndexType &indextype = mfab0.ixType();

        const int rank3 = int(indextype.cellCentered(0)) +
                          int(indextype.cellCentered(1)) +
                          int(indextype.cellCentered(2));
        if (rank3 != 0 && rank3 != 3)
          continue;

        const std::array<int, dim> &nghosts3 = groupdata0.nghostzones;
        const plane_mesh_props_t mesh_props{
            {nghosts3[axis_a], nghosts3[axis_b]}};

        std::vector<std::string> meshnames;
        std::vector<std::vector<std::string> > varnames_per_var(numvars);

        for (const auto &patchdata : ghext->patchdata) {
          if (!patchdata.is_cartesian)
            continue;
          for (const auto &leveldata : patchdata.leveldata) {
            const amrex::Geometry &geom =
                patchdata.amrcore->Geom(leveldata.level);
            const int slice_idx = snap_to_grid_index(plane, geom, indextype);
            if (slice_idx < 0)
              continue;

            const auto &gdata = *leveldata.groupdata.at(gi);
            if (gdata.mfab.empty())
              continue;
            const amrex::MultiFab &mfab = *gdata.mfab[tl];
            const amrex::DistributionMapping &dm = mfab.DistributionMap();
            const int ncomponents = dm.size();
            for (int c = 0; c < ncomponents; ++c) {
              const amrex::Box &fabbox = mfab.fabbox(c);
              if (slice_idx < fabbox.smallEnd(plane.normal_axis) ||
                  slice_idx > fabbox.bigEnd(plane.normal_axis))
                continue;

              const int proc = dm[c];
              const std::string proc_filename =
                  make_plane_subdirname(output_file, plane.tag,
                                        cctk_iteration) +
                  "/" +
                  make_plane_filename(output_file, plane.tag, cctk_iteration,
                                      proc / ioproc_every);
              const std::string meshname =
                  proc_filename + ":" +
                  make_plane_meshname(plane.tag, mesh_props.nghosts,
                                      patchdata.patch, leveldata.level, c);
              meshnames.push_back(meshname);

              for (int vi = 0; vi < numvars; ++vi) {
                const std::string varname =
                    proc_filename + ":" +
                    make_plane_varname(gi, vi, plane.tag, patchdata.patch,
                                       leveldata.level, c);
                varnames_per_var[vi].push_back(varname);
              }
            }
          }
        }

        if (meshnames.empty())
          continue;

        if (!meta_have_meshes.count(mesh_props)) {
          const std::string multimeshname =
              make_plane_meshname(plane.tag, mesh_props.nghosts);
          std::vector<const char *> meshname_ptrs;
          meshname_ptrs.reserve(meshnames.size());
          for (const auto &s : meshnames)
            meshname_ptrs.push_back(s.c_str());

          const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
          assert(optlist);
          int cycle = cctk_iteration;
          ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
          assert(!ierr);
          double dtime = cctk_time;
          ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
          assert(!ierr);
          int quadmesh = DB_QUADMESH;
          ierr = DBAddOption(optlist.get(), DBOPT_MB_BLOCK_TYPE, &quadmesh);
          assert(!ierr);

          ierr = DBPutMultimesh(metafile.get(), multimeshname.c_str(),
                                int(meshname_ptrs.size()), meshname_ptrs.data(),
                                nullptr, optlist.get());
          assert(!ierr);
          meta_have_meshes.insert(mesh_props);
        }

        const std::string multimeshname =
            make_plane_meshname(plane.tag, mesh_props.nghosts);
        const DB::ptr<DBoptlist> mv_optlist = DB::make(DBMakeOptlist(10));
        assert(mv_optlist);
        int cycle = cctk_iteration;
        ierr = DBAddOption(mv_optlist.get(), DBOPT_CYCLE, &cycle);
        assert(!ierr);
        double dtime = cctk_time;
        ierr = DBAddOption(mv_optlist.get(), DBOPT_DTIME, &dtime);
        assert(!ierr);
        ierr = DBAddOption(mv_optlist.get(), DBOPT_MMESH_NAME,
                           const_cast<char *>(multimeshname.c_str()));
        assert(!ierr);
        int vartype_scalar = DB_VARTYPE_SCALAR;
        ierr =
            DBAddOption(mv_optlist.get(), DBOPT_TENSOR_RANK, &vartype_scalar);
        assert(!ierr);
        int quadvar = DB_QUADVAR;
        ierr = DBAddOption(mv_optlist.get(), DBOPT_MB_BLOCK_TYPE, &quadvar);
        assert(!ierr);

        for (int vi = 0; vi < numvars; ++vi) {
          const std::string multivarname =
              make_plane_varname(gi, vi, plane.tag);
          std::vector<const char *> varname_ptrs;
          varname_ptrs.reserve(varnames_per_var[vi].size());
          for (const auto &s : varnames_per_var[vi])
            varname_ptrs.push_back(s.c_str());
          ierr = DBPutMultivar(metafile.get(), multivarname.c_str(),
                               int(varname_ptrs.size()), varname_ptrs.data(),
                               nullptr, mv_optlist.get());
          assert(!ierr);
        }
      }

      {
        const std::string visitname = output_dir + "/" + output_file + "." +
                                      plane.tag + ".silo_planes.visit";
        std::ofstream visit(visitname, std::ios::app);
        assert(visit.good());
        visit << make_plane_filename(output_file, plane.tag, cctk_iteration)
              << "\n";
      }

      {
        output_file_description_t ofd;
        ofd.filename = metafilename;
        ofd.description =
            "2D plane CarpetX HDF5 Silo output (" + plane.tag + ")";
        ofd.writer_thorn = CCTK_THORNSTRING;
        for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
          if (!output_group.at(gi))
            continue;
          if (CCTK_GroupTypeI(gi) != CCTK_GF)
            continue;
          const int numvars = CCTK_NumVarsInGroupI(gi);
          assert(numvars >= 0);
          const int firstvar = CCTK_FirstVarIndexI(gi);
          for (int vi = 0; vi < numvars; ++vi)
            ofd.variables.push_back(CCTK_FullVarName(firstvar + vi));
        }
        ofd.iterations = {cctk_iteration};
        ofd.output_directions.push_back(axis_a);
        ofd.output_directions.push_back(axis_b);
        ofd.format_name = "CarpetX/Silo/HDF5";
        ofd.format_version = {1, 0, 0};
        OutputMeta_RegisterOutputFile(std::move(ofd));
      }
    }
  }
}

} // namespace CarpetX

#endif // HAVE_CAPABILITY_Silo
