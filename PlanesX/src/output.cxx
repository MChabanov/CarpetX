#include "carpetx_params.hxx"
#include "openpmd_planes.hxx"
#include "planes.hxx"
#include "silo_planes.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstring>
#include <string>
#include <vector>

namespace PlanesX {

namespace {

// These mirror CarpetX's io.cxx helpers, which are file-local there.

std::vector<bool> find_groups(const char *const method,
                              const char *const out_vars) {
  DECLARE_CCTK_PARAMETERS;

  std::vector<bool> enabled(CCTK_NumGroups(), false);
  const auto callback{
      [](const int index, const char *const optstring, void *const arg) {
        std::vector<bool> &enabled = *static_cast<std::vector<bool> *>(arg);
        enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
      }};
  CCTK_TraverseString(out_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
  if (verbose) {
    CCTK_VINFO("%s output for groups:", method);
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
      if (enabled.at(gi))
        CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
  }
  return enabled;
}

std::string get_parameter_filename() {
  std::vector<char> buf(10000);
  int ilen = CCTK_ParameterFilename(buf.size(), buf.data());
  assert(ilen < int(buf.size() - 1));
  std::string parfilename(buf.data());
  // Remove directory prefix, if any
  auto slash = parfilename.rfind('/');
  if (slash != std::string::npos)
    parfilename = parfilename.substr(slash + 1);
  // Remove suffix, if it is there
  auto suffix = parfilename.rfind('.');
  if (suffix != std::string::npos && parfilename.substr(suffix) == ".par")
    parfilename = parfilename.substr(0, suffix);
  return parfilename;
}

std::string get_simulation_name() {
  std::string name = get_parameter_filename();
  const size_t last_slash = name.rfind('/');
  if (last_slash != std::string::npos && last_slash < name.length())
    name = name.substr(last_slash + 1);
  const size_t last_dot = name.rfind('.');
  if (last_dot != std::string::npos && last_dot > 0)
    name = name.substr(0, last_dot);
  return name;
}

} // namespace

extern "C" void PlanesX_Output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  {
    const int carpetx_openpmd_every =
        int(get_carpetx_int_param("out_openpmd_every"));
    const int openpmd_every =
        carpetx_openpmd_every == -1 ? out_every : carpetx_openpmd_every;
    const int every =
        openpmd_planes_every == -1 ? openpmd_every : openpmd_planes_every;
    if (every > 0 && cctk_iteration % every == 0 &&
        strlen(openpmd_planes) != 0) {
      const std::vector<bool> group_enabled =
          find_groups("openPMD planes", openpmd_plane_vars);
#ifdef HAVE_CAPABILITY_openPMD_api
      const std::vector<plane_spec_t> planes = parse_planes(
          openpmd_planes, planes_int_precision, planes_frac_precision);
      const std::string simulation_name = get_simulation_name();
      OutputOpenPMDPlanes(cctkGH, group_enabled, planes, out_dir,
                          simulation_name);
#else
      CCTK_VERROR("openPMD is not enabled. The parameter "
                  "PlanesX::openpmd_planes must be empty.");
#endif
    }
  }

  {
    const int carpetx_silo_every = int(get_carpetx_int_param("out_silo_every"));
    const int silo_every =
        carpetx_silo_every == -1 ? out_every : carpetx_silo_every;
    const int every = silo_planes_every == -1 ? silo_every : silo_planes_every;
    if (every > 0 && cctk_iteration % every == 0 && strlen(silo_planes) != 0) {
      const std::vector<bool> group_enabled =
          find_groups("Silo planes", silo_plane_vars);
#ifdef HAVE_CAPABILITY_Silo
      const std::vector<plane_spec_t> planes = parse_planes(
          silo_planes, planes_int_precision, planes_frac_precision);
      const std::string simulation_name = get_simulation_name();
      OutputSiloPlanes(cctkGH, group_enabled, planes, out_dir, simulation_name);
#else
      CCTK_VERROR("Silo is not enabled. The parameter "
                  "PlanesX::silo_planes must be empty.");
#endif
    }
  }
}

} // namespace PlanesX
