#include "prolongate_3d_rf2_vecpot.hxx"
#include "prolongate_3d_rf2_vecpot_kernels.hxx"

#include <AMReX_GpuLaunch.H>

#include <cassert>

namespace CarpetX {
using namespace amrex;
using namespace CarpetX::vecpot;

namespace {

// Index type of the coarse object that owns a class of fine edges
// (PROLONG_A.md section 5.2).  CarpetX's `indextype` convention is the
// opposite of AMReX's, so these are written in AMReX terms: 1 == NODE.

// A_d edges: CELL in d, NODE in the other two.
inline IndexType edge_ixtype(const int d) {
  return IndexType(IntVect::TheNodeVector() - IntVect::TheDimensionVector(d));
}
// f-faces: NODE in f, CELL in the other two.
inline IndexType face_ixtype(const int f) {
  return IndexType(IntVect::TheDimensionVector(f));
}

inline Array4<const int> mask_array(const IArrayBox *const p) {
  return p ? p->const_array() : Array4<const int>{};
}

} // namespace

prolongate_3d_rf2_vecpot prolongate_vecpot_3d_rf2;

prolongate_3d_rf2_vecpot::~prolongate_3d_rf2_vecpot() {}

Box prolongate_3d_rf2_vecpot::CoarseBox(const Box &fine, const int ratio) {
  return CoarseBox(fine, IntVect(ratio));
}

Box prolongate_3d_rf2_vecpot::CoarseBox(const Box &fine, const IntVect &ratio) {
  for (int d = 0; d < dim; ++d)
    assert(ratio[d] == 2);
  // Work on a cell-centred box that certainly contains `fine`.  The extra cell
  // on the high side absorbs the nodal-to-cell conversion: a nodal target
  // index reaches one coarse index past the coarsened cell box.
  Box cc = enclosedCells(surroundingNodes(fine));
  for (int d = 0; d < dim; ++d)
    cc.growHi(d, 1);
  Box crse = amrex::coarsen(cc, ratio);
  // Stencil radius: Stage 2 reads coarse B one coarse cell out, and B is a
  // curl of A, so A is read two coarse cells out.  See PROLONG_A.md 5.3 for
  // why 2 and not 3.
  crse.grow(2);
  return crse;
}

void prolongate_3d_rf2_vecpot::interp(
    const FArrayBox & /*crse*/, int /*crse_comp*/, FArrayBox & /*fine*/,
    int /*fine_comp*/, int /*ncomp*/, const Box & /*fine_region*/,
    const IntVect & /*ratio*/, const Geometry & /*crse_geom*/,
    const Geometry & /*fine_geom*/, Vector<BCRec> const & /*bcr*/,
    int /*actual_comp*/, int /*actual_state*/, RunOn /*gpu_or_cpu*/) {
  CCTK_ERROR("prolongate_3d_rf2_vecpot couples the three components of the "
             "vector potential and cannot act on a single FArrayBox. Call "
             "interp_vecpot instead.");
}

void prolongate_3d_rf2_vecpot::interp_arr(
    Array<FArrayBox *, dim> const & /*crse*/, int /*crse_comp*/,
    Array<FArrayBox *, dim> const & /*fine*/, int /*fine_comp*/,
    int /*ncomp*/, const Box & /*fine_region*/, const IntVect & /*ratio*/,
    Array<IArrayBox *, dim> const & /*solve_mask*/,
    const Geometry & /*crse_geom*/, const Geometry & /*fine_geom*/,
    Vector<Array<BCRec, dim> > const & /*bcr*/, int /*actual_comp*/,
    int /*actual_state*/, RunOn /*runon*/) {
  CCTK_ERROR("prolongate_3d_rf2_vecpot needs seven solve masks (one per class "
             "of owning coarse object), which AMReX's interp_arr signature "
             "cannot carry. Call interp_vecpot instead.");
}

void prolongate_3d_rf2_vecpot::interp_vecpot(
    Array<const FArrayBox *, dim> const &crse, const int crse_comp,
    Array<FArrayBox *, dim> const &fine, const int fine_comp, const int ncomp,
    const Box &fine_region_cc, const IntVect &ratio, const vecpot_masks &masks,
    const GpuArray<CCTK_REAL, dim> &coarse_dx, const RunOn runon) const {

  for (int d = 0; d < dim; ++d)
    assert(ratio[d] == 2);
  assert(fine_region_cc.ixType().cellCentered());

  // One thread per coarse object; the target is expressed in coarse cells.
  const Box c_region = amrex::coarsen(fine_region_cc, ratio);

  const GpuArray<Array4<const CCTK_REAL>, dim> ac{
      crse[0]->const_array(crse_comp), crse[1]->const_array(crse_comp),
      crse[2]->const_array(crse_comp)};
  const GpuArray<Array4<CCTK_REAL>, dim> af{fine[0]->array(fine_comp),
                                            fine[1]->array(fine_comp),
                                            fine[2]->array(fine_comp)};
  const GpuArray<Array4<const int>, dim> me{mask_array(masks.edge[0]),
                                            mask_array(masks.edge[1]),
                                            mask_array(masks.edge[2])};
  const GpuArray<Array4<const int>, dim> mf{mask_array(masks.face[0]),
                                            mask_array(masks.face[1]),
                                            mask_array(masks.face[2])};
  const Array4<const int> mc = mask_array(masks.cell);

  const GpuArray<CCTK_REAL, dim> &hc = coarse_dx;

  // ---- Launch 1: fine edges lying ON a coarse edge (Stage 1) --------------
  AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(
      runon, amrex::convert(c_region, edge_ixtype(0)), bx0,
      {
        AMREX_LOOP_3D(bx0, i, j, k, {
          for (int n = 0; n < ncomp; ++n)
            surface_interp_n<0>(ac, af, me[0], idx3{{i, j, k}}, n, hc);
        });
      },
      amrex::convert(c_region, edge_ixtype(1)), bx1,
      {
        AMREX_LOOP_3D(bx1, i, j, k, {
          for (int n = 0; n < ncomp; ++n)
            surface_interp_n<4>(ac, af, me[1], idx3{{i, j, k}}, n, hc);
        });
      },
      amrex::convert(c_region, edge_ixtype(2)), bx2,
      {
        AMREX_LOOP_3D(bx2, i, j, k, {
          for (int n = 0; n < ncomp; ++n)
            surface_interp_n<8>(ac, af, me[2], idx3{{i, j, k}}, n, hc);
        });
      });

  // ---- Launch 2: fine edges in the interior of a coarse face (Stage 2) ----
  //
  // Face f owns the mid-line edges of two components: A_{(f+2)%3} with
  // parities (0,1) and A_{(f+1)%3} with parities (1,0).  For the z-face that
  // is A_y and A_x, matching the draft's four mid-line edges per face.
  AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(
      runon, amrex::convert(c_region, face_ixtype(0)), fx,
      {
        AMREX_LOOP_3D(fx, i, j, k, {
          for (int n = 0; n < ncomp; ++n) {
            surface_interp_n<9>(ac, af, mf[0], idx3{{i, j, k}}, n, hc);
            surface_interp_n<6>(ac, af, mf[0], idx3{{i, j, k}}, n, hc);
          }
        });
      },
      amrex::convert(c_region, face_ixtype(1)), fy,
      {
        AMREX_LOOP_3D(fy, i, j, k, {
          for (int n = 0; n < ncomp; ++n) {
            surface_interp_n<1>(ac, af, mf[1], idx3{{i, j, k}}, n, hc);
            surface_interp_n<10>(ac, af, mf[1], idx3{{i, j, k}}, n, hc);
          }
        });
      },
      amrex::convert(c_region, face_ixtype(2)), fz,
      {
        AMREX_LOOP_3D(fz, i, j, k, {
          for (int n = 0; n < ncomp; ++n) {
            surface_interp_n<5>(ac, af, mf[2], idx3{{i, j, k}}, n, hc);
            surface_interp_n<2>(ac, af, mf[2], idx3{{i, j, k}}, n, hc);
          }
        });
      });

  // ---- Launch 3: fine edges in the interior of a coarse cell (Stage 3) ----
  //
  // Reads only the destination array, so it is correct whether the surface
  // edges above were prolonged or copied from an already refined neighbour.
  // Must follow launches 1 and 2: it reads their output at the neighbouring
  // coarse index.
  AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, c_region, ncomp, i, j, k, n, {
    interior_interp<0>(af, mc, idx3{{i, j, k}}, n, hc);
    interior_interp<1>(af, mc, idx3{{i, j, k}}, n, hc);
    interior_interp<2>(af, mc, idx3{{i, j, k}}, n, hc);
  });
}

} // namespace CarpetX
