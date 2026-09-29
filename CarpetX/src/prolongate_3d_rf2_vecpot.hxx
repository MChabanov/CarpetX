#ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_VECPOT_HXX
#define CARPETX_CARPETX_PROLONGATE_3D_RF2_VECPOT_HXX

// The standalone unit test in CarpetX/test/unit/ compiles this header and
// prolongate_3d_rf2_vecpot.cxx against a small AMReX/Cactus shim, so that the
// operator can be exercised without an Einstein Toolkit build.  That is the
// only purpose of this guard; CarpetX itself never defines the macro.
#ifndef CARPETX_VECPOT_STANDALONE_TEST
#include "driver.hxx"
#endif

#include <AMReX_Interpolater.H>

namespace CarpetX {

// Flux-preserving prolongation for an edge-staggered vector potential.
//
// Operates on all three components at once, so it is driven through
// `interp_vecpot` rather than through `amrex::Interpolater::interp`.  See
// PROLONG_A.md for the construction, the verification status, and the
// comparison with `amrex::FaceDivFree` (the face-centred B analogue).
//
// Deviation from AMReX's `interp_arr` convention, and the reason for it:
// `FaceDivFree` needs exactly one mask per direction, because a coarse face is
// its only class of "known" object.  Edge data has four classes -- coarse
// edges, the two coarse face types that carry mid-line edges, and coarse cells
// -- so seven masks are needed (PROLONG_A.md section 8.2).  AMReX's
// `interp_arr` signature carries only three, and nothing inside AMReX can call
// us anyway (its Array FillPatch path asserts `ixType().nodeCentered(d)`,
// which edge data fails).  CarpetX therefore owns the call site and uses the
// wider entry point below; the virtual overrides exist only to fail loudly.

struct vecpot_masks {
  // Null entries mean "solve everywhere".  Indexed on the owning object's
  // index type: edge[d] on the A_d edge type, face[f] on the f-face type,
  // cell on cell centres.
  const amrex::IArrayBox *edge[dim] = {nullptr, nullptr, nullptr};
  const amrex::IArrayBox *face[dim] = {nullptr, nullptr, nullptr};
  const amrex::IArrayBox *cell = nullptr;
};

class prolongate_3d_rf2_vecpot final : public amrex::Interpolater {
public:
  virtual ~prolongate_3d_rf2_vecpot() override;

  // NOTE: takes and returns CELL-CENTRED boxes, matching the interp_arr
  // convention rather than the per-group index type used by the scalar
  // operators.  Grows by 2 coarse cells (PROLONG_A.md section 5.3).
  virtual amrex::Box CoarseBox(const amrex::Box &fine, int ratio) override;
  virtual amrex::Box CoarseBox(const amrex::Box &fine,
                               const amrex::IntVect &ratio) override;

  // Pure virtual in the base, but meaningless here: this operator couples the
  // three components.  Fails loudly.
  virtual void interp(const amrex::FArrayBox &crse, int crse_comp,
                      amrex::FArrayBox &fine, int fine_comp, int ncomp,
                      const amrex::Box &fine_region,
                      const amrex::IntVect &ratio,
                      const amrex::Geometry &crse_geom,
                      const amrex::Geometry &fine_geom,
                      amrex::Vector<amrex::BCRec> const &bcr, int actual_comp,
                      int actual_state, amrex::RunOn gpu_or_cpu) override;

  // Likewise: AMReX's three-mask signature cannot express what this operator
  // needs.  Use `interp_vecpot`.
  virtual void
  interp_arr(amrex::Array<amrex::FArrayBox *, dim> const &crse, int crse_comp,
             amrex::Array<amrex::FArrayBox *, dim> const &fine, int fine_comp,
             int ncomp, const amrex::Box &fine_region,
             const amrex::IntVect &ratio,
             amrex::Array<amrex::IArrayBox *, dim> const &solve_mask,
             const amrex::Geometry &crse_geom, const amrex::Geometry &fine_geom,
             amrex::Vector<amrex::Array<amrex::BCRec, dim> > const &bcr,
             int actual_comp, int actual_state, amrex::RunOn runon) override;

  // The real entry point.
  //
  // `crse` and `fine` hold A_x, A_y, A_z in that order.  `fine_region_cc` is
  // the target region as a CELL-CENTRED box (as AMReX passes it to
  // `interp_arr`); each component's own target is derived from it.  `fine`
  // must already hold the fine level's own values wherever `masks` says
  // "known", and must be snapped out to whole coarse cells.
  //
  // `coarse_dx` is the coarse level's cell size, i.e. `geom.CellSizeArray()`.
  // The cell size is all the operator needs from the geometry; taking it
  // directly keeps the operator a pure function of its inputs.
  void interp_vecpot(amrex::Array<const amrex::FArrayBox *, dim> const &crse,
                     int crse_comp,
                     amrex::Array<amrex::FArrayBox *, dim> const &fine,
                     int fine_comp, int ncomp,
                     const amrex::Box &fine_region_cc,
                     const amrex::IntVect &ratio, const vecpot_masks &masks,
                     const amrex::GpuArray<CCTK_REAL, dim> &coarse_dx,
                     amrex::RunOn runon) const;
};

extern prolongate_3d_rf2_vecpot prolongate_vecpot_3d_rf2;

// Self-test, run once from CarpetX_Startup.  Mirrors
// CarpetX/test/unit/test_vecpot_kernels.cxx, which is the same check without
// the Cactus dependency.
void test_prolongate_3d_rf2_vecpot();

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_VECPOT_HXX
