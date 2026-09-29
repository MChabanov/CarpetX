#ifndef CARPETX_CARPETX_FILLPATCH_HXX
#define CARPETX_CARPETX_FILLPATCH_HXX

#include "driver.hxx"
#include "task_manager.hxx"

#include <functional>

namespace CarpetX {

// Sync
void FillPatch_Sync(task_manager &tasks2,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    amrex::MultiFab &mfab, const amrex::Geometry &geom);

// Prolongate ghosts from coarse level, optionally with same-level sync.
// When do_sync=true, also performs FillBoundary (same-level ghost exchange).
// When do_sync=false, only performs coarse-to-fine interpolation.
//
// Optional time-blend: when cmfab_old != nullptr and w_new != 1, the coarse
// patch is filled as w_new*cmfab + (1-w_new)*cmfab_old before interpolation.
// This is used to time-interpolate the coarse source when the fine level is
// mid-subcycle (misaligned in time with the coarse level).
void FillPatch_Prolongate(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs,
    bool do_sync, const amrex::MultiFab *cmfab_old = nullptr,
    CCTK_REAL w_new = 1);

// Prolongate and sync ghosts (same-level exchange + coarse-to-fine
// interpolation)
inline void FillPatch_ProlongateGhosts(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs,
    const amrex::MultiFab *cmfab_old = nullptr, CCTK_REAL w_new = 1) {
  FillPatch_Prolongate(tasks2, tasks3, groupdata, coarsegroupdata, mfab, cmfab,
                       fgeom, cgeom, mapper, bcrecs, /*do_sync=*/true,
                       cmfab_old, w_new);
}

// Prolongate only (coarse-to-fine interpolation, no same-level exchange)
inline void FillPatch_ProlongateOnly(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs) {
  FillPatch_Prolongate(tasks2, tasks3, groupdata, coarsegroupdata, mfab, cmfab,
                       fgeom, cgeom, mapper, bcrecs, /*do_sync=*/false);
}

// The back half of a coarse-to-fine ghost fill, shared by
// `FillPatch_Prolongate` (temporary buffers) and the subcycling RK fill
// `FillRKBoundary` (persistent buffers). The caller owns both buffers and has
// already filled `crse_patch` (zero ghosts, FPinfo::ba_crse_patch on
// FPinfo::dm_patch) with the coarse state; `fine_patch` (zero ghosts,
// FPinfo::ba_fine_patch on the same distribution) is overwritten. Both must
// stay alive until `Prolongate_Finish` has returned.
//
// Coarse boundary conditions on `crse_patch`, spatial interpolation into
// `fine_patch` with `mapper`, then the start of the copy into the ghosts of
// `mfab`.
void Prolongate_Start(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, amrex::MultiFab &crse_patch,
    amrex::MultiFab &fine_patch, const amrex::Geometry &fgeom,
    const amrex::Geometry &cgeom, amrex::Interpolater *mapper,
    const amrex::Vector<amrex::BCRec> &bcrecs);

// Finish the copy into the ghosts of `mfab`, then apply the fine boundary
// conditions (after the prolongation, because symmetry boundary conditions
// might require prolongated points).
void Prolongate_Finish(const GHExt::PatchData::LevelData::GroupData &groupdata,
                       amrex::MultiFab &mfab);

#warning "TODO: Restrict"

// Prolongate and sync interior. Expects coarse mfab prolongated and
// synced. ("InterpFromCoarseLevel")
void FillPatch_NewLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &cgeom, const amrex::Geometry &fgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

// Vector-potential triples.
//
// The flux-preserving prolongation of an edge-staggered A_i couples the three
// components, so they are prolonged in one pass rather than group by group.
// `groups` holds the three group indices in x, y, z order, as resolved by
// get_group_vector_potential; the caller dispatches once, on the x member.
//
// Two things differ from the scalar paths and are the reason these need their
// own entry points at all.  The three coarse patches must all have their
// boundary conditions applied before any of them is read, because prolonging
// A_x reads A_y and A_z.  And the fine target must be snapped out to whole
// coarse cells, because the operator's interior closure is defined per coarse
// cell and needs all six of that cell's faces; the snapped patches may
// overlap, which is harmless since the operator is a pure function of the
// coarse data and so the overlaps agree.
//
// ("InterpFromCoarseLevel" for three coupled components.)
void FillPatch_NewLevel_vecpot(
    const GHExt::PatchData::LevelData &leveldata,
    const GHExt::PatchData::LevelData &coarseleveldata,
    const std::array<int, dim> &groups, int tl, const amrex::Geometry &cgeom,
    const amrex::Geometry &fgeom);

// Prolongate the ghosts of a vector-potential triple from the coarse level,
// optionally with the same-level sync first (as FillPatch_ProlongateGhosts
// versus FillPatch_ProlongateOnly).
//
// Unlike the scalar path this preserves the fine level's own data: an edge the
// fine grid evolves is not part of the prolongation at all.  That is arranged
// exactly as AMReX arranges it for FaceDivFree -- the destination patch is
// snapped to whole coarse cells and pre-filled from the fine level, and a
// solve mask per class of owning coarse object tells the operator which
// entries it may write.  The operator's interior closure reads the destination
// rather than the coarse data, so copied and prolonged faces need no branch.
//
// Runs synchronously rather than through the task managers.  The scalar path
// overlaps its coarse-patch copy with other groups' work; doing the same here
// means interleaving seven mask builds and three patch copies, which is worth
// doing only once the operator has been exercised in a real run.
void FillPatch_Prolongate_vecpot(
    const GHExt::PatchData::LevelData &leveldata,
    const GHExt::PatchData::LevelData &coarseleveldata,
    const std::array<int, dim> &groups, int tl, const amrex::Geometry &fgeom,
    const amrex::Geometry &cgeom, bool do_sync);

// As FillPatch_RemakeLevel, for a triple.  `fmfab` is the level's previous
// data: it fixes the coarse-fine footprint, pre-fills the destination and
// defines the solve masks, exactly as the level's own data does during a ghost
// fill, and is then copied over the result.
void FillPatch_RemakeLevel_vecpot(
    const GHExt::PatchData::LevelData &leveldata,
    const GHExt::PatchData::LevelData &coarseleveldata,
    const std::array<int, dim> &groups, int tl,
    const std::array<const amrex::MultiFab *, dim> &fmfab,
    const amrex::Geometry &cgeom, const amrex::Geometry &fgeom);

// ("FillPatchTwoLevels")
void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::MultiFab &fmfab, const amrex::Geometry &cgeom,
    const amrex::Geometry &fgeom, amrex::Interpolater *mapper,
    const amrex::Vector<amrex::BCRec> &bcrecs);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_FILLPATCH_HXX
