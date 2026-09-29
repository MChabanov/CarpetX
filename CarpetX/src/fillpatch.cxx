#include "fillpatch.hxx"
#include "prolongate_3d_rf2_vecpot.hxx"
#include "schedule.hxx"

#include <utility>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Version.H>

namespace CarpetX {

using namespace amrex;
#if AMREX_RELEASE_NUMBER >= 240500
using namespace amrex::detail;
#endif

// The code in this file is written in a "coroutine style"
// <https://en.wikipedia.org/wiki/Coroutine>. That is, each function
// returns another function that describes what to do next. This
// allows the caller to interleave many function calls, for example to
// schedule many calls to `MPI_Irecv` and `MPI_Isend` simultaneously.
//
// This programming style is obviously quite tedious. C++20 will have
// special support for this via `co_yield` etc.
// <https://en.cppreference.com/w/cpp/coroutine>, and the functions in
// this file will then look like normal functions.
//
// Coroutines were popularized in the "Modula" language in the 1980s.
// Welcome to the future, C++, you're only 40 years behind.

void FillPatch_Sync(task_manager &tasks2,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    MultiFab &mfab, const Geometry &geom) {
  assert(!groupdata.mfab.empty());
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           geom.periodicity());
  tasks2.submit_serially([&groupdata, &mfab]() {
    mfab.FillBoundary_finish();
    groupdata.apply_boundary_conditions(mfab);
  });
}

void FillPatch_Prolongate(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const Geometry &fgeom,
    const Geometry &cgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs, const bool do_sync,
    const MultiFab *const cmfab_old, const CCTK_REAL w_new) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());
  const IntVect &nghosts = mfab.nGrowVect();
  if (nghosts.max() == 0)
    return;

  // Time-blend: when an old coarse snapshot is provided and the weight is not
  // unity, the coarse patch is filled as w_new*cmfab + (1-w_new)*cmfab_old
  // before the spatial interpolation.
  const bool do_blend = cmfab_old && w_new != CCTK_REAL(1);

  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      mfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  // Same-level ghost exchange (only when do_sync)
  if (do_sync)
    mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                             fgeom.periodicity());

  if (fpc.ba_crse_patch.empty()) {
    // There is no coarser level for our boundaries, i.e. there is no
    // prolongation.

    // Early-exit cleanup (only when do_sync)
    if (do_sync) {
      tasks2.submit_serially([&groupdata, &mfab]() {
        mfab.FillBoundary_finish();
        groupdata.apply_boundary_conditions(mfab);
      });
    }
    return;
  }

  // Prolongate from the next coarser level. Apply the boundary
  // conditions after the prolongation is done (because symmetry
  // boundary conditions might require prolongated points).

  // Copy parts of coarse grid into temporary buffer
  MultiFab *const mfab_crse_patch_ptr =
      new MultiFab(make_mf_crse_patch<MultiFab>(fpc, ncomps));
  MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;
  mf_set_domain_bndry(mfab_crse_patch, cgeom);

  // This is not local
  mfab_crse_patch.ParallelCopy_nowait(
      cmfab, 0, 0, ncomps, IntVect{0} /* don't use coarse ghosts */,
      mfab_crse_patch.nGrowVect(), cgeom.periodicity());

  // Optionally copy the old coarse snapshot into a sibling buffer, to be
  // blended with the new one before interpolation. Same lifetime/delete
  // pattern as mfab_crse_patch; its copy is overlapped with the one above.
  MultiFab *mfab_crse_patch_old_ptr = nullptr;
  if (do_blend) {
    mfab_crse_patch_old_ptr =
        new MultiFab(make_mf_crse_patch<MultiFab>(fpc, ncomps));
    mf_set_domain_bndry(*mfab_crse_patch_old_ptr, cgeom);
    mfab_crse_patch_old_ptr->ParallelCopy_nowait(
        *cmfab_old, 0, 0, ncomps, IntVect{0} /* don't use coarse ghosts */,
        mfab_crse_patch_old_ptr->nGrowVect(), cgeom.periodicity());
  }

  tasks2.submit_serially([&tasks3, &groupdata, &coarsegroupdata, &mfab, &cgeom,
                          &fgeom, mapper, &bcrecs, &fpc, mfab_crse_patch_ptr,
                          mfab_crse_patch_old_ptr, w_new, do_sync]() {
    const int ncomps = mfab.nComp();
    MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;

    // Finish same-level sync (only when do_sync)
    if (do_sync)
      mfab.FillBoundary_finish();

    // Finish copying parts of coarse grid into temporary buffer
    mfab_crse_patch.ParallelCopy_finish();

    // Blend the new and old coarse snapshots in place, in physical time, before
    // applying boundary conditions: applying the (linear) coarse BC after the
    // blend is equivalent to blending the BCs and keeps a single code path.
    if (mfab_crse_patch_old_ptr) {
      MultiFab &mfab_crse_patch_old = *mfab_crse_patch_old_ptr;
      mfab_crse_patch_old.ParallelCopy_finish();
      const CCTK_REAL w_old = CCTK_REAL(1) - w_new;
      // dst = w_new*new + w_old*old  (elementwise, in place)
      MultiFab::LinComb(mfab_crse_patch, w_new, mfab_crse_patch, 0, w_old,
                        mfab_crse_patch_old, 0, 0, ncomps, IntVect{0});
      delete mfab_crse_patch_old_ptr;
    }

    MultiFab *const mfab_fine_patch_ptr =
        new MultiFab(make_mf_fine_patch<MultiFab>(fpc, ncomps));
    MultiFab &mfab_fine_patch = *mfab_fine_patch_ptr;

    // Coarse boundary conditions, interpolation into the fine buffer, and
    // start of the copy into the destination
    Prolongate_Start(groupdata, coarsegroupdata, mfab, mfab_crse_patch,
                     mfab_fine_patch, fgeom, cgeom, mapper, bcrecs);

    delete mfab_crse_patch_ptr;

    tasks3.submit_serially([&groupdata, &mfab, mfab_fine_patch_ptr]() {
      // Finish the copy into the destination, fine boundary conditions
      Prolongate_Finish(groupdata, mfab);

      delete mfab_fine_patch_ptr;
    });
  });
}

void Prolongate_Start(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, MultiFab &crse_patch, MultiFab &fine_patch,
    const Geometry &fgeom, const Geometry &cgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());
  const IntVect &nghosts = mfab.nGrowVect();
  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  assert(crse_patch.nComp() == ncomps);
  assert(fine_patch.nComp() == ncomps);
  // `FillPatchInterp` walks the fine buffer and indexes the coarse buffer box
  // by box
  assert(crse_patch.DistributionMap() == fine_patch.DistributionMap());
  assert(crse_patch.size() == fine_patch.size());

  coarsegroupdata.apply_boundary_conditions(crse_patch);

  // Interpolate coarse buffer into fine buffer (in space, local)
  FillPatchInterp(fine_patch, 0, crse_patch, 0, ncomps,
                  IntVect{0} /* don't add any new ghosts */, cgeom, fgeom,
                  grow(convert(fgeom.Domain(), mfab.ixType()), nghosts), ratio,
                  mapper, bcrecs, 0);

  // Copy fine buffer into destination
  mfab.ParallelCopy_nowait(
      fine_patch, 0, 0, ncomps,
      IntVect{0} /* don't use any ghosts from the buffer */, nghosts);
}

void Prolongate_Finish(const GHExt::PatchData::LevelData::GroupData &groupdata,
                       MultiFab &mfab) {
  // Finish copying fine buffer into destination
  mfab.ParallelCopy_finish();

  // Apply symmetry and boundary conditions
  groupdata.apply_boundary_conditions(mfab);
}

void FillPatch_NewLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const Geometry &cgeom,
    const Geometry &fgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());
  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect &nghosts = mfab.nGrowVect();
  // const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const BoxArray &ba = mfab.boxArray();
  const DistributionMapping &dm = mfab.DistributionMap();

  const IndexType &ixtype = ba.ixType();
  assert(ixtype == cmfab.boxArray().ixType());

  // Suffix `_g` is for "with ghosts added"
  Box fdomain_g(amrex::convert(fgeom.Domain(), mfab.ixType()));
  for (int d = 0; d < dim; ++d)
    if (fgeom.isPeriodic(d))
      fdomain_g.grow(d, nghosts[d]);

  const int nboxes = ba.size();
  BoxArray cba_g(nboxes);
  for (int i = 0; i < nboxes; ++i) {
    Box box = amrex::convert(amrex::grow(ba[i], nghosts), ixtype);
    box &= fdomain_g;
    cba_g.set(i, coarsener.doit(box));
  }
  MultiFab cmfab_g(cba_g, dm, ncomps, 0);
  mf_set_domain_bndry(cmfab_g, cgeom);

  cmfab_g.ParallelCopy(cmfab, 0, 0, ncomps, cgeom.periodicity());

  coarsegroupdata.apply_boundary_conditions(cmfab_g);

  FillPatchInterp(mfab, 0, cmfab_g, 0, ncomps, nghosts, cgeom, fgeom, fdomain_g,
                  ratio, mapper, bcrecs, 0);

  groupdata.apply_boundary_conditions(mfab);
}

void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const MultiFab &fmfab,
    const Geometry &cgeom, const Geometry &fgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());
  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect &nghosts = mfab.nGrowVect();
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      fmfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  if (!fpc.ba_crse_patch.empty()) {
    MultiFab mfab_crse_patch = make_mf_crse_patch<MultiFab>(fpc, ncomps);
    mf_set_domain_bndry(mfab_crse_patch, cgeom);

    mfab_crse_patch.ParallelCopy(
        cmfab, 0, 0, ncomps, IntVect{0} /* don't use coarse ghosts */,
        mfab_crse_patch.nGrowVect(), cgeom.periodicity());
    coarsegroupdata.apply_boundary_conditions(mfab_crse_patch);

    MultiFab mfab_fine_patch = make_mf_fine_patch<MultiFab>(fpc, ncomps);

    // In space, local
    FillPatchInterp(mfab_fine_patch, 0, mfab_crse_patch, 0, ncomps,
                    IntVect{0} /* don't add any new ghosts */, cgeom, fgeom,
                    grow(convert(fgeom.Domain(), mfab.ixType()), nghosts),
                    ratio, mapper, bcrecs, 0);

    mfab.ParallelCopy_nowait(
        mfab_fine_patch, 0, 0, ncomps,
        IntVect{0} /* don't use any ghosts from the buffer */, nghosts);
    mfab.ParallelCopy_finish();
  }

  mfab.ParallelCopy(fmfab, 0, 0, ncomps, IntVect{0} /* don't use old ghosts */,
                    nghosts, fgeom.periodicity());
  groupdata.apply_boundary_conditions(mfab);
}


////////////////////////////////////////////////////////////////////////////////
// Vector-potential triples.  See fillpatch.hxx for the contract.

namespace {

// The three components of a triple, gathered from a level.
struct vecpot_level {
  std::array<GHExt::PatchData::LevelData::GroupData *, dim> gd;
  std::array<MultiFab *, dim> mfab;
};

vecpot_level gather(const GHExt::PatchData::LevelData &leveldata,
                    const std::array<int, dim> &groups, const int tl) {
  vecpot_level r;
  for (int d = 0; d < dim; ++d) {
    r.gd[d] = leveldata.groupdata.at(groups[d]).get();
    assert(!r.gd[d]->mfab.empty());
    r.mfab[d] = r.gd[d]->mfab.at(tl).get();
  }
  // The three MultiFabs share a decomposition and differ only in index type;
  // that is what makes one coarse patch geometry usable for all of them.
  for (int d = 1; d < dim; ++d) {
    assert(r.mfab[d]->nComp() == r.mfab[0]->nComp());
    assert(r.mfab[d]->nGrowVect() == r.mfab[0]->nGrowVect());
    assert(r.mfab[d]->size() == r.mfab[0]->size());
  }
  return r;
}

// Index type of the coarse object that owns a class of fine edges
// (PROLONG_A.md section 5.2).  AMReX convention here: 1 == NODE.
IndexType vecpot_edge_ixtype(const int d) {
  return IndexType(IntVect::TheNodeVector() - IntVect::TheDimensionVector(d));
}
IndexType vecpot_face_ixtype(const int f) {
  return IndexType(IntVect::TheDimensionVector(f));
}

// Shared core of the coarse-to-fine fill for a triple.
//
// `dst` is written.  `fsrc` is the fine data that already exists and must be
// preserved: it fixes the coarse-fine footprint, pre-fills the destination
// patch, and defines the solve masks.  For a ghost fill `fsrc` is `dst`
// itself; for a regrid it is the level's previous MultiFabs, and it is empty
// on a level that has just been created.
//
// The two structural requirements of the operator are met here.  All three
// coarse patches have their boundary conditions applied before any of them is
// read, because prolonging A_x reads A_y and A_z.  And the destination patch
// is snapped out to whole coarse cells, because the interior closure is
// defined per coarse cell and needs all six of that cell's faces; the snapped
// patches may overlap, which is harmless since the operator is a pure function
// of its inputs and so the overlaps agree.
void vecpot_fill_from_coarse(
    const vecpot_level &dst,
    const std::array<const MultiFab *, dim> &fsrc,
    const std::array<const GHExt::PatchData::LevelData::GroupData *, dim> &cgd,
    const std::array<const MultiFab *, dim> &cmfab, const Geometry &fgeom,
    const Geometry &cgeom) {

  const int ncomps = dst.mfab[0]->nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect nghosts = dst.mfab[0]->nGrowVect();

  const InterpolaterBoxCoarsener &coarsener =
      prolongate_vecpot_3d_rf2.BoxCoarsener(ratio);

  // One index-type independent decomposition for all three components: the
  // recipe AMReX uses in its own Array-valued FillPatch.
  const BoxArray dba_cc =
      amrex::convert(dst.mfab[0]->boxArray(), IntVect::TheZeroVector());
  const BoxArray sba_cc =
      amrex::convert(fsrc[0]->boxArray(), IntVect::TheZeroVector());
  const MultiFab dst_cc_dummy(dba_cc, dst.mfab[0]->DistributionMap(), ncomps,
                              nghosts, MFInfo().SetAlloc(false));
  const MultiFab src_cc_dummy(sba_cc, fsrc[0]->DistributionMap(), ncomps,
                              nghosts, MFInfo().SetAlloc(false));
  const EB2::IndexSpace *const index_space = nullptr;
  const FabArrayBase::FPinfo &fpc =
      FabArrayBase::TheFPinfo(src_cc_dummy, dst_cc_dummy, nghosts, coarsener,
                              fgeom, cgeom, index_space);

  if (fpc.ba_crse_patch.empty())
    return; // nothing needs the coarser level

  // Cell-centred footprint of the destination patch, snapped to whole coarse
  // cells; see above.
  const BoxArray rba_cc =
      amrex::refine(amrex::coarsen(fpc.ba_fine_patch, ratio), ratio);

  // A mask entry is 1 where the operator must solve and 0 where the fine level
  // already owns the value.  "Known" is the existing fine data's box array,
  // coarsened and converted to the owning object's index type, so that an
  // object lying on the boundary of the fine region counts as known -- which
  // is what the copy-from-fine rule requires, nodal data belonging to the
  // boxes it bounds.
  const auto build_mask = [&](const IndexType &ixt) {
    iMultiFab mask = make_mf_crse_mask<iMultiFab>(fpc, ncomps, ixt, ratio);
    const MultiFab known(amrex::convert(amrex::coarsen(sba_cc, ratio), ixt),
                         fsrc[0]->DistributionMap(), ncomps, nghosts,
                         MFInfo().SetAlloc(false));
    const MultiFab solution(mask.boxArray(), mask.DistributionMap(), ncomps, 0,
                            MFInfo().SetAlloc(false));
    const FabArrayBase::CPC cpc(solution, IntVect::TheZeroVector(), known,
                                IntVect::TheZeroVector(), cgeom.periodicity());
    mask.setVal(1);                 // to solve
    mask.setVal(0, cpc, 0, ncomps); // already known
    return mask;
  };

  std::array<MultiFab, dim> cpatch, rpatch;
  std::array<iMultiFab, dim> medge, mface;
  for (int d = 0; d < dim; ++d) {
    const IndexType ixt = dst.mfab[d]->boxArray().ixType();

    cpatch[d] = make_mf_crse_patch<MultiFab>(fpc, ncomps, ixt);
    mf_set_domain_bndry(cpatch[d], cgeom);
    cpatch[d].ParallelCopy(*cmfab[d], 0, 0, ncomps,
                           IntVect{0} /* don't use coarse ghosts */,
                           cpatch[d].nGrowVect(), cgeom.periodicity());

    rpatch[d].define(amrex::convert(rba_cc, ixt), fpc.dm_patch, ncomps, 0);
    mf_set_domain_bndry(rpatch[d], fgeom);
    // Pre-fill with the fine level's own values: those edges are evolved by
    // the fine grid and must not be prolonged over.
    rpatch[d].ParallelCopy(*fsrc[d], 0, 0, ncomps, nghosts, IntVect{0},
                           fgeom.periodicity());

    medge[d] = build_mask(vecpot_edge_ixtype(d));
    mface[d] = build_mask(vecpot_face_ixtype(d));
  }
  iMultiFab mcell = build_mask(IndexType::TheCellType());

  for (int d = 0; d < dim; ++d)
    cgd[d]->apply_boundary_conditions(cpatch[d]);

  const GpuArray<CCTK_REAL, dim> coarse_dx = cgeom.CellSizeArray();

  for (MFIter mfi(rpatch[0]); mfi.isValid(); ++mfi) {
    vecpot_masks masks;
    for (int d = 0; d < dim; ++d) {
      masks.edge[d] = &medge[d][mfi];
      masks.face[d] = &mface[d][mfi];
    }
    masks.cell = &mcell[mfi];
    const Array<const FArrayBox *, dim> cfab{&cpatch[0][mfi], &cpatch[1][mfi],
                                             &cpatch[2][mfi]};
    const Array<FArrayBox *, dim> ffab{&rpatch[0][mfi], &rpatch[1][mfi],
                                       &rpatch[2][mfi]};
    prolongate_vecpot_3d_rf2.interp_vecpot(cfab, 0, ffab, 0, ncomps,
                                           rba_cc[mfi.index()], ratio, masks,
                                           coarse_dx, RunOn::Gpu);
  }

  for (int d = 0; d < dim; ++d)
    // ToGhost, not ParallelCopy: the snapped patch overlaps the destination's
    // valid region, which the fine level owns.
    dst.mfab[d]->ParallelCopyToGhost(rpatch[d], 0, 0, ncomps, IntVect{0},
                                     nghosts);
}

} // namespace

void FillPatch_NewLevel_vecpot(
    const GHExt::PatchData::LevelData &leveldata,
    const GHExt::PatchData::LevelData &coarseleveldata,
    const std::array<int, dim> &groups, const int tl, const Geometry &cgeom,
    const Geometry &fgeom) {

  const vecpot_level fl = gather(leveldata, groups, tl);
  std::array<const GHExt::PatchData::LevelData::GroupData *, dim> cgd;
  std::array<const MultiFab *, dim> cmfab;
  for (int d = 0; d < dim; ++d) {
    cgd[d] = coarseleveldata.groupdata.at(groups[d]).get();
    cmfab[d] = cgd[d]->mfab.at(tl).get();
  }

  const int ncomps = fl.mfab[0]->nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect nghosts = fl.mfab[0]->nGrowVect();
  const DistributionMapping &dm = fl.mfab[0]->DistributionMap();

  const InterpolaterBoxCoarsener &coarsener =
      prolongate_vecpot_3d_rf2.BoxCoarsener(ratio);

  // A new level has no fine data at all, so there is no coarse-fine footprint
  // to compute and no mask: every point is prolonged.  That makes this the
  // simplest of the paths, and the natural one to bring up first.
  const BoxArray ba_cc =
      amrex::convert(fl.mfab[0]->boxArray(), IntVect::TheZeroVector());

  Box fdomain_g_cc = amrex::convert(fgeom.Domain(), IntVect::TheZeroVector());
  for (int d = 0; d < dim; ++d)
    if (fgeom.isPeriodic(d))
      fdomain_g_cc.grow(d, nghosts[d]);

  const int nboxes = ba_cc.size();
  BoxArray fba_snap(nboxes), cba_g(nboxes);
  for (int i = 0; i < nboxes; ++i) {
    Box t = amrex::grow(ba_cc[i], nghosts);
    t &= fdomain_g_cc;
    fba_snap.set(i, amrex::refine(amrex::coarsen(t, ratio), ratio));
    cba_g.set(i, coarsener.doit(t));
  }

  std::array<MultiFab, dim> cpatch, fpatch;
  for (int d = 0; d < dim; ++d) {
    const IndexType ixt = fl.mfab[d]->boxArray().ixType();
    cpatch[d].define(amrex::convert(cba_g, ixt), dm, ncomps, 0);
    mf_set_domain_bndry(cpatch[d], cgeom);
    cpatch[d].ParallelCopy(*cmfab[d], 0, 0, ncomps, cgeom.periodicity());
    fpatch[d].define(amrex::convert(fba_snap, ixt), dm, ncomps, 0);
  }

  for (int d = 0; d < dim; ++d)
    cgd[d]->apply_boundary_conditions(cpatch[d]);

  const GpuArray<CCTK_REAL, dim> coarse_dx = cgeom.CellSizeArray();
  const vecpot_masks masks; // nothing to preserve

  for (MFIter mfi(fpatch[0]); mfi.isValid(); ++mfi) {
    const Array<const FArrayBox *, dim> cfab{&cpatch[0][mfi], &cpatch[1][mfi],
                                             &cpatch[2][mfi]};
    const Array<FArrayBox *, dim> ffab{&fpatch[0][mfi], &fpatch[1][mfi],
                                       &fpatch[2][mfi]};
    prolongate_vecpot_3d_rf2.interp_vecpot(cfab, 0, ffab, 0, ncomps,
                                           fba_snap[mfi.index()], ratio, masks,
                                           coarse_dx, RunOn::Gpu);
  }

  for (int d = 0; d < dim; ++d) {
    fl.mfab[d]->ParallelCopy(fpatch[d], 0, 0, ncomps,
                             IntVect{0} /* no ghosts from the buffer */,
                             nghosts);
    fl.gd[d]->apply_boundary_conditions(*fl.mfab[d]);
  }
}

void FillPatch_Prolongate_vecpot(
    const GHExt::PatchData::LevelData &leveldata,
    const GHExt::PatchData::LevelData &coarseleveldata,
    const std::array<int, dim> &groups, const int tl, const Geometry &fgeom,
    const Geometry &cgeom, const bool do_sync) {

  const vecpot_level fl = gather(leveldata, groups, tl);
  if (fl.mfab[0]->nGrowVect().max() == 0)
    return;

  std::array<const GHExt::PatchData::LevelData::GroupData *, dim> cgd;
  std::array<const MultiFab *, dim> cmfab;
  std::array<const MultiFab *, dim> fsrc;
  for (int d = 0; d < dim; ++d) {
    cgd[d] = coarseleveldata.groupdata.at(groups[d]).get();
    cmfab[d] = cgd[d]->mfab.at(tl).get();
    fsrc[d] = fl.mfab[d]; // the fine data to preserve is the level's own
  }

  // Same-level exchange first, so that the pre-fill inside sees the fine
  // level's ghost data as well as its interior.
  if (do_sync)
    for (int d = 0; d < dim; ++d)
      fl.mfab[d]->FillBoundary(fgeom.periodicity());

  vecpot_fill_from_coarse(fl, fsrc, cgd, cmfab, fgeom, cgeom);

  for (int d = 0; d < dim; ++d)
    fl.gd[d]->apply_boundary_conditions(*fl.mfab[d]);
}

void FillPatch_RemakeLevel_vecpot(
    const GHExt::PatchData::LevelData &leveldata,
    const GHExt::PatchData::LevelData &coarseleveldata,
    const std::array<int, dim> &groups, const int tl,
    const std::array<const amrex::MultiFab *, dim> &fmfab,
    const Geometry &cgeom, const Geometry &fgeom) {

  const vecpot_level fl = gather(leveldata, groups, tl);
  std::array<const GHExt::PatchData::LevelData::GroupData *, dim> cgd;
  std::array<const MultiFab *, dim> cmfab;
  for (int d = 0; d < dim; ++d) {
    cgd[d] = coarseleveldata.groupdata.at(groups[d]).get();
    cmfab[d] = cgd[d]->mfab.at(tl).get();
  }

  const int ncomps = fl.mfab[0]->nComp();
  const IntVect nghosts = fl.mfab[0]->nGrowVect();

  // The level's previous data is what must be preserved: it fixes the
  // coarse-fine footprint, pre-fills the destination patch and defines the
  // masks, exactly as the level's own data does during a ghost fill.
  vecpot_fill_from_coarse(fl, fmfab, cgd, cmfab, fgeom, cgeom);

  for (int d = 0; d < dim; ++d) {
    fl.mfab[d]->ParallelCopy(*fmfab[d], 0, 0, ncomps,
                             IntVect{0} /* don't use old ghosts */, nghosts,
                             fgeom.periodicity());
    fl.gd[d]->apply_boundary_conditions(*fl.mfab[d]);
  }
}

} // namespace CarpetX
