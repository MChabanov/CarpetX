// Standalone test of the real interpolater, CarpetX/src/prolongate_3d_rf2_vecpot.cxx,
// compiled against the AMReX/Cactus shim in shim/.  Where
// test_vecpot_kernels.cxx checks the arithmetic, this checks the glue: the
// owning-type iteration boxes, the index-type conversions, CoarseBox, the
// three launches and their ordering, and the mask handling.
//
// See CarpetX/test/unit/README.md.

#define CARPETX_VECPOT_STANDALONE_TEST 1
#include "shim/amrex_shim.hxx"

#include "../../src/prolongate_3d_rf2_vecpot.cxx"

#include <cmath>
#include <random>

using namespace amrex;
using CarpetX::vecpot_masks;
using CarpetX::prolongate_vecpot_3d_rf2;

static int failures = 0;
static void check(const char *name, double err, double tol) {
  const bool ok = std::isfinite(err) && err <= tol;
  std::printf("  %-58s %10.3e  %s\n", name, err, ok ? "ok" : "FAIL");
  if (!ok)
    ++failures;
}

static constexpr int NC = 8; // coarse cells in the target region
static constexpr int CG = 4; // coarse ghosts of the allocation

static IndexType edge_t(int d) {
  return IndexType(IntVect::TheNodeVector() - IntVect::TheDimensionVector(d));
}

// Verification curls, written out longhand and independently of the kernels.
template <typename F>
static double Bfine(int C, const F &fx, const F &fy, const F &fz, int i, int j,
                    int k, const double h[3]) {
  if (C == 0)
    return (fy(i, j, k, 0) - fy(i, j, k + 1, 0)) / h[2] +
           (fz(i, j + 1, k, 0) - fz(i, j, k, 0)) / h[1];
  if (C == 1)
    return (fz(i, j, k, 0) - fz(i + 1, j, k, 0)) / h[0] +
           (fx(i, j, k + 1, 0) - fx(i, j, k, 0)) / h[2];
  return (fx(i, j, k, 0) - fx(i, j + 1, k, 0)) / h[1] +
         (fy(i + 1, j, k, 0) - fy(i, j, k, 0)) / h[0];
}

int main() {
  const double hx = 0.7, hy = 1.3, hz = 0.35;
  const GpuArray<double, 3> coarse_dx{{hx, hy, hz}};
  const double hc[3] = {hx, hy, hz};
  const double hf[3] = {hx / 2, hy / 2, hz / 2};

  std::printf("vecpot interpolater (real code, shimmed AMReX)\n\n");

  // ---- CoarseBox ---------------------------------------------------------
  {
    std::printf("[1] CoarseBox\n");
    const Box fine(IntVect(0, 0, 0), IntVect(2 * NC - 1, 2 * NC - 1, 2 * NC - 1));
    const Box c = prolongate_vecpot_3d_rf2.CoarseBox(fine, 2);
    // The stencil reaches base coarse index NC (nodal) plus 2.
    const bool lo_ok = c.smallEnd(0) <= -2;
    const bool hi_ok = c.bigEnd(0) >= NC + 1;
    std::printf("  coarse box for fine [0,%d]^3 : [%d,%d]^3   %s\n", 2 * NC - 1,
                c.smallEnd(0), c.bigEnd(0), (lo_ok && hi_ok) ? "ok" : "FAIL");
    if (!(lo_ok && hi_ok))
      ++failures;
  }

  // ---- allocate --------------------------------------------------------
  const IntVect clo(-CG, -CG, -CG), chi(NC + CG, NC + CG, NC + CG);
  const IntVect flo(-2 * CG, -2 * CG, -2 * CG),
      fhi(2 * (NC + CG), 2 * (NC + CG), 2 * (NC + CG));

  FArrayBox cxb(Box(clo, chi, edge_t(0)), 1), cyb(Box(clo, chi, edge_t(1)), 1),
      czb(Box(clo, chi, edge_t(2)), 1);
  FArrayBox fxb(Box(flo, fhi, edge_t(0)), 1), fyb(Box(flo, fhi, edge_t(1)), 1),
      fzb(Box(flo, fhi, edge_t(2)), 1);

  std::mt19937_64 rng(20260930);
  std::normal_distribution<double> nd(0.0, 1.0);
  for (auto *f : {&cxb, &cyb, &czb})
    for (auto &v : f->data)
      v = nd(rng);

  Array<const FArrayBox *, 3> crse{{&cxb, &cyb, &czb}};
  Array<FArrayBox *, 3> fineA{{&fxb, &fyb, &fzb}};

  // target: fine cells [0, 2*NC-1]^3, expressed cell-centred
  const Box fine_region_cc(IntVect(0, 0, 0),
                           IntVect(2 * NC - 1, 2 * NC - 1, 2 * NC - 1));

  vecpot_masks nomask;
  prolongate_vecpot_3d_rf2.interp_vecpot(crse, 0, fineA, 0, 1, fine_region_cc,
                                         IntVect(2), nomask,
                                         coarse_dx, RunOn::Cpu);

  auto FX = fxb.array(), FY = fyb.array(), FZ = fzb.array();
  auto CX = cxb.const_array(), CY = cyb.const_array(), CZ = czb.const_array();

  auto Bc = [&](int C, int i, int j, int k) {
    if (C == 0)
      return (CY(i, j, k, 0) - CY(i, j, k + 1, 0)) / hc[2] +
             (CZ(i, j + 1, k, 0) - CZ(i, j, k, 0)) / hc[1];
    if (C == 1)
      return (CZ(i, j, k, 0) - CZ(i + 1, j, k, 0)) / hc[0] +
             (CX(i, j, k + 1, 0) - CX(i, j, k, 0)) / hc[2];
    return (CX(i, j, k, 0) - CX(i, j + 1, k, 0)) / hc[1] +
           (CY(i + 1, j, k, 0) - CY(i, j, k, 0)) / hc[0];
  };
  auto Bf = [&](int C, int i, int j, int k) {
    return Bfine(C, FX, FY, FZ, i, j, k, hf);
  };

  const int lo = 2, hi = NC - 3; // coarse cells safely inside

  // ---- restriction compatibility ----------------------------------------
  {
    std::printf("\n[2] identities on the prolongated field\n");
    double e = 0;
    for (int C = 0; C < 3; ++C) {
      const int f1 = (C + 1) % 3, f2 = (C + 2) % 3;
      for (int i = lo; i <= hi; ++i)
        for (int j = lo; j <= hi; ++j)
          for (int k = lo; k <= hi; ++k) {
            int x[3] = {i, j, k};
            double s = 0;
            for (int q1 = 0; q1 < 2; ++q1)
              for (int q2 = 0; q2 < 2; ++q2) {
                int f[3];
                f[C] = 2 * x[C];
                f[f1] = 2 * x[f1] + q1;
                f[f2] = 2 * x[f2] + q2;
                s += Bf(C, f[0], f[1], f[2]);
              }
            e = std::max(e, std::fabs(0.25 * s - Bc(C, i, j, k)));
          }
    }
    check("coarse flux equals the mean of its four fine fluxes", e, 1e-12);

    double ediv = 0;
    for (int i = 2 * lo; i <= 2 * hi; ++i)
      for (int j = 2 * lo; j <= 2 * hi; ++j)
        for (int k = 2 * lo; k <= 2 * hi; ++k)
          ediv = std::max(
              ediv, std::fabs((Bf(0, i + 1, j, k) - Bf(0, i, j, k)) / hf[0] +
                              (Bf(1, i, j + 1, k) - Bf(1, i, j, k)) / hf[1] +
                              (Bf(2, i, j, k + 1) - Bf(2, i, j, k)) / hf[2]));
    check("fine div B", ediv, 1e-10);
  }

  // ---- every targeted fine edge was written -------------------------------
  {
    std::printf("\n[3] coverage of the target region\n");
    // Re-run into a poisoned destination and count untouched points.
    for (auto *f : {&fxb, &fyb, &fzb})
      for (auto &v : f->data)
        v = 1e300;
    prolongate_vecpot_3d_rf2.interp_vecpot(crse, 0, fineA, 0, 1, fine_region_cc,
                                           IntVect(2), nomask,
                                           coarse_dx, RunOn::Cpu);
    long unwritten = 0, total = 0;
    for (int d = 0; d < 3; ++d) {
      auto A = (d == 0) ? fxb.array() : (d == 1) ? fyb.array() : fzb.array();
      int hiv[3] = {2 * NC, 2 * NC, 2 * NC};
      hiv[d] = 2 * NC - 1; // cell-centred in its own direction
      for (int i = 0; i <= hiv[0]; ++i)
        for (int j = 0; j <= hiv[1]; ++j)
          for (int k = 0; k <= hiv[2]; ++k) {
            ++total;
            if (A(i, j, k, 0) == 1e300)
              ++unwritten;
          }
    }
    std::printf("  %-58s %10ld  %s\n", "fine edges in target left unwritten",
                unwritten, unwritten == 0 ? "ok" : "FAIL");
    std::printf("  %-58s %10ld\n", "fine edges in target", total);
    if (unwritten)
      ++failures;
  }

  // ---- pure gauge --------------------------------------------------------
  {
    std::printf("\n[4] pure gauge\n");
    FArrayBox lam(Box(clo, chi, IndexType::TheNodeType()), 1);
    for (auto &v : lam.data)
      v = nd(rng);
    auto L = lam.array();
    for (int i = -CG; i <= NC + CG - 1; ++i)
      for (int j = -CG; j <= NC + CG - 1; ++j)
        for (int k = -CG; k <= NC + CG - 1; ++k) {
          cxb.array()(i, j, k, 0) = (L(i + 1, j, k, 0) - L(i, j, k, 0)) / hc[0];
          cyb.array()(i, j, k, 0) = (L(i, j + 1, k, 0) - L(i, j, k, 0)) / hc[1];
          czb.array()(i, j, k, 0) = (L(i, j, k + 1, 0) - L(i, j, k, 0)) / hc[2];
        }
    prolongate_vecpot_3d_rf2.interp_vecpot(crse, 0, fineA, 0, 1, fine_region_cc,
                                           IntVect(2), nomask,
                                           coarse_dx, RunOn::Cpu);
    double e = 0;
    for (int C = 0; C < 3; ++C)
      for (int i = 2 * lo; i <= 2 * hi; ++i)
        for (int j = 2 * lo; j <= 2 * hi; ++j)
          for (int k = 2 * lo; k <= 2 * hi; ++k)
            e = std::max(e, std::fabs(Bf(C, i, j, k)));
    check("fine B from a pure gauge coarse field", e, 1e-10);
  }

  // ---- masks -------------------------------------------------------------
  {
    std::printf("\n[5] solve masks\n");
    for (auto *f : {&cxb, &cyb, &czb})
      for (auto &v : f->data)
        v = nd(rng);
    // Everything known -> nothing may be written.
    IArrayBox me0(Box(clo, chi, edge_t(0)), 1), me1(Box(clo, chi, edge_t(1)), 1),
        me2(Box(clo, chi, edge_t(2)), 1);
    IArrayBox mf0(Box(clo, chi, IndexType(IntVect::TheDimensionVector(0))), 1),
        mf1(Box(clo, chi, IndexType(IntVect::TheDimensionVector(1))), 1),
        mf2(Box(clo, chi, IndexType(IntVect::TheDimensionVector(2))), 1);
    IArrayBox mcell(Box(clo, chi, IndexType::TheCellType()), 1);
    vecpot_masks allknown;
    allknown.edge[0] = &me0;
    allknown.edge[1] = &me1;
    allknown.edge[2] = &me2;
    allknown.face[0] = &mf0;
    allknown.face[1] = &mf1;
    allknown.face[2] = &mf2;
    allknown.cell = &mcell;
    // masks default to 0 == known
    for (auto *f : {&fxb, &fyb, &fzb})
      for (auto &v : f->data)
        v = -12345.0;
    prolongate_vecpot_3d_rf2.interp_vecpot(crse, 0, fineA, 0, 1, fine_region_cc,
                                           IntVect(2), allknown,
                                           coarse_dx, RunOn::Cpu);
    double e = 0;
    for (auto *f : {&fxb, &fyb, &fzb})
      for (auto &v : f->data)
        e = std::max(e, std::fabs(v - (-12345.0)));
    check("all-known mask: nothing written", e, 0);

    // All unknown -> identical to the no-mask run.
    for (auto *m : {&me0, &me1, &me2, &mf0, &mf1, &mf2, &mcell})
      for (auto &v : m->data)
        v = 1;
    for (auto *f : {&fxb, &fyb, &fzb})
      for (auto &v : f->data)
        v = 0;
    prolongate_vecpot_3d_rf2.interp_vecpot(crse, 0, fineA, 0, 1, fine_region_cc,
                                           IntVect(2), allknown,
                                           coarse_dx, RunOn::Cpu);
    std::vector<double> ref[3] = {fxb.data, fyb.data, fzb.data};
    for (auto *f : {&fxb, &fyb, &fzb})
      for (auto &v : f->data)
        v = 0;
    prolongate_vecpot_3d_rf2.interp_vecpot(crse, 0, fineA, 0, 1, fine_region_cc,
                                           IntVect(2), nomask,
                                           coarse_dx, RunOn::Cpu);
    double e2 = 0;
    const std::vector<double> *got[3] = {&fxb.data, &fyb.data, &fzb.data};
    for (int d = 0; d < 3; ++d)
      for (size_t n = 0; n < ref[d].size(); ++n)
        e2 = std::max(e2, std::fabs(ref[d][n] - (*got[d])[n]));
    check("all-unknown mask == null mask", e2, 0);
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
