#include "prolongate_3d_rf2_vecpot.hxx"

#include <cmath>
#include <random>
#include <vector>

// Self-test for the vector-potential prolongation, run once at startup when
// the operator is selected.  It is a compact version of the standalone tests
// in CarpetX/test/unit/, which cover the same ground in more detail and can be
// run without an Einstein Toolkit build.
//
// The discrete curls used for verification are written out longhand, not taken
// from the kernel header, so that an error in the cyclic direction algebra
// cannot cancel against itself.

namespace CarpetX {
using namespace amrex;

namespace {

constexpr int nc = 6;  // coarse cells in the target region
constexpr int ng = 4;  // coarse ghosts of the scratch allocation

IndexType edge_ixtype_t(const int d) {
  return IndexType(IntVect::TheNodeVector() - IntVect::TheDimensionVector(d));
}

struct curls {
  Array4<const CCTK_REAL> x, y, z;
  CCTK_REAL h[3];
  CCTK_REAL operator()(const int C, const int i, const int j,
                       const int k) const {
    if (C == 0)
      return (y(i, j, k) - y(i, j, k + 1)) / h[2] +
             (z(i, j + 1, k) - z(i, j, k)) / h[1];
    if (C == 1)
      return (z(i, j, k) - z(i + 1, j, k)) / h[0] +
             (x(i, j, k + 1) - x(i, j, k)) / h[2];
    return (x(i, j, k) - x(i, j + 1, k)) / h[1] +
           (y(i + 1, j, k) - y(i, j, k)) / h[0];
  }
};

} // namespace

void test_prolongate_3d_rf2_vecpot() {
  const CCTK_REAL hx = 0.7, hy = 1.3, hz = 0.35; // non-unit aspect ratio
  const GpuArray<CCTK_REAL, dim> coarse_dx{hx, hy, hz};

  const IntVect clo(-ng), chi(nc + ng);
  const IntVect flo(-2 * ng), fhi(2 * (nc + ng));

  std::vector<FArrayBox> cf, ff;
  cf.reserve(dim);
  ff.reserve(dim);
  for (int d = 0; d < dim; ++d) {
    cf.emplace_back(Box(clo, chi, edge_ixtype_t(d)), 1);
    ff.emplace_back(Box(flo, fhi, edge_ixtype_t(d)), 1);
  }

  std::mt19937_64 rng(20260930);
  std::normal_distribution<CCTK_REAL> nd(0.0, 1.0);

  Array<const FArrayBox *, dim> crse{&cf[0], &cf[1], &cf[2]};
  Array<FArrayBox *, dim> fine{&ff[0], &ff[1], &ff[2]};
  const Box region(IntVect(0), IntVect(2 * nc - 1)); // cell-centred

  const CCTK_REAL hc[3] = {hx, hy, hz};
  const CCTK_REAL hf[3] = {hx / 2, hy / 2, hz / 2};
  const int lo = 2, hi = nc - 3;

  CCTK_REAL worst_flux = 0, worst_div = 0, worst_gauge = 0;

  for (int pass = 0; pass < 2; ++pass) {
    if (pass == 0) {
      // random coarse data
      for (int d = 0; d < dim; ++d)
        for (auto *p = cf[d].dataPtr(); p != cf[d].dataPtr() + cf[d].size(); ++p)
          *p = nd(rng);
    } else {
      // pure gauge: A = grad lambda
      FArrayBox lam(Box(clo, chi, IndexType::TheNodeType()), 1);
      for (auto *p = lam.dataPtr(); p != lam.dataPtr() + lam.size(); ++p)
        *p = nd(rng);
      const auto L = lam.const_array();
      const auto AX = cf[0].array(), AY = cf[1].array(), AZ = cf[2].array();
      for (int i = -ng; i <= nc + ng - 1; ++i)
        for (int j = -ng; j <= nc + ng - 1; ++j)
          for (int k = -ng; k <= nc + ng - 1; ++k) {
            AX(i, j, k) = (L(i + 1, j, k) - L(i, j, k)) / hx;
            AY(i, j, k) = (L(i, j + 1, k) - L(i, j, k)) / hy;
            AZ(i, j, k) = (L(i, j, k + 1) - L(i, j, k)) / hz;
          }
    }

    prolongate_vecpot_3d_rf2.interp_vecpot(crse, 0, fine, 0, 1, region,
                                           IntVect(2), vecpot_masks{},
                                           coarse_dx, RunOn::Cpu);

    const curls Bc{cf[0].const_array(), cf[1].const_array(),
                   cf[2].const_array(), {hc[0], hc[1], hc[2]}};
    const curls Bf{ff[0].const_array(), ff[1].const_array(),
                   ff[2].const_array(), {hf[0], hf[1], hf[2]}};

    // Flux conservation: each coarse face equals the mean of its four fine
    // faces (draft Eq. 7), i.e. the prolongation is compatible with the
    // restriction average_down_edges performs.
    for (int C = 0; C < dim; ++C) {
      const int f1 = (C + 1) % 3, f2 = (C + 2) % 3;
      for (int i = lo; i <= hi; ++i)
        for (int j = lo; j <= hi; ++j)
          for (int k = lo; k <= hi; ++k) {
            const int x[3] = {i, j, k};
            CCTK_REAL s = 0;
            for (int q1 = 0; q1 < 2; ++q1)
              for (int q2 = 0; q2 < 2; ++q2) {
                int f[3];
                f[C] = 2 * x[C];
                f[f1] = 2 * x[f1] + q1;
                f[f2] = 2 * x[f2] + q2;
                s += Bf(C, f[0], f[1], f[2]);
              }
            using std::fabs;
            worst_flux = std::max(worst_flux,
                                  fabs(0.25 * s - Bc(C, i, j, k)));
          }
    }

    // The fine field is a discrete curl, so its divergence vanishes
    // identically; this checks the bookkeeping, not the operator.
    for (int i = 2 * lo; i <= 2 * hi; ++i)
      for (int j = 2 * lo; j <= 2 * hi; ++j)
        for (int k = 2 * lo; k <= 2 * hi; ++k) {
          using std::fabs;
          worst_div = std::max(
              worst_div,
              fabs((Bf(0, i + 1, j, k) - Bf(0, i, j, k)) / hf[0] +
                   (Bf(1, i, j + 1, k) - Bf(1, i, j, k)) / hf[1] +
                   (Bf(2, i, j, k + 1) - Bf(2, i, j, k)) / hf[2]));
        }

    // A pure gauge coarse field must produce no fine magnetic field at all.
    if (pass == 1)
      for (int C = 0; C < dim; ++C)
        for (int i = 2 * lo; i <= 2 * hi; ++i)
          for (int j = 2 * lo; j <= 2 * hi; ++j)
            for (int k = 2 * lo; k <= 2 * hi; ++k) {
              using std::fabs;
              worst_gauge = std::max(worst_gauge, fabs(Bf(C, i, j, k)));
            }
  }

  const CCTK_REAL tol = 1.0e-10;
  if (!(worst_flux <= tol) || !(worst_div <= tol) || !(worst_gauge <= tol))
    CCTK_VERROR("prolongate_3d_rf2_vecpot self-test failed: flux conservation "
                "%.3g, div B %.3g, pure gauge %.3g (tolerance %.3g)",
                double(worst_flux), double(worst_div), double(worst_gauge),
                double(tol));

  CCTK_VINFO("prolongate_3d_rf2_vecpot self-test passed "
             "(flux %.3g, div B %.3g, gauge %.3g)",
             double(worst_flux), double(worst_div), double(worst_gauge));
}

} // namespace CarpetX
