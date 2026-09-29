// Standalone unit test for the vector-potential prolongation arithmetic in
// CarpetX/src/prolongate_3d_rf2_vecpot_kernels.hxx.
//
// This test has no AMReX and no Cactus dependency; it can be built and run
// with a bare C++17 compiler.  See CarpetX/test/unit/README.md.
//
// It is the C++ port of the numerical verification recorded in PROLONG_A.md
// section 2.7, and it checks, at a deliberately non-unit cell aspect ratio:
//
//   1  the four owning-type boxes partition the fine edges exactly
//   2  restriction compatibility, draft Eq. (7)
//   3  the Toth surface fluxes, draft Eq. (11)
//   4  the Toth interior fluxes, draft Eq. (15)
//   5  a pure gauge field prolongs to zero magnetic field
//   6  gauge covariance: prolonged A is the gradient of trilinear lambda
//   7  the fine field is exactly divergence free
//   8  masked (already refined) edges are never overwritten
//   9  a coarse face with mixed provenance still conserves flux
//  10  second-order accuracy of B on a smooth field
//
// The discrete curls used for *verification* are written out longhand below,
// independently of the cyclic templates in the header, so that a mistake in
// the cyclic direction algebra cannot cancel against itself.

#include "../../src/prolongate_3d_rf2_vecpot_kernels.hxx"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace CarpetX::vecpot;

////////////////////////////////////////////////////////////////////////////////
// A minimal Array4-alike

template <typename T> struct Field {
  std::vector<T> data;
  int lo[3]{}, hi[3]{};
  int ncomp = 1;

  Field() = default;
  Field(int l, int h, int nc = 1) : ncomp(nc) {
    for (int d = 0; d < 3; ++d) {
      lo[d] = l;
      hi[d] = h;
    }
    data.assign(size_t(n(0)) * n(1) * n(2) * ncomp, T(0));
  }
  int n(int d) const { return hi[d] - lo[d] + 1; }
  size_t off(int i, int j, int k, int c) const {
    const size_t a = size_t(i - lo[0]);
    const size_t b = size_t(j - lo[1]);
    const size_t e = size_t(k - lo[2]);
    if (i < lo[0] || i > hi[0] || j < lo[1] || j > hi[1] || k < lo[2] ||
        k > hi[2]) {
      std::fprintf(stderr, "index (%d,%d,%d) out of range [%d,%d]\n", i, j, k,
                   lo[0], hi[0]);
      std::abort();
    }
    return ((size_t(c) * n(2) + e) * n(1) + b) * n(0) + a;
  }
  T &operator()(int i, int j, int k, int c = 0) { return data[off(i, j, k, c)]; }
  const T &operator()(int i, int j, int k, int c = 0) const {
    return data[off(i, j, k, c)];
  }
};

// A read-only view with the accessor signature the kernels expect, plus the
// null-testability that the mask guard relies on.
template <typename T> struct View {
  Field<T> *f = nullptr;
  T &operator()(int i, int j, int k, int c) const { return (*f)(i, j, k, c); }
  explicit operator bool() const { return f != nullptr; }
};

using DView = View<double>;
using IView = View<int>;

struct Arrays {
  DView a[3];
  const DView &operator[](int d) const { return a[d]; }
};

////////////////////////////////////////////////////////////////////////////////
// Geometry

static constexpr int NC = 8;  // coarse cells per direction in the target box
static constexpr int CG = 4;  // coarse ghost width of the allocation
static double hc[3] = {0.7, 1.3, 0.35}; // deliberately non-unit aspect ratio

static double hf(int d) { return hc[d] / 2; }

struct Grids {
  Field<double> cx, cy, cz; // coarse A_x, A_y, A_z
  Field<double> fx, fy, fz; // fine   A_x, A_y, A_z
  Grids()
      : cx(-CG, NC + CG), cy(-CG, NC + CG), cz(-CG, NC + CG),
        fx(-2 * CG, 2 * (NC + CG)), fy(-2 * CG, 2 * (NC + CG)),
        fz(-2 * CG, 2 * (NC + CG)) {}
  Arrays coarse() { return Arrays{{DView{&cx}, DView{&cy}, DView{&cz}}}; }
  Arrays fine() { return Arrays{{DView{&fx}, DView{&fy}, DView{&fz}}}; }
};

// Verification curls, written out longhand (see file header).
static double Bx_of(const Field<double> &ay, const Field<double> &az, int i,
                    int j, int k, const double h[3]) {
  return (ay(i, j, k) - ay(i, j, k + 1)) / h[2] +
         (az(i, j + 1, k) - az(i, j, k)) / h[1];
}
static double By_of(const Field<double> &az, const Field<double> &ax, int i,
                    int j, int k, const double h[3]) {
  return (az(i, j, k) - az(i + 1, j, k)) / h[0] +
         (ax(i, j, k + 1) - ax(i, j, k)) / h[2];
}
static double Bz_of(const Field<double> &ax, const Field<double> &ay, int i,
                    int j, int k, const double h[3]) {
  return (ax(i, j, k) - ax(i, j + 1, k)) / h[1] +
         (ay(i + 1, j, k) - ay(i, j, k)) / h[0];
}
static double Bc(const Grids &g, int C, int i, int j, int k) {
  const double h[3] = {hc[0], hc[1], hc[2]};
  return C == 0 ? Bx_of(g.cy, g.cz, i, j, k, h)
                : C == 1 ? By_of(g.cz, g.cx, i, j, k, h)
                         : Bz_of(g.cx, g.cy, i, j, k, h);
}
static double Bfn(const Grids &g, int C, int i, int j, int k) {
  const double h[3] = {hf(0), hf(1), hf(2)};
  return C == 0 ? Bx_of(g.fy, g.fz, i, j, k, h)
                : C == 1 ? By_of(g.fz, g.fx, i, j, k, h)
                         : Bz_of(g.fx, g.fy, i, j, k, h);
}

////////////////////////////////////////////////////////////////////////////////
// The three launches, over the owning-type boxes (PROLONG_A.md section 5.2)

struct Box {
  int lo[3], hi[3];
};

// Visit every (component D, parities P1 P2, owning index x) the prolongation of
// cell box `cb` is responsible for.  Both the prolongation driver and the
// partition check go through this, so they cannot disagree.
template <typename F> static void foreach_owned(const Box &cb, F &&f) {
  for (int D = 0; D < 3; ++D) {
    const int b1 = (D + 1) % 3, b2 = (D + 2) % 3;
    for (int P1 = 0; P1 < 2; ++P1)
      for (int P2 = 0; P2 < 2; ++P2) {
        // Owning type: CELL in D; CELL in b1 iff P1; CELL in b2 iff P2.
        // convert(cell box, that type) extends by one wherever it is NODE.
        int lo[3], hi[3];
        for (int d = 0; d < 3; ++d) {
          lo[d] = cb.lo[d];
          hi[d] = cb.hi[d];
        }
        if (!P1)
          hi[b1] += 1;
        if (!P2)
          hi[b2] += 1;
        for (int i = lo[0]; i <= hi[0]; ++i)
          for (int j = lo[1]; j <= hi[1]; ++j)
            for (int k = lo[2]; k <= hi[2]; ++k)
              f(D, P1, P2, idx3{{i, j, k}});
      }
  }
}

// Masks, one per owning object type, indexed exactly like its launch.
struct Masks {
  Field<int> *m[3][2][2] = {};
  const IView get(int D, int P1, int P2) const {
    Field<int> *p = m[D][P1][P2];
    return p ? IView{p} : IView{};
  }
};

static void prolong(const Box &cb, Grids &g, const Masks &mk) {
  Arrays ac = g.coarse(), af = g.fine();
  // Launches 1 and 2: everything that reads coarse data.
  foreach_owned(cb, [&](int D, int P1, int P2, const idx3 &x) {
    if (P1 && P2)
      return; // launch 3
    const IView mask = mk.get(D, P1, P2);
    switch (D * 4 + P1 * 2 + P2) {
    case 0 * 4 + 0 * 2 + 0: surface_interp<0,0,0>(ac, af, mask, x, 0, hc); break;
    case 0 * 4 + 1 * 2 + 0: surface_interp<0,1,0>(ac, af, mask, x, 0, hc); break;
    case 0 * 4 + 0 * 2 + 1: surface_interp<0,0,1>(ac, af, mask, x, 0, hc); break;
    case 1 * 4 + 0 * 2 + 0: surface_interp<1,0,0>(ac, af, mask, x, 0, hc); break;
    case 1 * 4 + 1 * 2 + 0: surface_interp<1,1,0>(ac, af, mask, x, 0, hc); break;
    case 1 * 4 + 0 * 2 + 1: surface_interp<1,0,1>(ac, af, mask, x, 0, hc); break;
    case 2 * 4 + 0 * 2 + 0: surface_interp<2,0,0>(ac, af, mask, x, 0, hc); break;
    case 2 * 4 + 1 * 2 + 0: surface_interp<2,1,0>(ac, af, mask, x, 0, hc); break;
    case 2 * 4 + 0 * 2 + 1: surface_interp<2,0,1>(ac, af, mask, x, 0, hc); break;
    default: std::abort();
    }
  });
  // Launch 3: the cell interiors, reading only the fine array.
  foreach_owned(cb, [&](int D, int P1, int P2, const idx3 &x) {
    if (!(P1 && P2))
      return;
    const IView mask = mk.get(D, P1, P2);
    switch (D) {
    case 0: interior_interp<0>(af, mask, x, 0, hc); break;
    case 1: interior_interp<1>(af, mask, x, 0, hc); break;
    case 2: interior_interp<2>(af, mask, x, 0, hc); break;
    default: std::abort();
    }
  });
}

////////////////////////////////////////////////////////////////////////////////
// Harness

static int failures = 0;
static void check(const char *name, double err, double tol) {
  const bool ok = std::isfinite(err) && err <= tol;
  std::printf("  %-58s %10.3e  %s\n", name, err, ok ? "ok" : "FAIL");
  if (!ok)
    ++failures;
}

static std::mt19937_64 rng(20260929);
static double rnd() {
  std::normal_distribution<double> d(0.0, 1.0);
  return d(rng);
}

static void fill_random(Grids &g) {
  for (auto *f : {&g.cx, &g.cy, &g.cz})
    for (auto &v : f->data)
      v = rnd();
}

static Box target() { return Box{{0, 0, 0}, {NC - 1, NC - 1, NC - 1}}; }

// The sub-box of coarse cells whose full +-2 stencil and whose fine surface are
// safely inside the prolongated region.
static Box inner() { return Box{{2, 2, 2}, {NC - 3, NC - 3, NC - 3}}; }

////////////////////////////////////////////////////////////////////////////////

int main() {
  std::printf("vecpot prolongation kernels: dX,dY,dZ = %g, %g, %g\n\n", hc[0],
              hc[1], hc[2]);

  const Box cb = target(), ib = inner();

  // ---- 1  exact partition of the fine edges -------------------------------
  {
    std::printf("[1] owning-type boxes partition the fine edges\n");
    std::map<std::string, int> count;
    long total = 0;
    foreach_owned(cb, [&](int D, int P1, int P2, const idx3 &x) {
      const int b1 = (D + 1) % 3, b2 = (D + 2) % 3;
      idx3 f{{0, 0, 0}};
      f[D] = 2 * x[D];
      f[b1] = 2 * x[b1] + P1;
      f[b2] = 2 * x[b2] + P2;
      for (int h = 0; h < 2; ++h) {
        char buf[64];
        std::snprintf(buf, sizeof buf, "%d:%d,%d,%d", D, f[0] + (D == 0 ? h : 0),
                      f[1] + (D == 1 ? h : 0), f[2] + (D == 2 ? h : 0));
        ++count[buf];
        ++total;
      }
    });
    int multi = 0;
    for (const auto &kv : count)
      if (kv.second != 1)
        ++multi;
    // Expected: per component, a snapped patch of NC^3 coarse cells holds
    // (2*NC) * (2*NC+1) * (2*NC+1) fine edges.
    const long expect = 3L * (2 * NC) * (2 * NC + 1) * (2 * NC + 1);
    std::printf("  %-58s %10ld  %s\n", "distinct fine edges written",
                long(count.size()), count.size() == size_t(expect) ? "ok" : "FAIL");
    if (count.size() != size_t(expect))
      ++failures;
    check("edges written more than once", multi, 0);
    check("writes minus distinct edges", double(total - long(count.size())), 0);
  }

  // ---- random coarse data, plain prolongation -----------------------------
  Grids g;
  fill_random(g);
  Masks nomask;
  prolong(cb, g, nomask);

  // ---- 2  restriction compatibility, Eq. (7) ------------------------------
  {
    std::printf("\n[2] restriction compatibility (Eq. 7) and (Eq. 11), (Eq. 15)\n");
    double e7 = 0;
    for (int C = 0; C < 3; ++C) {
      const int f1 = (C + 1) % 3, f2 = (C + 2) % 3;
      for (int i = ib.lo[0]; i <= ib.hi[0]; ++i)
        for (int j = ib.lo[1]; j <= ib.hi[1]; ++j)
          for (int k = ib.lo[2]; k <= ib.hi[2]; ++k) {
            const idx3 x{{i, j, k}};
            for (int side = 0; side < 2; ++side) {
              idx3 xs = shifted(x, C, side);
              double s = 0;
              for (int q1 = 0; q1 < 2; ++q1)
                for (int q2 = 0; q2 < 2; ++q2) {
                  idx3 f{{0, 0, 0}};
                  f[C] = 2 * xs[C];
                  f[f1] = 2 * xs[f1] + q1;
                  f[f2] = 2 * xs[f2] + q2;
                  s += Bfn(g, C, f[0], f[1], f[2]);
                }
              e7 = std::max(e7, std::fabs(0.25 * s - Bc(g, C, xs[0], xs[1], xs[2])));
            }
          }
    }
    check("coarse flux equals the mean of its four fine fluxes", e7, 1e-12);

    // ---- 3  Toth surface fluxes, Eq. (11) ---------------------------------
    double e11 = 0;
    for (int C = 0; C < 3; ++C) {
      const int f1 = (C + 1) % 3, f2 = (C + 2) % 3;
      for (int i = ib.lo[0]; i <= ib.hi[0]; ++i)
        for (int j = ib.lo[1]; j <= ib.hi[1]; ++j)
          for (int k = ib.lo[2]; k <= ib.hi[2]; ++k) {
            const idx3 x{{i, j, k}};
            for (int side = 0; side < 2; ++side) {
              const idx3 xs = shifted(x, C, side);
              const idx3 p1 = shifted(xs, f1, +1), m1 = shifted(xs, f1, -1);
              const idx3 p2 = shifted(xs, f2, +1), m2 = shifted(xs, f2, -1);
              const double b0 = Bc(g, C, xs[0], xs[1], xs[2]);
              const double s1 = Bc(g, C, p1[0], p1[1], p1[2]) -
                                Bc(g, C, m1[0], m1[1], m1[2]);
              const double s2 = Bc(g, C, p2[0], p2[1], p2[2]) -
                                Bc(g, C, m2[0], m2[1], m2[2]);
              for (int q1 = 0; q1 < 2; ++q1)
                for (int q2 = 0; q2 < 2; ++q2) {
                  idx3 f{{0, 0, 0}};
                  f[C] = 2 * xs[C];
                  f[f1] = 2 * xs[f1] + q1;
                  f[f2] = 2 * xs[f2] + q2;
                  const double ref =
                      b0 + (2 * q1 - 1) / 8.0 * s1 + (2 * q2 - 1) / 8.0 * s2;
                  e11 = std::max(e11,
                                 std::fabs(Bfn(g, C, f[0], f[1], f[2]) - ref));
                }
            }
          }
    }
    check("fine surface fluxes match Toth Eq. (11)", e11, 1e-12);

    // ---- 4  Toth interior fluxes, Eq. (15) --------------------------------
    double e15 = 0;
    for (int C = 0; C < 3; ++C) {
      const int c1 = (C + 1) % 3, c2 = (C + 2) % 3;
      for (int i = ib.lo[0]; i <= ib.hi[0]; ++i)
        for (int j = ib.lo[1]; j <= ib.hi[1]; ++j)
          for (int k = ib.lo[2]; k <= ib.hi[2]; ++k) {
            const idx3 x{{i, j, k}};
            // moments of the cell's own surface fluxes
            double S1 = 0, S2 = 0;
            for (int eC = 0; eC < 2; ++eC)
              for (int s = 0; s < 2; ++s)
                for (int e = 0; e < 2; ++e) {
                  idx3 fa{{0, 0, 0}};
                  fa[C] = 2 * x[C] + eC;
                  fa[c1] = 2 * (x[c1] + s);
                  fa[c2] = 2 * x[c2] + e;
                  S1 += (2 * eC - 1) * (2 * s - 1) *
                        Bfn(g, c1, fa[0], fa[1], fa[2]);
                  idx3 fb{{0, 0, 0}};
                  fb[C] = 2 * x[C] + eC;
                  fb[c1] = 2 * x[c1] + e;
                  fb[c2] = 2 * (x[c2] + s);
                  S2 += (2 * eC - 1) * (2 * s - 1) *
                        Bfn(g, c2, fb[0], fb[1], fb[2]);
                }
            for (int q1 = 0; q1 < 2; ++q1)
              for (int q2 = 0; q2 < 2; ++q2) {
                idx3 fl{{0, 0, 0}}, fr{{0, 0, 0}}, fm{{0, 0, 0}};
                for (int d = 0; d < 3; ++d)
                  fl[d] = fr[d] = fm[d] = 0;
                fl[C] = 2 * x[C];
                fr[C] = 2 * x[C] + 2;
                fm[C] = 2 * x[C] + 1;
                fl[c1] = fr[c1] = fm[c1] = 2 * x[c1] + q1;
                fl[c2] = fr[c2] = fm[c2] = 2 * x[c2] + q2;
                const double ref =
                    0.5 * (Bfn(g, C, fl[0], fl[1], fl[2]) +
                           Bfn(g, C, fr[0], fr[1], fr[2])) +
                    hc[C] / (8 * hc[c1]) * S1 + hc[C] / (8 * hc[c2]) * S2;
                e15 =
                    std::max(e15, std::fabs(Bfn(g, C, fm[0], fm[1], fm[2]) - ref));
              }
          }
    }
    check("fine interior fluxes match Toth Eq. (15)", e15, 1e-12);
  }

  // ---- 7  divergence ------------------------------------------------------
  {
    std::printf("\n[3] divergence and gauge\n");
    double ediv = 0;
    for (int i = 2 * ib.lo[0]; i <= 2 * ib.hi[0] + 1; ++i)
      for (int j = 2 * ib.lo[1]; j <= 2 * ib.hi[1] + 1; ++j)
        for (int k = 2 * ib.lo[2]; k <= 2 * ib.hi[2] + 1; ++k)
          ediv = std::max(ediv, std::fabs((Bfn(g, 0, i + 1, j, k) -
                                           Bfn(g, 0, i, j, k)) / hf(0) +
                                          (Bfn(g, 1, i, j + 1, k) -
                                           Bfn(g, 1, i, j, k)) / hf(1) +
                                          (Bfn(g, 2, i, j, k + 1) -
                                           Bfn(g, 2, i, j, k)) / hf(2)));
    check("fine div B (random coarse data)", ediv, 1e-10);
  }

  // ---- 5, 6  pure gauge ---------------------------------------------------
  {
    Grids gg;
    Field<double> lam(-CG, NC + CG);
    for (auto &v : lam.data)
      v = rnd();
    for (int i = -CG; i <= NC + CG - 1; ++i)
      for (int j = -CG; j <= NC + CG - 1; ++j)
        for (int k = -CG; k <= NC + CG - 1; ++k) {
          gg.cx(i, j, k) = (lam(i + 1, j, k) - lam(i, j, k)) / hc[0];
          gg.cy(i, j, k) = (lam(i, j + 1, k) - lam(i, j, k)) / hc[1];
          gg.cz(i, j, k) = (lam(i, j, k + 1) - lam(i, j, k)) / hc[2];
        }
    prolong(cb, gg, nomask);

    double eB = 0;
    for (int C = 0; C < 3; ++C)
      for (int i = 2 * ib.lo[0]; i <= 2 * ib.hi[0]; ++i)
        for (int j = 2 * ib.lo[1]; j <= 2 * ib.hi[1]; ++j)
          for (int k = 2 * ib.lo[2]; k <= 2 * ib.hi[2]; ++k)
            eB = std::max(eB, std::fabs(Bfn(gg, C, i, j, k)));
    check("fine B from a pure gauge coarse field", eB, 1e-10);

    // trilinear interpolant of lambda on the fine vertices
    Field<double> lf(-2 * CG, 2 * (NC + CG - 1));
    for (int I = lf.lo[0]; I <= lf.hi[0]; ++I)
      for (int J = lf.lo[1]; J <= lf.hi[1]; ++J)
        for (int K = lf.lo[2]; K <= lf.hi[2]; ++K) {
          const int i = I >> 1, oi = I & 1, j = J >> 1, oj = J & 1, k = K >> 1,
                    ok = K & 1;
          if (i + oi > NC + CG || j + oj > NC + CG || k + ok > NC + CG)
            continue;
          double s = 0;
          for (int a = 0; a <= oi; ++a)
            for (int b = 0; b <= oj; ++b)
              for (int c = 0; c <= ok; ++c)
                s += lam(i + a, j + b, k + c);
          lf(I, J, K) = s / ((oi + 1) * (oj + 1) * (ok + 1));
        }
    double eg = 0;
    for (int i = 2 * ib.lo[0]; i <= 2 * ib.hi[0]; ++i)
      for (int j = 2 * ib.lo[1]; j <= 2 * ib.hi[1]; ++j)
        for (int k = 2 * ib.lo[2]; k <= 2 * ib.hi[2]; ++k) {
          eg = std::max(eg, std::fabs(gg.fx(i, j, k) -
                                      (lf(i + 1, j, k) - lf(i, j, k)) / hf(0)));
          eg = std::max(eg, std::fabs(gg.fy(i, j, k) -
                                      (lf(i, j + 1, k) - lf(i, j, k)) / hf(1)));
          eg = std::max(eg, std::fabs(gg.fz(i, j, k) -
                                      (lf(i, j, k + 1) - lf(i, j, k)) / hf(2)));
        }
    check("prolonged A equals grad of trilinear lambda (Eq. 31)", eg, 1e-10);
  }

  // ---- 8, 9  masked edges and mixed faces ---------------------------------
  {
    std::printf("\n[4] copy-from-fine rule\n");
    Grids gm;
    // A "fine truth" whose restriction is the coarse data, so that the copied
    // halves satisfy (C2) exactly.
    Grids ft;
    for (auto *f : {&ft.fx, &ft.fy, &ft.fz})
      for (auto &v : f->data)
        v = rnd();
    for (int i = -CG; i <= NC + CG - 1; ++i)
      for (int j = -CG; j <= NC + CG - 1; ++j)
        for (int k = -CG; k <= NC + CG - 1; ++k) {
          gm.cx(i, j, k) = 0.5 * (ft.fx(2 * i, 2 * j, 2 * k) +
                                  ft.fx(2 * i + 1, 2 * j, 2 * k));
          gm.cy(i, j, k) = 0.5 * (ft.fy(2 * i, 2 * j, 2 * k) +
                                  ft.fy(2 * i, 2 * j + 1, 2 * k));
          gm.cz(i, j, k) = 0.5 * (ft.fz(2 * i, 2 * j, 2 * k) +
                                  ft.fz(2 * i, 2 * j, 2 * k + 1));
        }

    // Mark one coarse A_y edge as already covered by fine data, pre-fill the
    // destination there, and check the prolongation leaves it alone.
    Field<int> m_y00(-CG, NC + CG);
    for (auto &v : m_y00.data)
      v = 1;
    const int i0 = 4, j0 = 4, k0 = 4;
    m_y00(i0, j0, k0) = 0; // known
    Masks mk;
    mk.m[1][0][0] = &m_y00;
    for (int h = 0; h < 2; ++h)
      gm.fy(2 * i0, 2 * j0 + h, 2 * k0) = ft.fy(2 * i0, 2 * j0 + h, 2 * k0);

    prolong(cb, gm, mk);

    double ekeep = 0;
    for (int h = 0; h < 2; ++h)
      ekeep = std::max(ekeep, std::fabs(gm.fy(2 * i0, 2 * j0 + h, 2 * k0) -
                                        ft.fy(2 * i0, 2 * j0 + h, 2 * k0)));
    check("masked fine edges are not overwritten", ekeep, 0);

    double ediff = 0;
    for (int h = 0; h < 2; ++h)
      ediff = std::max(ediff, std::fabs(ft.fy(2 * i0, 2 * j0 + h, 2 * k0) -
                                        gm.cy(i0, j0, k0)));
    std::printf("  %-58s %10.3e  (must be nonzero)\n",
                "copied halves differ from the equal split by", ediff);
    if (!(ediff > 1e-3))
      ++failures;

    // The four coarse faces bounded by that edge still conserve flux.
    double eseam = 0;
    const int faces[4][4] = {{2, i0 - 1, j0, k0}, {2, i0, j0, k0},
                             {0, i0, j0, k0 - 1}, {0, i0, j0, k0}};
    for (const auto &fc : faces) {
      const int C = fc[0];
      const int f1 = (C + 1) % 3, f2 = (C + 2) % 3;
      const idx3 x{{fc[1], fc[2], fc[3]}};
      double s = 0;
      for (int q1 = 0; q1 < 2; ++q1)
        for (int q2 = 0; q2 < 2; ++q2) {
          idx3 f{{0, 0, 0}};
          f[C] = 2 * x[C];
          f[f1] = 2 * x[f1] + q1;
          f[f2] = 2 * x[f2] + q2;
          s += Bfn(gm, C, f[0], f[1], f[2]);
        }
      eseam = std::max(eseam, std::fabs(0.25 * s - Bc(gm, C, x[0], x[1], x[2])));
    }
    check("mixed-provenance coarse faces still conserve flux", eseam, 1e-12);

    double ediv = 0;
    for (int i = 2 * ib.lo[0]; i <= 2 * ib.hi[0]; ++i)
      for (int j = 2 * ib.lo[1]; j <= 2 * ib.hi[1]; ++j)
        for (int k = 2 * ib.lo[2]; k <= 2 * ib.hi[2]; ++k)
          ediv = std::max(ediv, std::fabs((Bfn(gm, 0, i + 1, j, k) -
                                           Bfn(gm, 0, i, j, k)) / hf(0) +
                                          (Bfn(gm, 1, i, j + 1, k) -
                                           Bfn(gm, 1, i, j, k)) / hf(1) +
                                          (Bfn(gm, 2, i, j, k + 1) -
                                           Bfn(gm, 2, i, j, k)) / hf(2)));
    check("fine div B with mixed provenance", ediv, 1e-10);
  }

  // ---- 10  second-order accuracy -----------------------------------------
  {
    std::printf("\n[5] accuracy\n");
    // A_i = smooth analytic field; measure |B_fine - B_exact| at two
    // resolutions and check the ratio.
    auto Ax = [](double x, double y, double z) {
      return std::sin(1.1 * y) * std::cos(0.7 * z) + 0.3 * x * x;
    };
    auto Ay = [](double x, double y, double z) {
      return std::cos(0.9 * z) * std::sin(1.3 * x) + 0.2 * y * y;
    };
    auto Az = [](double x, double y, double z) {
      return std::sin(0.8 * x) * std::cos(1.2 * y) + 0.1 * z * z;
    };
    // Probe Bz = dAy/dx - dAx/dy against its analytic value.  One component
    // is enough here; the identities above already cover all three.
    double err[2];
    for (int lvl = 0; lvl < 2; ++lvl) {
      const double sc = lvl == 0 ? 1.0 : 0.5;
      double hsave[3];
      for (int d = 0; d < 3; ++d) {
        hsave[d] = hc[d];
        hc[d] *= sc;
      }
      Grids gs;
      for (int i = -CG; i <= NC + CG; ++i)
        for (int j = -CG; j <= NC + CG; ++j)
          for (int k = -CG; k <= NC + CG; ++k) {
            if (i <= NC + CG - 1)
              gs.cx(i, j, k) = Ax((i + 0.5) * hc[0], j * hc[1], k * hc[2]);
            if (j <= NC + CG - 1)
              gs.cy(i, j, k) = Ay(i * hc[0], (j + 0.5) * hc[1], k * hc[2]);
            if (k <= NC + CG - 1)
              gs.cz(i, j, k) = Az(i * hc[0], j * hc[1], (k + 0.5) * hc[2]);
          }
      prolong(cb, gs, nomask);
      // exact Bz on a fine z-face centre: dAy/dx - dAx/dy
      double e = 0;
      for (int i = 2 * ib.lo[0]; i <= 2 * ib.hi[0]; ++i)
        for (int j = 2 * ib.lo[1]; j <= 2 * ib.hi[1]; ++j)
          for (int k = 2 * ib.lo[2]; k <= 2 * ib.hi[2]; ++k) {
            const double X = (i + 0.5) * hf(0), Y = (j + 0.5) * hf(1),
                         Z = k * hf(2);
            const double dAy_dx = std::cos(0.9 * Z) * 1.3 * std::cos(1.3 * X);
            const double dAx_dy = 1.1 * std::cos(1.1 * Y) * std::cos(0.7 * Z);
            e = std::max(e, std::fabs(Bfn(gs, 2, i, j, k) - (dAy_dx - dAx_dy)));
          }
      err[lvl] = e;
      for (int d = 0; d < 3; ++d)
        hc[d] = hsave[d];
    }
    const double order = std::log2(err[0] / err[1]);
    std::printf("  %-58s %10.3e\n", "|Bz - exact| at h", err[0]);
    std::printf("  %-58s %10.3e\n", "|Bz - exact| at h/2", err[1]);
    std::printf("  %-58s %10.3f  %s\n", "observed convergence order", order,
                order > 1.8 ? "ok" : "FAIL");
    if (!(order > 1.8))
      ++failures;
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
              failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
