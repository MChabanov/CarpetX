#ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_VECPOT_KERNELS_HXX
#define CARPETX_CARPETX_PROLONGATE_3D_RF2_VECPOT_KERNELS_HXX

// Arithmetic for the flux-preserving prolongation of an edge-staggered vector
// potential.  See PROLONG_A.md for the derivation and for the verification
// status of every formula below.
//
// This header is deliberately free of AMReX and Cactus dependencies.  It is
// templated on the array accessor so that exactly the same code is exercised
// by the standalone unit test in CarpetX/test/unit/ (which can be built with a
// bare C++ compiler) and by the interpolater in prolongate_3d_rf2_vecpot.cxx.
//
// Conventions
// -----------
// * Components are stored, not circulations.  A_d is the d-component of the
//   vector potential and lives on edges parallel to d.  Index types, in the
//   CarpetX convention (0 = vertex/node, 1 = cell):
//
//       A_x {cvv}   A_y {vcv}   A_z {vvc}
//       B^x {vcc}   B^y {cvc}   B^z {ccv}
//
// * Cyclic helpers: for edge direction D the two transverse directions are
//   B1 = (D+1)%3 and B2 = (D+2)%3.  The "+" correction is along B1 and uses
//   B^{B2}; the "-" correction is along B2 and uses B^{B1}.
//
// * Every coarse lookup uses the index of the *coarse object that owns* the
//   fine edges being written, unshifted.  That is what makes the index types
//   line up automatically (PROLONG_A.md section 5.2): when the parity along a
//   transverse direction is 1 the owning type is CELL there, which is exactly
//   what the B component needs.
//
// * There are no cell-aspect-ratio weights anywhere in this file, and there
//   must never be any.  Their absence is the whole point of the edge
//   formulation; compare amrex::facediv_int, which is the same operator for
//   face-staggered B and carries weights like dx*(2*dz*dz+dx*dx)/(8*dy*...).

namespace CarpetX {
namespace vecpot {

#if defined(CCTK_DEVICE) && defined(CCTK_HOST)
#define CARPETX_VECPOT_HD CCTK_DEVICE CCTK_HOST
#else
#define CARPETX_VECPOT_HD
#endif

////////////////////////////////////////////////////////////////////////////////
// Index helpers

/// A 3-index addressable by a (possibly compile-time) direction.
struct idx3 {
  int v[3];
  CARPETX_VECPOT_HD constexpr int operator[](const int d) const { return v[d]; }
  CARPETX_VECPOT_HD constexpr int &operator[](const int d) { return v[d]; }
};

/// A copy of `x` shifted by `s` in direction `d`.
CARPETX_VECPOT_HD constexpr idx3 shifted(idx3 x, const int d, const int s) {
  x[d] += s;
  return x;
}

/// Read accessor `a` at index `x`, component `n`.
template <typename A>
CARPETX_VECPOT_HD inline auto at(const A &a, const idx3 &x, const int n) {
  return a(x[0], x[1], x[2], n);
}

////////////////////////////////////////////////////////////////////////////////
// Discrete curl (discrete Stokes)

/// B^C on a face of type C, from the three edge arrays `a` with spacing `h`.
///
/// Cyclically, with c1 = (C+1)%3 and c2 = (C+2)%3,
///
///     B^C(x) = [ A_c1(x) - A_c1(x + e_c2) ] / h[c2]
///            + [ A_c2(x + e_c1) - A_c2(x) ] / h[c1]
///
/// which for C = z reads Bz = (Ax(i,j,k) - Ax(i,j+1,k))/hY
///                          + (Ay(i+1,j,k) - Ay(i,j,k))/hX.
/// These are Eqs. (5a)-(5c) of the draft divided by the face areas.  Used on
/// the coarse grid with the coarse spacing, and on the fine grid with h/2.
template <int C, typename AA, typename HH>
CARPETX_VECPOT_HD inline auto curl(const AA &a, const idx3 &x, const int n,
                                   const HH &h) {
  constexpr int c1 = (C + 1) % 3;
  constexpr int c2 = (C + 2) % 3;
  return (at(a[c1], x, n) - at(a[c1], shifted(x, c2, +1), n)) / h[c2] +
         (at(a[c2], shifted(x, c1, +1), n) - at(a[c2], x, n)) / h[c1];
}

////////////////////////////////////////////////////////////////////////////////
// Stages 1 and 2: fine edges owned by a coarse edge or by a coarse face
//
// One template, nine instantiations (3 directions x {(0,0),(1,0),(0,1)}).
// (P1,P2) = (0,0) is Stage 1, pure injection.  The other two are Stage 2.
//
//     a_D = <A_D> + P1 * (h[B1]/16) * < B^{B2}(+e_B1) - B^{B2}(-e_B1) >
//                 - P2 * (h[B2]/16) * < B^{B1}(+e_B2) - B^{B1}(-e_B2) >
//
// where <.> on the A term averages over every transverse direction whose
// parity is 1, and <.> on each B term averages over the *other* transverse
// direction if its parity is 1.  That second average is what turns the 1/16
// into the 1/32 that appears when both parities are set.

template <int D, int P1, int P2, typename AA, typename HH>
CARPETX_VECPOT_HD inline auto surface_value(const AA &ac, const idx3 &x,
                                            const int n, const HH &hc) {
  static_assert(D >= 0 && D < 3, "");
  static_assert(P1 == 0 || P1 == 1, "");
  static_assert(P2 == 0 || P2 == 1, "");
  constexpr int b1 = (D + 1) % 3;
  constexpr int b2 = (D + 2) % 3;

  using T = decltype(at(ac[0], x, n));

  // <A_D>: mean over the 2^(P1+P2) coarse edges the fine edge sits between.
  T sum_a = T(0);
  for (int d1 = 0; d1 <= P1; ++d1)
    for (int d2 = 0; d2 <= P2; ++d2)
      sum_a += at(ac[D], shifted(shifted(x, b1, d1), b2, d2), n);
  T res = sum_a / T((P1 + 1) * (P2 + 1));

  // + P1 * (h[B1]/16) * < B^{B2}(+e_B1) - B^{B2}(-e_B1) >, averaged over B2.
  if (P1) {
    T s = T(0);
    for (int d2 = 0; d2 <= P2; ++d2) {
      const idx3 y = shifted(x, b2, d2);
      s += curl<b2>(ac, shifted(y, b1, +1), n, hc) -
           curl<b2>(ac, shifted(y, b1, -1), n, hc);
    }
    res += hc[b1] / T(16) * (s / T(P2 + 1));
  }

  // - P2 * (h[B2]/16) * < B^{B1}(+e_B2) - B^{B1}(-e_B2) >, averaged over B1.
  if (P2) {
    T s = T(0);
    for (int d1 = 0; d1 <= P1; ++d1) {
      const idx3 y = shifted(x, b1, d1);
      s += curl<b1>(ac, shifted(y, b2, +1), n, hc) -
           curl<b1>(ac, shifted(y, b2, -1), n, hc);
    }
    res -= hc[b2] / T(16) * (s / T(P1 + 1));
  }

  return res;
}

/// Fine index of the first of the two fine edges owned by coarse index `x`.
template <int D, int P1, int P2>
CARPETX_VECPOT_HD constexpr idx3 surface_fine_index(const idx3 &x) {
  constexpr int b1 = (D + 1) % 3;
  constexpr int b2 = (D + 2) % 3;
  idx3 f{{0, 0, 0}};
  f[D] = 2 * x[D];
  f[b1] = 2 * x[b1] + P1;
  f[b2] = 2 * x[b2] + P2;
  return f;
}

/// Write the two fine A_D edges owned by coarse index `x`, honouring `mask`.
///
/// The equal-split gauge of Stage 1 gives both halves the same value, so the
/// value is computed once.  `mask` may be a null accessor; the guard mirrors
/// amrex::facediv_face_interp.
template <int D, int P1, int P2, typename AA, typename AF, typename MK,
          typename HH>
CARPETX_VECPOT_HD inline void surface_interp(const AA &ac, const AF &af,
                                             const MK &mask, const idx3 &x,
                                             const int n, const HH &hc) {
  if (mask)
    if (!mask(x[0], x[1], x[2], n))
      return;
  const auto v = surface_value<D, P1, P2>(ac, x, n, hc);
  idx3 f = surface_fine_index<D, P1, P2>(x);
  af[D](f[0], f[1], f[2], n) = v;
  f[D] += 1;
  af[D](f[0], f[1], f[2], n) = v;
}

/// As `surface_interp`, but with the triple (D, P1, P2) packed into a single
/// template argument as D*4 + P1*2 + P2.
///
/// Launch macros such as AMREX_LOOP_3D and AMREX_HOST_DEVICE_PARALLEL_FOR_4D
/// take a fixed number of arguments and split their body on top-level commas.
/// Angle brackets do not group, so `surface_interp<D, P1, P2>(...)` inside such
/// a macro is split into three arguments and fails to compile.  Braces do not
/// group either, which is why idx3 must be built inside the call parentheses.
template <int DP, typename AA, typename AF, typename MK, typename HH>
CARPETX_VECPOT_HD inline void surface_interp_n(const AA &ac, const AF &af,
                                               const MK &mask, const idx3 &x,
                                               const int n, const HH &hc) {
  surface_interp<DP / 4, (DP / 2) % 2, DP % 2>(ac, af, mask, x, n, hc);
}

////////////////////////////////////////////////////////////////////////////////
// Stage 3: fine edges in the interior of a coarse cell
//
// Provenance-agnostic: this reads ONLY the fine (destination) array, so it is
// correct whether the surrounding surface edges were prolonged by the stencil
// above or copied from an already refined neighbour.  That is the draft's
// Eq. (33) moment recovery, and it is what makes the copy rule work without a
// branch.  It also matches the signature of amrex::facediv_int, which likewise
// takes no coarse array.
//
//     a_D = <A_D>_fine
//         + (h[B1]/16) * sum_{s,e,f} (2f-1) B^{B2}_fine(...)
//         - (h[B2]/16) * sum_{s,e,f} (2f-1) B^{B1}_fine(...)
//
// with h the COARSE spacing; the fine curls use h/2.

template <int D, typename AF, typename HH>
CARPETX_VECPOT_HD inline auto interior_value(const AF &af, const idx3 &x,
                                             const int n, const HH &hc) {
  static_assert(D >= 0 && D < 3, "");
  constexpr int b1 = (D + 1) % 3;
  constexpr int b2 = (D + 2) % 3;

  using T = decltype(at(af[0], x, n));

  T hf[3];
  for (int d = 0; d < 3; ++d)
    hf[d] = hc[d] / T(2);

  // <A_D>: mean of the eight fine halves on the four bounding coarse edges.
  // By (C2) this equals the mean of the four coarse edge values, but reading
  // the fine array instead keeps the formula independent of whether those
  // edges were copied (PROLONG_A.md section 2.5).
  T sum_a = T(0);
  for (int h = 0; h < 2; ++h)
    for (int d1 = 0; d1 < 2; ++d1)
      for (int d2 = 0; d2 < 2; ++d2) {
        idx3 f{{0, 0, 0}};
        f[D] = 2 * x[D] + h;
        f[b1] = 2 * (x[b1] + d1);
        f[b2] = 2 * (x[b2] + d2);
        sum_a += at(af[D], f, n);
      }
  T res = sum_a / T(8);

  // + (h[B1]/16) * moment of B^{B2} over this cell's two B2-faces.
  {
    T s = T(0);
    for (int sb = 0; sb < 2; ++sb)
      for (int e = 0; e < 2; ++e)
        for (int f1 = 0; f1 < 2; ++f1) {
          idx3 f{{0, 0, 0}};
          f[D] = 2 * x[D] + e;
          f[b1] = 2 * x[b1] + f1;
          f[b2] = 2 * (x[b2] + sb);
          s += T(2 * f1 - 1) * curl<b2>(af, f, n, hf);
        }
    res += hc[b1] / T(16) * s;
  }

  // - (h[B2]/16) * moment of B^{B1} over this cell's two B1-faces.
  {
    T s = T(0);
    for (int sb = 0; sb < 2; ++sb)
      for (int e = 0; e < 2; ++e)
        for (int f2 = 0; f2 < 2; ++f2) {
          idx3 f{{0, 0, 0}};
          f[D] = 2 * x[D] + e;
          f[b1] = 2 * (x[b1] + sb);
          f[b2] = 2 * x[b2] + f2;
          s += T(2 * f2 - 1) * curl<b1>(af, f, n, hf);
        }
    res -= hc[b2] / T(16) * s;
  }

  return res;
}

/// Write the two fine A_D edges in the interior of coarse cell `x`.
template <int D, typename AF, typename MK, typename HH>
CARPETX_VECPOT_HD inline void interior_interp(const AF &af, const MK &mask,
                                              const idx3 &x, const int n,
                                              const HH &hc) {
  if (mask)
    if (!mask(x[0], x[1], x[2], n))
      return;
  const auto v = interior_value<D>(af, x, n, hc);
  constexpr int b1 = (D + 1) % 3;
  constexpr int b2 = (D + 2) % 3;
  idx3 f{{0, 0, 0}};
  f[D] = 2 * x[D];
  f[b1] = 2 * x[b1] + 1;
  f[b2] = 2 * x[b2] + 1;
  af[D](f[0], f[1], f[2], n) = v;
  f[D] += 1;
  af[D](f[0], f[1], f[2], n) = v;
}

} // namespace vecpot
} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_VECPOT_KERNELS_HXX
