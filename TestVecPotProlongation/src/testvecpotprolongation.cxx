#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

// Integration test for prolongation_type="vecpot": the flux-preserving
// prolongation of an edge-staggered vector potential.  See PROLONG_A.md.
//
// The headline check is "pure gauge".  Initial data is the exact discrete
// gradient of a smooth scalar, A_i = (lambda(x + e_i h) - lambda(x)) / h, so
// the magnetic field is zero to machine precision on every level.  A
// prolongation that merely interpolates A_i would produce a nonzero B in the
// refinement-boundary ghost zones; this one must not.  Look at max |Bvec_*|.
//
// The "smooth" data has a nonzero B and is there for convergence and for flux
// conservation across the refinement boundary.
//
// divB is identically zero for any edge-staggered A, prolonged or not, so it
// checks this thorn's own bookkeeping rather than the operator.

namespace TestVecPotProlongation {

// The gauge scalar.  Only its values matter, never its derivatives: the
// initial data below differences it exactly.
template <typename T>
static inline CCTK_DEVICE CCTK_HOST T lambda(const T A, const T kx, const T ky,
                                             const T kz, const T x, const T y,
                                             const T z) {
  using std::cos, std::sin;
  return A * sin(kx * x) * cos(ky * y) * sin(kz * z);
}

// A smooth vector potential with a nonzero curl.
template <typename T>
static inline CCTK_DEVICE CCTK_HOST T smooth_Ax(const T A, const T kx,
                                                const T ky, const T kz,
                                                const T x, const T y,
                                                const T z) {
  using std::cos, std::sin;
  return A * sin(ky * y) * cos(kz * z);
}
template <typename T>
static inline CCTK_DEVICE CCTK_HOST T smooth_Ay(const T A, const T kx,
                                                const T ky, const T kz,
                                                const T x, const T y,
                                                const T z) {
  using std::cos, std::sin;
  return A * sin(kz * z) * cos(kx * x);
}
template <typename T>
static inline CCTK_DEVICE CCTK_HOST T smooth_Az(const T A, const T kx,
                                                const T ky, const T kz,
                                                const T x, const T y,
                                                const T z) {
  using std::cos, std::sin;
  return A * sin(kx * x) * cos(ky * y);
}

extern "C" void TestVecPotProlongation_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestVecPotProlongation_Initial;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL A = amplitude;
  const CCTK_REAL k0 = kx, k1 = ky, k2 = kz;
  const bool gauge = CCTK_EQUALS(initial_condition, "pure gauge");

  // A_d lives on the edge parallel to d, between the two vertices at
  // p.X -+ DI[d] DX[d] / 2.  Differencing lambda between exactly those two
  // vertices makes A the discrete gradient, not merely its approximation, so
  // the discrete curl vanishes identically rather than to truncation order.
  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (gauge) {
          const CCTK_REAL h = p.DX[0];
          const CCTK_REAL xm = p.X[0] - h / 2, xp = p.X[0] + h / 2;
          Avec_x(p.I) = (lambda(A, k0, k1, k2, xp, p.X[1], p.X[2]) -
                         lambda(A, k0, k1, k2, xm, p.X[1], p.X[2])) /
                        h;
        } else {
          Avec_x(p.I) = smooth_Ax(A, k0, k1, k2, p.X[0], p.X[1], p.X[2]);
        }
      });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (gauge) {
          const CCTK_REAL h = p.DX[1];
          const CCTK_REAL ym = p.X[1] - h / 2, yp = p.X[1] + h / 2;
          Avec_y(p.I) = (lambda(A, k0, k1, k2, p.X[0], yp, p.X[2]) -
                         lambda(A, k0, k1, k2, p.X[0], ym, p.X[2])) /
                        h;
        } else {
          Avec_y(p.I) = smooth_Ay(A, k0, k1, k2, p.X[0], p.X[1], p.X[2]);
        }
      });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (gauge) {
          const CCTK_REAL h = p.DX[2];
          const CCTK_REAL zm = p.X[2] - h / 2, zp = p.X[2] + h / 2;
          Avec_z(p.I) = (lambda(A, k0, k1, k2, p.X[0], p.X[1], zp) -
                         lambda(A, k0, k1, k2, p.X[0], p.X[1], zm)) /
                        h;
        } else {
          Avec_z(p.I) = smooth_Az(A, k0, k1, k2, p.X[0], p.X[1], p.X[2]);
        }
      });
}

extern "C" void TestVecPotProlongation_Curl(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestVecPotProlongation_Curl;

  // The same discrete Stokes relation the operator is built on:
  //   B^x = (dAz/dy - dAy/dz), and cyclically.
  grid.loop_int_device<0, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Bvec_x(p.I) = (Avec_z(p.I + p.DI[1]) - Avec_z(p.I)) / p.DX[1] -
                      (Avec_y(p.I + p.DI[2]) - Avec_y(p.I)) / p.DX[2];
      });

  grid.loop_int_device<1, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Bvec_y(p.I) = (Avec_x(p.I + p.DI[2]) - Avec_x(p.I)) / p.DX[2] -
                      (Avec_z(p.I + p.DI[0]) - Avec_z(p.I)) / p.DX[0];
      });

  grid.loop_int_device<1, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Bvec_z(p.I) = (Avec_y(p.I + p.DI[0]) - Avec_y(p.I)) / p.DX[0] -
                      (Avec_x(p.I + p.DI[1]) - Avec_x(p.I)) / p.DX[1];
      });
}

extern "C" void TestVecPotProlongation_Div(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestVecPotProlongation_Div;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        divB(p.I) = (Bvec_x(p.I + p.DI[0]) - Bvec_x(p.I)) / p.DX[0] +
                    (Bvec_y(p.I + p.DI[1]) - Bvec_y(p.I)) / p.DX[1] +
                    (Bvec_z(p.I + p.DI[2]) - Bvec_z(p.I)) / p.DX[2];
      });
}

} // namespace TestVecPotProlongation
