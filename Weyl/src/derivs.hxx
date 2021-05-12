#ifndef DERIVS_HXX
#define DERIVS_HXX

#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace Weyl {
using namespace Arith;
using namespace Loop;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

constexpr int deriv_order = 4;

////////////////////////////////////////////////////////////////////////////////

#if 0
template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T
deriv1d(const T *restrict const var, const ptrdiff_t di, const T dx) {
  if constexpr (deriv_order == 2)
    return -1 / T(2) * (var[-di] - var[+di]) / dx;
  if constexpr (deriv_order == 4)
    return (1 / T(12) * (var[-2 * di] - var[+2 * di]) //
            - 2 / T(3) * (var[-di] - var[+di])) /
           dx;
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T deriv1d_upwind(
    const T *restrict const var, const ptrdiff_t di, const T sign, const T dx) {
  // arXiv:1111.2177 [gr-qc], (71)
  if constexpr (deriv_order == 2) {
    // if (sign)
    //   // +     [ 0   -1   +1    0    0]
    //   // + 1/2 [+1   -2   +1    0    0]
    //   //       [+1/2 -2   +3/2  0    0]
    //   return (1 / T(2) * var[-2 * di] //
    //           - 2 * var[-di]          //
    //           + 3 / T(2) * var[0]) /
    //          dx;
    // else
    //   // +     [ 0    0   -1   +1    0  ]
    //   // - 1/2 [ 0    0   +1   -2   +1  ]
    //   //       [ 0    0   -3/2 +2   -1/2]
    //   return (-3 / T(2) * var[0] //
    //           + 2 * var[+di]     //
    //           - 1 / T(2) * var[+2 * di]) /
    //          dx;

    // + 1/2 [+1/2 -2   +3/2  0    0  ]
    // + 1/2 [ 0    0   -3/2 +2   -1/2]
    //       [+1/4 -1    0   +1   -1/4]
    const T symm = 1 / T(4) * (var[-2 * di] - var[2 * di]) //
                   - (var[-di] - var[di]);
    // + 1/2 [+1/2 -2   +3/2  0    0  ]
    // - 1/2 [ 0    0   -3/2 +2   -1/2]
    //       [+1/4 -1   +3/2 -1   +1/4]
    const T anti = 1 / T(4) * (var[-2 * di] + var[2 * di]) //
                   - (var[-di] + var[di])                  //
                   + 3 / T(2) * var[0];
    return (symm - flipsign(anti, sign)) / dx;
  }
  if constexpr (deriv_order == 4) {
    // A fourth order stencil for a first derivative, shifted by one grid point

    // if (sign)
    //   return (-1 / T(12) * var[-3 * di] //
    //           + 1 / T(2) * var[-2 * di] //
    //           - 3 / T(2) * var[-di]     //
    //           + 5 / T(6) * var[0]       //
    //           + 1 / T(4) * var[+di]) /
    //          dx;
    // else
    //   return (-1 / T(4) * var[-di]      //
    //           - 5 / T(6) * var[0]       //
    //           + 3 / T(2) * var[+di]     //
    //           - 1 / T(2) * var[+2 * di] //
    //           + 1 / T(12) * var[+3 * di]) /
    //          dx;

    const T symm = -1 / T(24) * (var[-3 * di] - var[3 * di]) //
                   + 1 / T(4) * (var[-2 * di] - var[2 * di]) //
                   - 7 / T(8) * (var[-di] - var[di]);
    const T anti = -1 / T(24) * (var[-3 * di] + var[3 * di]) //
                   + 1 / T(4) * (var[-2 * di] + var[2 * di]) //
                   - 5 / T(8) * (var[-di] + var[di])         //
                   + 5 / T(6) * var[0];
    return (symm - flipsign(anti, sign)) / dx;
  }
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T
deriv2_1d(const T *restrict const var, const ptrdiff_t di, const T dx) {
  if constexpr (deriv_order == 2)
    return ((var[-di] + var[+di]) //
            - 2 * var[0]) /
           pow2(dx);
  if constexpr (deriv_order == 4)
    return (-1 / T(12) * (var[-2 * di] + var[+2 * di]) //
            + 4 / T(3) * (var[-di] + var[+di])         //
            - 5 / T(2) * var[0]) /
           pow2(dx);
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T
deriv2_2d(const T *restrict const var, const ptrdiff_t di, const ptrdiff_t dj,
          const T dx, const T dy) {
  array<T, deriv_order + 1> arrx;
  T *const varx = &arrx[deriv_order / 2];
  for (int j = -deriv_order / 2; j <= deriv_order / 2; ++j)
    varx[j] = deriv1d(&var[j * dj], di, dx);
  return deriv1d(varx, 1, dy);
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T
deriv1d_diss(const T *restrict const var, const ptrdiff_t di, const T dx) {
  if constexpr (deriv_order == 2)
    return ((var[-2 * di] + var[+2 * di]) //
            - 4 * (var[-di] + var[+di])   //
            + 6 * var[0]) /
           dx;
  if constexpr (deriv_order == 4)
    return ((var[-3 * di] + var[+3 * di])       //
            - 6 * (var[-2 * di] + var[+2 * di]) //
            + 15 * (var[-di] + var[+di])        //
            - 20 * var[0]) /
           dx;
}
#endif

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv1d(const simdl<T> &mask, const T *restrict const var, const ptrdiff_t di,
        const T dx) {
  if constexpr (deriv_order == 2)
    return -1 / T(2) *
           (maskz_loadu(mask, &var[-di]) - maskz_loadu(mask, &var[+di])) / dx;
  if constexpr (deriv_order == 4)
    return (1 / T(12) *
                (maskz_loadu(mask, &var[-2 * di]) -
                 maskz_loadu(mask, &var[+2 * di])) //
            -
            2 / T(3) *
                (maskz_loadu(mask, &var[-di]) - maskz_loadu(mask, &var[+di]))) /
           dx;
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv1d_upwind(const simdl<T> &mask, const T *restrict const var,
               const ptrdiff_t di, const simd<T> &vel, const T dx) {
  // arXiv:1111.2177 [gr-qc], (71)
  if constexpr (deriv_order == 2) {
    // if (sign)
    //   // +     [ 0   -1   +1    0    0]
    //   // + 1/2 [+1   -2   +1    0    0]
    //   //       [+1/2 -2   +3/2  0    0]
    //   return (1 / T(2) * var[-2 * di] //
    //           - 2 * var[-di]          //
    //           + 3 / T(2) * var[0]) /
    //          dx;
    // else
    //   // +     [ 0    0   -1   +1    0  ]
    //   // - 1/2 [ 0    0   +1   -2   +1  ]
    //   //       [ 0    0   -3/2 +2   -1/2]
    //   return (-3 / T(2) * var[0] //
    //           + 2 * var[+di]     //
    //           - 1 / T(2) * var[+2 * di]) /
    //          dx;

    // + 1/2 [+1/2 -2   +3/2  0    0  ]
    // + 1/2 [ 0    0   -3/2 +2   -1/2]
    //       [+1/4 -1    0   +1   -1/4]
    const simd<T> symm =
        1 / T(4) *
            (maskz_loadu(mask, &var[-2 * di]) -
             maskz_loadu(mask, &var[2 * di])) //
        - (maskz_loadu(mask, &var[-di]) - maskz_loadu(mask, &var[di]));
    // + 1/2 [+1/2 -2   +3/2  0    0  ]
    // - 1/2 [ 0    0   -3/2 +2   -1/2]
    //       [+1/4 -1   +3/2 -1   +1/4]
    const simd<T> anti =
        1 / T(4) *
            (maskz_loadu(mask, &var[-2 * di]) +
             maskz_loadu(mask, &var[2 * di]))                          //
        - (maskz_loadu(mask, &var[-di]) + maskz_loadu(mask, &var[di])) //
        + 3 / T(2) * maskz_loadu(mask, &var[0]);
    return (vel * symm - fabs(vel) * anti) / dx;
  }
  if constexpr (deriv_order == 4) {
    // A fourth order stencil for a first derivative, shifted by one grid point

    // if (sign)
    //   return (-1 / T(12) * var[-3 * di] //
    //           + 1 / T(2) * var[-2 * di] //
    //           - 3 / T(2) * var[-di]     //
    //           + 5 / T(6) * var[0]       //
    //           + 1 / T(4) * var[+di]) /
    //          dx;
    // else
    //   return (-1 / T(4) * var[-di]      //
    //           - 5 / T(6) * var[0]       //
    //           + 3 / T(2) * var[+di]     //
    //           - 1 / T(2) * var[+2 * di] //
    //           + 1 / T(12) * var[+3 * di]) /
    //          dx;

    const simd<T> symm =
        -1 / T(24) *
            (maskz_loadu(mask, &var[-3 * di]) -
             maskz_loadu(mask, &var[3 * di])) //
        + 1 / T(4) *
              (maskz_loadu(mask, &var[-2 * di]) -
               maskz_loadu(mask, &var[2 * di])) //
        -
        7 / T(8) * (maskz_loadu(mask, &var[-di]) - maskz_loadu(mask, &var[di]));
    const simd<T> anti =
        -1 / T(24) *
            (maskz_loadu(mask, &var[-3 * di]) +
             maskz_loadu(mask, &var[3 * di])) //
        + 1 / T(4) *
              (maskz_loadu(mask, &var[-2 * di]) +
               maskz_loadu(mask, &var[2 * di])) //
        - 5 / T(8) *
              (maskz_loadu(mask, &var[-di]) + maskz_loadu(mask, &var[di])) //
        + 5 / T(6) * maskz_loadu(mask, &var[0]);
    return (symm - fabs(vel) * anti) / dx;
  }
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv2_1d(const simdl<T> &mask, const T *restrict const var, const ptrdiff_t di,
          const T dx) {
  if constexpr (deriv_order == 2)
    return ((maskz_loadu(mask, &var[-di]) + maskz_loadu(mask, &var[+di])) //
            - 2 * maskz_loadu(mask, &var[0])) /
           pow2(dx);
  if constexpr (deriv_order == 4)
    return (-1 / T(12) *
                (maskz_loadu(mask, &var[-2 * di]) +
                 maskz_loadu(mask, &var[+2 * di])) //
            +
            4 / T(3) *
                (maskz_loadu(mask, &var[-di]) + maskz_loadu(mask, &var[+di])) //
            - 5 / T(2) * maskz_loadu(mask, &var[0])) /
           pow2(dx);
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv2_2d(const int vavail, const simdl<T> &mask, const T *restrict const var,
          const ptrdiff_t di, const ptrdiff_t dj, const T dx, const T dy) {
  constexpr size_t vsize = tuple_size_v<simd<T> >;
  if (di == 1) {
    assert(vavail > 0);
    constexpr int maxnpoints = deriv_order + 1 + vsize - 1;
    const int npoints = deriv_order + 1 + min(int(vsize), vavail) - 1;
    alignas(alignof(simd<T>)) array<T, maxnpoints> arrx;
    for (int i = 0; i < npoints; i += vsize) {
      const simdl<T> mask1 = mask_for_loop_tail<simdl<T> >(i, npoints);
      mask_storea(mask1, &arrx[i], deriv1d(mask1, &var[i], dj, dx));
    }
    const T *const varx = (T *)(&arrx[deriv_order / 2]);
    return deriv1d(mask, varx, 1, dy);
  } else {
    assert(dj != 1);
    array<simd<T>, deriv_order + 1> arrx;
    for (int j = -deriv_order / 2; j <= deriv_order / 2; ++j)
      arrx[deriv_order / 2 + j] = deriv1d(mask, &var[j * dj], di, dx);
    const T *const varx = (T *)(&arrx[deriv_order / 2]);
    return deriv1d(mask, varx, vsize, dy);
  }
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv1d_diss(const simdl<T> &mask, const T *restrict const var,
             const ptrdiff_t di, const T dx) {
  if constexpr (deriv_order == 2)
    return ((maskz_loadu(mask, &var[-2 * di]) +
             maskz_loadu(mask, &var[+2 * di])) //
            -
            4 * (maskz_loadu(mask, &var[-di]) + maskz_loadu(mask, &var[+di])) //
            + 6 * maskz_loadu(mask, &var[0])) /
           dx;
  if constexpr (deriv_order == 4)
    return ((maskz_loadu(mask, &var[-3 * di]) +
             maskz_loadu(mask, &var[+3 * di])) //
            - 6 * (maskz_loadu(mask, &var[-2 * di]) +
                   maskz_loadu(mask, &var[+2 * di])) //
            + 15 * (maskz_loadu(mask, &var[-di]) +
                    maskz_loadu(mask, &var[+di])) //
            - 20 * maskz_loadu(mask, &var[0])) /
           dx;
}

////////////////////////////////////////////////////////////////////////////////

#if 0
template <int dir, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T deriv(const GF3D2<const T> &gf_,
                                  const vect<int, dim> &I,
                                  const vec<T, D, UP> &dx) {
  static_assert(dir >= 0 && dir < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir));
  return deriv1d(&gf_(I), di, dx(dir));
}

template <int dir, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T deriv_upwind(const GF3D2<const T> &gf_,
                                         const vect<int, dim> &I,
                                         const vec<T, D, UP> &sign,
                                         const vec<T, D, UP> &dx) {
  static_assert(dir >= 0 && dir < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir));
  return deriv1d_upwind(&gf_(I), di, sign(dir), dx(dir));
}

template <int dir1, int dir2, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST enable_if_t<(dir1 == dir2), T>
deriv2(const GF3D2<const T> &gf_, const vect<int, dim> &I,
       const vec<T, D, UP> &dx) {
  static_assert(dir1 >= 0 && dir1 < D, "");
  static_assert(dir2 >= 0 && dir2 < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir1));
  return deriv2_1d(&gf_(I), di, dx(dir1));
}

template <int dir1, int dir2, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST enable_if_t<(dir1 != dir2), T>
deriv2(const GF3D2<const T> &gf_, const vect<int, dim> &I,
       const vec<T, D, UP> &dx) {
  static_assert(dir1 >= 0 && dir1 < D, "");
  static_assert(dir2 >= 0 && dir2 < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir1));
  const ptrdiff_t dj = gf_.delta(DI(dir2));
  return deriv2_2d(&gf_(I), di, dj, dx(dir1), dx(dir2));
}

template <int dir, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T deriv_diss(const GF3D2<const T> &gf_,
                                       const vect<int, dim> &I,
                                       const vec<T, D, UP> &dx) {
  static_assert(dir >= 0 && dir < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir));
  return deriv1d_diss(&gf_(I), di, dx(dir));
}
#endif

////////////////////////////////////////////////////////////////////////////////

template <int dir, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv(const simdl<T> &mask, const GF3D2<const T> &gf_, const vect<int, dim> &I,
      const vec<T, D, UP> &dx) {
  static_assert(dir >= 0 && dir < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir));
  return deriv1d(mask, &gf_(I), di, dx(dir));
}

template <int dir, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv_upwind(const simdl<T> &mask, const GF3D2<const T> &gf_,
             const vect<int, dim> &I, const vec<simd<T>, D, UP> &vel,
             const vec<T, D, UP> &dx) {
  static_assert(dir >= 0 && dir < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir));
  return deriv1d_upwind(mask, &gf_(I), di, vel(dir), dx(dir));
}

template <int dir1, int dir2, typename T, int D>
inline ARITH_INLINE
    ARITH_DEVICE ARITH_HOST enable_if_t<(dir1 == dir2), simd<T> >
    deriv2(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
           const vect<int, dim> &I, const vec<T, D, UP> &dx) {
  static_assert(dir1 >= 0 && dir1 < D, "");
  static_assert(dir2 >= 0 && dir2 < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir1));
  return deriv2_1d(mask, &gf_(I), di, dx(dir1));
}

template <int dir1, int dir2, typename T, int D>
inline ARITH_INLINE
    ARITH_DEVICE ARITH_HOST enable_if_t<(dir1 != dir2), simd<T> >
    deriv2(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
           const vect<int, dim> &I, const vec<T, D, UP> &dx) {
  static_assert(dir1 >= 0 && dir1 < D, "");
  static_assert(dir2 >= 0 && dir2 < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir1));
  const ptrdiff_t dj = gf_.delta(DI(dir2));
  return deriv2_2d(vavail, mask, &gf_(I), di, dj, dx(dir1), dx(dir2));
}

template <int dir, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv_diss(const simdl<T> &mask, const GF3D2<const T> &gf_,
           const vect<int, dim> &I, const vec<T, D, UP> &dx) {
  static_assert(dir >= 0 && dir < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir));
  return deriv1d_diss(mask, &gf_(I), di, dx(dir));
}

////////////////////////////////////////////////////////////////////////////////

#if 0
template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, dim, DN> deriv(const GF3D2<const T> &gf_,
                                                const vect<int, dim> &I,
                                                const vec<T, dim, UP> &dx) {
  return vec<T, dim, DN>(deriv<0>(gf_, I, dx), deriv<1>(gf_, I, dx),
                         deriv<2>(gf_, I, dx));
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, dim, DN>
deriv_upwind(const GF3D2<const T> &gf_, const vect<int, dim> &I,
             const vec<T, dim, UP> &dir, const vec<T, dim, UP> &dx) {
  return vec<T, D, DN>(deriv_upwind<0>(gf_, I, dir, dx),
                       deriv_upwind<1>(gf_, I, dir, dx),
                       deriv_upwind<2>(gf_, I, dir, dx));
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST smat<T, dim, DN, DN>
deriv2(const GF3D2<const T> &gf_, const vect<int, dim> &I,
       const vec<T, dim, UP> &dx) {
  return smat<T, dim, DN, DN>(
      deriv2<0, 0>(gf_, I, dx), deriv2<0, 1>(gf_, I, dx),
      deriv2<0, 2>(gf_, I, dx), deriv2<1, 1>(gf_, I, dx),
      deriv2<1, 2>(gf_, I, dx), deriv2<2, 2>(gf_, I, dx));
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T diss(const GF3D2<const T> &gf_,
                                 const vect<int, dim> &I,
                                 const vec<T, dim, UP> &dx) {
  // arXiv:gr-qc/0610128, (63), with r=2
  constexpr int diss_order = deriv_order + 2;
  constexpr int sign = diss_order % 4 == 0 ? -1 : +1;
  return sign / T(pown(2, deriv_order + 2)) *
         (deriv_diss<0>(gf_, I, dx)   //
          + deriv_diss<1>(gf_, I, dx) //
          + deriv_diss<2>(gf_, I, dx));
}
#endif

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<simd<T>, dim, DN>
deriv(const simdl<T> &mask, const GF3D2<const T> &gf_, const vect<int, dim> &I,
      const vec<T, dim, UP> &dx) {
  return {deriv<0>(mask, gf_, I, dx), deriv<1>(mask, gf_, I, dx),
          deriv<2>(mask, gf_, I, dx)};
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv_upwind(const simdl<T> &mask, const GF3D2<const T> &gf_,
             const vect<int, dim> &I, const vec<simd<T>, dim, UP> &vel,
             const vec<T, dim, UP> &dx) {
  return deriv_upwind<0>(mask, gf_, I, vel, dx) +
         deriv_upwind<1>(mask, gf_, I, vel, dx) +
         deriv_upwind<2>(mask, gf_, I, vel, dx);
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST smat<simd<T>, dim, DN, DN>
deriv2(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
       const vect<int, dim> &I, const vec<T, dim, UP> &dx) {
  return {deriv2<0, 0>(vavail, mask, gf_, I, dx),
          deriv2<0, 1>(vavail, mask, gf_, I, dx),
          deriv2<0, 2>(vavail, mask, gf_, I, dx),
          deriv2<1, 1>(vavail, mask, gf_, I, dx),
          deriv2<1, 2>(vavail, mask, gf_, I, dx),
          deriv2<2, 2>(vavail, mask, gf_, I, dx)};
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
diss(const simdl<T> &mask, const GF3D2<const T> &gf_, const vect<int, dim> &I,
     const vec<T, dim, UP> &dx) {
  // arXiv:gr-qc/0610128, (63), with r=2
  constexpr int diss_order = deriv_order + 2;
  constexpr int sign = diss_order % 4 == 0 ? -1 : +1;
  return sign / T(pown(2, deriv_order + 2)) *
         (deriv_diss<0>(mask, gf_, I, dx)   //
          + deriv_diss<1>(mask, gf_, I, dx) //
          + deriv_diss<2>(mask, gf_, I, dx));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH, const GF3D2<const T> &gf1,
            const GF3D5<T> &gf0, const vec<GF3D5<T>, dim, DN> &dgf0,
            const GF3D5layout &layout0) {
  DECLARE_CCTK_ARGUMENTS;

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const vec<CCTK_REAL, dim, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);
        const auto val = gf1(mask, p.I);
        gf0.store(mask, index0, val);
        const auto dval = deriv(mask, gf1, p.I, dx);
        dgf0.store(mask, index0, dval);
      });
}

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH, const GF3D2<const T> &gf1,
             const GF3D5<T> &gf0, const vec<GF3D5<T>, dim, DN> &dgf0,
             const smat<GF3D5<T>, dim, DN, DN> &ddgf0,
             const GF3D5layout &layout0) {
  DECLARE_CCTK_ARGUMENTS;

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const vec<CCTK_REAL, dim, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);
        const auto val = gf1(mask, p.I);
        gf0.store(mask, index0, val);
        const auto dval = deriv(mask, gf1, p.I, dx);
        dgf0.store(mask, index0, dval);
        const auto ddval = deriv2(p.imax - p.i, mask, gf1, p.I, dx);
        ddgf0.store(mask, index0, ddval);
      });
}

template <typename T, dnup_t dnup>
void calc_derivs(const cGH *restrict const cctkGH,
                 const vec<GF3D2<const T>, dim, dnup> &gf0_,
                 const vec<GF3D5<T>, dim, dnup> &gf_,
                 const vec<vec<GF3D5<T>, dim, DN>, dim, dnup> &dgf_,
                 const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    calc_derivs(cctkGH, gf0_(a), gf_(a), dgf_(a), layout);
}

template <typename T, dnup_t dnup>
void calc_derivs2(const cGH *restrict const cctkGH,
                  const vec<GF3D2<const T>, dim, dnup> &gf0_,
                  const vec<GF3D5<T>, dim, dnup> &gf_,
                  const vec<vec<GF3D5<T>, dim, DN>, dim, dnup> &dgf_,
                  const vec<smat<GF3D5<T>, dim, DN, DN>, dim, dnup> &ddgf_,
                  const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    calc_derivs2(cctkGH, gf0_(a), gf_(a), dgf_(a), ddgf_(a), layout);
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
void calc_derivs(const cGH *restrict const cctkGH,
                 const smat<GF3D2<const T>, dim, dnup1, dnup2> &gf0_,
                 const smat<GF3D5<T>, dim, dnup1, dnup2> &gf_,
                 const smat<vec<GF3D5<T>, dim, DN>, dim, dnup1, dnup2> &dgf_,
                 const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b), layout);
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
void calc_derivs2(
    const cGH *restrict const cctkGH,
    const smat<GF3D2<const T>, dim, dnup1, dnup2> &gf0_,
    const smat<GF3D5<T>, dim, dnup1, dnup2> &gf_,
    const smat<vec<GF3D5<T>, dim, DN>, dim, dnup1, dnup2> &dgf_,
    const smat<smat<GF3D5<T>, dim, DN, DN>, dim, dnup1, dnup2> &ddgf_,
    const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs2(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b), ddgf_(a, b),
                   layout);
}

} // namespace Z4c

#endif // #ifndef DERIVS_HXX
