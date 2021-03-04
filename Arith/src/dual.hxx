#ifndef DUAL_HXX
#define DUAL_HXX

#ifdef __CUDACC__
#define DUAL_DEVICE __device__
#define DUAL_HOST __host__
#else
#define DUAL_DEVICE
#define DUAL_HOST
#endif

#include "vect.hxx"

#include <cassert>
#include <cmath>
#include <functional>
#include <ostream>

namespace Arith {
using namespace std;

template <typename T, typename U = T> struct dual;

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
abs(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
cbrt(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
cos(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
exp(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
fabs(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
isnan(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
pow(const Arith::dual<T, U> &x, int n);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
pow2(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
sin(const Arith::dual<T, U> &x);
template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
sqrt(const Arith::dual<T, U> &x);

template <typename T, typename U> struct dual {
  T val;
  U eps;

  CCTK_ATTRIBUTE_ALWAYS_INLINE dual(const dual &) = default;
  CCTK_ATTRIBUTE_ALWAYS_INLINE dual(dual &&) = default;
  CCTK_ATTRIBUTE_ALWAYS_INLINE dual &operator=(const dual &) = default;
  CCTK_ATTRIBUTE_ALWAYS_INLINE dual &operator=(dual &&) = default;

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual()
      : val(), eps() {}
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual(const T &x)
      : val(x), eps() {}
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual(const T &x,
                                                                    const U &y)
      : val(x), eps(y) {}

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator+(const dual &x) {
    return {+x.val, +x.eps};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator-(const dual &x) {
    return {-x.val, -x.eps};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator+(const dual &x, const dual &y) {
    return {x.val + y.val, x.eps + y.eps};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator-(const dual &x, const dual &y) {
    return {x.val - y.val, x.eps - y.eps};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator+(const dual &x, const T &y) {
    return {x.val + y, x.eps};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator-(const dual &x, const T &y) {
    return {x.val - y, x.eps};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator*(const dual &x, const T &y) {
    return {x.val * y, x.eps * y};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator*(const T &x, const dual &y) {
    return {x * y.val, x * y.eps};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator/(const dual &x, const T &y) {
    return {x.val / y, x.eps / y};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator*(const dual &x, const dual &y) {
    return {x.val * y.val, x.val * y.eps + x.eps * y.val};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual
  operator/(const dual &x, const dual &y) {
    return {x.val / y.val, (x.eps * y.val - x.val * y.eps) / (y.val * y.val)};
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual &
  operator+=(const dual &x) {
    return *this = *this + x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual &
  operator-=(const dual &x) {
    return *this = *this - x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual &
  operator*=(const T &x) {
    return *this = *this * x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual &
  operator/=(const T &x) {
    return *this = *this / x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual &
  operator*=(const dual &x) {
    return *this = *this * x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual &
  operator/=(const dual &x) {
    return *this = *this / x;
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator==(const dual &x, const dual &y) {
    return x.val == y.val;
  };
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator<(const dual &x, const dual &y) {
    return x.val < y.val;
  };

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator!=(const dual &x, const dual &y) {
    return !(x == y);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator>(const dual &x, const dual &y) {
    return y < x;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator<=(const dual &x, const dual &y) {
    return !(x > y);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator>=(const dual &x, const dual &y) {
    return !(x < y);
  }

  friend ostream &operator<<(ostream &os, const dual &x) {
    return os << x.val << "+eps*" << x.val;
  }
};

template <typename T, typename U> struct zero<dual<T, U> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual<T, U>
  operator()() const {
    return dual<T, U>(zero<T>()(), zero<U>()());
  }
};

} // namespace Arith
namespace std {
template <typename T, typename U> struct equal_to<Arith::dual<T, U> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator()(const Arith::dual<T, U> &x, const Arith::dual<T, U> &y) const {
    return equal_to<T>()(x.val, y.val) && equal_to<U>()(x.eps, y.eps);
  }
};

template <typename T, typename U> struct less<Arith::dual<T, U> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
  operator()(const Arith::dual<T, U> &x, const Arith::dual<T, U> &y) const {
    if (less<T>(x.val, y.val))
      return true;
    if (less<T>(y.val, x.val))
      return false;
    if (less<U>(x.eps, y.eps))
      return true;
    if (less<U>(y.eps, x.eps))
      return false;
    return false;
  }
};
} // namespace std

namespace Arith {

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
abs(const Arith::dual<T, U> &x) {
  // return sqrt(pow2(x));
  using std::abs, std::copysign;
  return {abs(x.val), copysign(T{1}, x.val) * x.eps};
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
cbrt(const Arith::dual<T, U> &x) {
  using std::cbrt;
  const T r = cbrt(x.val);
  return {r, r / (3 * x.val) * x.eps};
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
cos(const Arith::dual<T, U> &x) {
  using std::cos, std::sin;
  return {cos(x.val), -sin(x.val) * x.eps};
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
exp(const Arith::dual<T, U> &x) {
  using std::exp;
  const T r = exp(x.val);
  return {r, r * x.eps};
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
fabs(const Arith::dual<T, U> &x) {
  using std::copysign, std::fabs;
  // return sqrt(pow2(x));
  return {fabs(x.val), copysign(T{1}, x.val) * x.eps};
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST bool
isnan(const Arith::dual<T, U> &x) {
  using std::isnan;
  return isnan(x.val) || x.eps.isnan().any();
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
pow(const Arith::dual<T, U> &x, const int n) {
  using std::pow;
  if (n == 0)
    return {1, U{}};
  return {pow(x.val, n), n * pow(x.val, n - 1) * x.eps};
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
pow2(const Arith::dual<T, U> &x) {
  return x * x;
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
sin(const Arith::dual<T, U> &x) {
  using std::cos, std::sin;
  return {sin(x.val), cos(x.val) * x.eps};
}

template <typename T, typename U>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST Arith::dual<T, U>
sqrt(const Arith::dual<T, U> &x) {
  using std::sqrt;
  const T r = sqrt(x.val);
  return {r, x.eps / (2 * r)};
}

////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U> struct one<dual<T, U> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE DUAL_DEVICE DUAL_HOST dual<T, U>
  operator()() const {
    return dual<T, U>(one<T>()(), zero<U>()());
  }
};

} // namespace Arith

#endif // #ifndef DUAL_HXX
