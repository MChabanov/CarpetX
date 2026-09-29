#ifndef CARPETX_TEST_UNIT_AMREX_SHIM_HXX
#define CARPETX_TEST_UNIT_AMREX_SHIM_HXX

// A minimal, *functional* stand-in for the slice of AMReX and Cactus that
// CarpetX/src/prolongate_3d_rf2_vecpot.cxx depends on, so that the real
// interpolater can be compiled and executed without an Einstein Toolkit build.
//
// It is deliberately small: its size is a statement of exactly how much AMReX
// surface the operator touches.  Semantics follow AMReX (index-type conversion
// grows/shrinks the high end; coarsening floors) and were read off
// amrex/Src/Base/AMReX_Box.H and AMReX_IndexType.H.
//
// This shim is for testing only and is never compiled into CarpetX.

#include <array>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>

#define AMREX_SPACEDIM 3
#define AMREX_GPU_DEVICE
#define AMREX_GPU_HOST_DEVICE
#define AMREX_FORCE_INLINE inline

namespace amrex {

using Real = double;

template <typename T, std::size_t N> using GpuArray = std::array<T, N>;
template <typename T, std::size_t N> using Array = std::array<T, N>;
template <typename T> using Vector = std::vector<T>;

enum class RunOn { Gpu, Cpu, Device, Host };

struct BCRec {
  int lo[3]{}, hi[3]{};
};

struct IntVect {
  int v[3]{};
  IntVect() = default;
  explicit IntVect(int a) : v{a, a, a} {}
  IntVect(int a, int b, int c) : v{a, b, c} {}
  int operator[](int d) const { return v[d]; }
  int &operator[](int d) { return v[d]; }
  static IntVect TheZeroVector() { return IntVect(0); }
  static IntVect TheUnitVector() { return IntVect(1); }
  static IntVect TheNodeVector() { return IntVect(1); }
  static IntVect TheDimensionVector(int d) {
    IntVect r(0);
    r[d] = 1;
    return r;
  }
  friend IntVect operator-(IntVect a, const IntVect &b) {
    for (int d = 0; d < 3; ++d)
      a[d] -= b[d];
    return a;
  }
  friend bool operator==(const IntVect &a, const IntVect &b) {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
  }
};

// 0 == CELL, 1 == NODE (AMReX convention).
struct IndexType {
  IntVect t;
  IndexType() = default;
  explicit IndexType(const IntVect &iv) : t(iv) {}
  bool cellCentered(int d) const { return t[d] == 0; }
  bool nodeCentered(int d) const { return t[d] == 1; }
  bool cellCentered() const {
    return t[0] == 0 && t[1] == 0 && t[2] == 0;
  }
  int ixType(int d) const { return t[d]; }
  static IndexType TheCellType() { return IndexType(IntVect(0)); }
  static IndexType TheNodeType() { return IndexType(IntVect(1)); }
  friend bool operator==(const IndexType &a, const IndexType &b) {
    return a.t == b.t;
  }
};

struct Box {
  IntVect lo, hi;
  IndexType typ;
  Box() = default;
  Box(const IntVect &l, const IntVect &h,
      const IndexType &t = IndexType::TheCellType())
      : lo(l), hi(h), typ(t) {}
  const IntVect &smallEnd() const { return lo; }
  const IntVect &bigEnd() const { return hi; }
  int smallEnd(int d) const { return lo[d]; }
  int bigEnd(int d) const { return hi[d]; }
  IndexType ixType() const { return typ; }
  Box &grow(int n) {
    for (int d = 0; d < 3; ++d) {
      lo[d] -= n;
      hi[d] += n;
    }
    return *this;
  }
  Box &grow(int d, int n) {
    lo[d] -= n;
    hi[d] += n;
    return *this;
  }
  Box &growHi(int d, int n) {
    hi[d] += n;
    return *this;
  }
  Box &convert(const IndexType &t) {
    for (int d = 0; d < 3; ++d) {
      if (t.nodeCentered(d) && typ.cellCentered(d))
        hi[d] += 1;
      else if (t.cellCentered(d) && typ.nodeCentered(d))
        hi[d] -= 1;
    }
    typ = t;
    return *this;
  }
};

inline Box convert(Box b, const IndexType &t) { return b.convert(t); }
inline Box enclosedCells(Box b) { return b.convert(IndexType::TheCellType()); }
inline Box surroundingNodes(Box b) {
  return b.convert(IndexType::TheNodeType());
}
inline int floordiv(int a, int r) { return a >= 0 ? a / r : -((-a + r - 1) / r); }
inline Box coarsen(const Box &b, const IntVect &r) {
  assert(b.ixType().cellCentered()); // all we need here
  IntVect l, h;
  for (int d = 0; d < 3; ++d) {
    l[d] = floordiv(b.lo[d], r[d]);
    h[d] = floordiv(b.hi[d], r[d]);
  }
  return Box(l, h, b.typ);
}
inline Box coarsen(const Box &b, int r) { return coarsen(b, IntVect(r)); }
inline Box refine(const Box &b, const IntVect &r) {
  assert(b.ixType().cellCentered());
  IntVect l, h;
  for (int d = 0; d < 3; ++d) {
    l[d] = b.lo[d] * r[d];
    h[d] = (b.hi[d] + 1) * r[d] - 1;
  }
  return Box(l, h, b.typ);
}
inline Box refine(const Box &b, int r) { return refine(b, IntVect(r)); }
inline Box grow(Box b, const IntVect &n) {
  for (int d = 0; d < 3; ++d) {
    b.lo[d] -= n[d];
    b.hi[d] += n[d];
  }
  return b;
}

template <typename T> struct Array4 {
  T *p = nullptr;
  int lo[3]{}, n[3]{};
  int begin_comp = 0;
  Array4() = default;
  Array4(T *pp, const int l[3], const int nn[3], int bc)
      : p(pp), begin_comp(bc) {
    for (int d = 0; d < 3; ++d) {
      lo[d] = l[d];
      n[d] = nn[d];
    }
  }
  explicit operator bool() const { return p != nullptr; }
  T &operator()(int i, int j, int k, int c = 0) const {
    const int a = i - lo[0], b = j - lo[1], e = k - lo[2];
    if (a < 0 || a >= n[0] || b < 0 || b >= n[1] || e < 0 || e >= n[2]) {
      std::fprintf(stderr, "Array4 index (%d,%d,%d) out of range\n", i, j, k);
      std::abort();
    }
    return p[((std::size_t(c + begin_comp) * n[2] + e) * n[1] + b) * n[0] + a];
  }
};

template <typename T> struct BaseFab {
  std::vector<T> data;
  Box bx;
  int ncomp = 1;
  BaseFab() = default;
  BaseFab(const Box &b, int nc) : bx(b), ncomp(nc) {
    data.assign(std::size_t(len(0)) * len(1) * len(2) * nc, T(0));
  }
  int len(int d) const { return bx.hi[d] - bx.lo[d] + 1; }
  T *dataPtr() { return data.data(); }
  const T *dataPtr() const { return data.data(); }
  std::size_t size() const { return data.size(); }
  const Box &box() const { return bx; }
  Array4<T> array(int c = 0) {
    const int l[3] = {bx.lo[0], bx.lo[1], bx.lo[2]};
    const int n[3] = {len(0), len(1), len(2)};
    return Array4<T>(data.data(), l, n, c);
  }
  Array4<const T> const_array(int c = 0) const {
    const int l[3] = {bx.lo[0], bx.lo[1], bx.lo[2]};
    const int n[3] = {len(0), len(1), len(2)};
    return Array4<const T>(data.data(), l, n, c);
  }
};

using FArrayBox = BaseFab<Real>;
using IArrayBox = BaseFab<int>;

struct Geometry {
  GpuArray<Real, 3> dx{{1, 1, 1}};
  GpuArray<Real, 3> CellSizeArray() const { return dx; }
};

class Interpolater {
public:
  virtual ~Interpolater() = default;
  virtual Box CoarseBox(const Box &fine, int ratio) = 0;
  virtual Box CoarseBox(const Box &fine, const IntVect &ratio) = 0;
  virtual void interp(const FArrayBox &crse, int crse_comp, FArrayBox &fine,
                      int fine_comp, int ncomp, const Box &fine_region,
                      const IntVect &ratio, const Geometry &crse_geom,
                      const Geometry &fine_geom, Vector<BCRec> const &bcr,
                      int actual_comp, int actual_state, RunOn runon) = 0;
  virtual void
  interp_arr(Array<FArrayBox *, 3> const &, int, Array<FArrayBox *, 3> const &,
             int, int, const Box &, const IntVect &,
             Array<IArrayBox *, 3> const &, const Geometry &, const Geometry &,
             Vector<Array<BCRec, 3> > const &, int, int, RunOn) {}
};

} // namespace amrex

// ---- launch macros: plain serial loops --------------------------------------

#define AMREX_LOOP_3D(bx, i, j, k, block)                                      \
  for (int k = (bx).lo[2]; k <= (bx).hi[2]; ++k)                               \
    for (int j = (bx).lo[1]; j <= (bx).hi[1]; ++j)                             \
      for (int i = (bx).lo[0]; i <= (bx).hi[0]; ++i) block

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(fl, a1, a2, a3, b1, b2, b3,   \
                                                 c1, c2, c3)                   \
  do {                                                                         \
    (void)(fl);                                                                \
    {                                                                          \
      const ::amrex::Box a2 = (a1);                                            \
      a3                                                                       \
    }                                                                          \
    {                                                                          \
      const ::amrex::Box b2 = (b1);                                            \
      b3                                                                       \
    }                                                                          \
    {                                                                          \
      const ::amrex::Box c2 = (c1);                                            \
      c3                                                                       \
    }                                                                          \
  } while (0)

#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(fl, bx, nc, i, j, k, n, block)  \
  do {                                                                         \
    (void)(fl);                                                                \
    const ::amrex::Box _b = (bx);                                              \
    for (int n = 0; n < (nc); ++n)                                             \
      for (int k = _b.lo[2]; k <= _b.hi[2]; ++k)                               \
        for (int j = _b.lo[1]; j <= _b.hi[1]; ++j)                             \
          for (int i = _b.lo[0]; i <= _b.hi[0]; ++i) block                     \
  } while (0)

// ---- the slice of Cactus the operator uses ----------------------------------

using CCTK_REAL = double;
#define CCTK_VINFO(...)                                                        \
  do {                                                                         \
    std::fprintf(stderr, __VA_ARGS__);                                         \
    std::fprintf(stderr, "\n");                                                \
  } while (0)
#define CCTK_VERROR(...)                                                       \
  do {                                                                         \
    std::fprintf(stderr, __VA_ARGS__);                                         \
    std::fprintf(stderr, "\n");                                                \
    std::abort();                                                              \
  } while (0)
#define CCTK_ERROR(msg)                                                        \
  do {                                                                         \
    std::fprintf(stderr, "CCTK_ERROR: %s\n", msg);                             \
    std::abort();                                                              \
  } while (0)

namespace CarpetX {
constexpr int dim = 3;
}

#endif
