#include "planes.hxx"

#include <cctk.h>

#include <AMReX_Box.H>
#include <AMReX_IntVect.H>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace PlanesX {

namespace {

long long pow10ll(int n) {
  long long r = 1;
  for (int i = 0; i < n; ++i)
    r *= 10;
  return r;
}

std::string trim(std::string s) {
  const auto ws = [](unsigned char c) { return std::isspace(c); };
  while (!s.empty() && ws(s.front()))
    s.erase(s.begin());
  while (!s.empty() && ws(s.back()))
    s.pop_back();
  return s;
}

} // namespace

std::string format_plane_tag(int normal_axis, CCTK_REAL rounded_elevation,
                             int int_precision, int frac_precision) {
  assert(normal_axis >= 0 && normal_axis < 3);
  assert(int_precision >= 1 && frac_precision >= 0);

  const char *plane_name;
  const char *axis_name;
  switch (normal_axis) {
  case 0:
    plane_name = "yz";
    axis_name = "x";
    break;
  case 1:
    plane_name = "xz";
    axis_name = "y";
    break;
  case 2:
    plane_name = "xy";
    axis_name = "z";
    break;
  default:
    assert(0);
    plane_name = "";
    axis_name = "";
  }

  const char *sign = rounded_elevation < 0 ? "neg" : "pos";
  const double mag = std::abs(double(rounded_elevation));
  const long long scale = pow10ll(frac_precision);
  const long long total = std::llround(mag * double(scale));
  const long long int_part = total / scale;
  const long long frac_part = total % scale;

  std::ostringstream buf;
  buf << plane_name << "_" << axis_name << "_" << sign
      << std::setw(int_precision) << std::setfill('0') << int_part;
  if (frac_precision > 0)
    buf << "p" << std::setw(frac_precision) << std::setfill('0') << frac_part;
  return buf.str();
}

std::vector<plane_spec_t> parse_planes(const std::string &spec,
                                       int int_precision, int frac_precision) {
  static std::set<std::string> warned_specs;

  std::vector<plane_spec_t> out;
  if (spec.empty())
    return out;

  assert(int_precision >= 1 && frac_precision >= 0);
  const long long scale = pow10ll(frac_precision);
  const long long max_int = pow10ll(int_precision) - 1;

  std::size_t pos = 0;
  while (pos < spec.size()) {
    const std::size_t comma = spec.find(',', pos);
    std::string item = trim(spec.substr(
        pos, comma == std::string::npos ? std::string::npos : comma - pos));
    pos = (comma == std::string::npos) ? spec.size() : comma + 1;
    if (item.empty())
      continue;

    const std::size_t colon = item.find(':');
    if (colon == std::string::npos) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "Plane spec '%s' missing ':' separator; skipping",
                 item.c_str());
      continue;
    }

    const std::string axes = trim(item.substr(0, colon));
    const std::string elev_str = trim(item.substr(colon + 1));

    int normal_axis;
    if (axes == "yz")
      normal_axis = 0;
    else if (axes == "xz")
      normal_axis = 1;
    else if (axes == "xy")
      normal_axis = 2;
    else {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "Plane spec '%s': unknown axes '%s' (expected xy/xz/yz); "
                 "skipping",
                 item.c_str(), axes.c_str());
      continue;
    }

    CCTK_REAL elev;
    try {
      std::size_t end = 0;
      elev = CCTK_REAL(std::stod(elev_str, &end));
      if (end != elev_str.size())
        throw std::invalid_argument("trailing characters");
    } catch (const std::exception &) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "Plane spec '%s': unparseable elevation '%s'; skipping",
                 item.c_str(), elev_str.c_str());
      continue;
    }

    const double sign_factor = elev < 0 ? -1.0 : 1.0;
    const double mag = std::abs(double(elev));
    const long long total = std::llround(mag * double(scale));
    const CCTK_REAL rounded =
        CCTK_REAL(sign_factor * (double(total) / double(scale)));

    const bool first_time = warned_specs.insert(item).second;
    if (first_time) {
      const double tol = 1e-12 * std::max(1.0, mag);
      if (std::abs(double(rounded - elev)) > tol)
        CCTK_VWARN(CCTK_WARN_ALERT,
                   "Plane elevation %.17g exceeds "
                   "planes_frac_precision=%d digits; degraded to %.17g.",
                   double(elev), frac_precision, double(rounded));
      if (total / scale > max_int)
        CCTK_VWARN(CCTK_WARN_ALERT,
                   "Plane elevation %g exceeds planes_int_precision=%d "
                   "digits; tag field widens (lexicographic sort across "
                   "planes may not match numeric sort).",
                   double(rounded), int_precision);
    }

    plane_spec_t p;
    p.normal_axis = normal_axis;
    p.elevation = rounded;
    p.tag =
        format_plane_tag(normal_axis, rounded, int_precision, frac_precision);
    out.push_back(std::move(p));
  }
  return out;
}

int snap_to_grid_index(const plane_spec_t &plane, const amrex::Geometry &geom,
                       const amrex::IndexType &indextype) {
  const int d = plane.normal_axis;
  const double x0 = geom.ProbLo()[d];
  const double dx = geom.CellSize()[d];
  const amrex::Box &domain = geom.Domain();
  const int ilo = domain.smallEnd(d);
  const int ihi = domain.bigEnd(d);

  const double r = (double(plane.elevation) - x0) / dx;
  int i_snap;
  int i_lo, i_hi;
  if (indextype.cellCentered(d)) {
    i_snap = ilo + int(std::lround(r - 0.5));
    i_lo = ilo;
    i_hi = ihi;
  } else {
    i_snap = ilo + int(std::lround(r));
    i_lo = ilo;
    i_hi = ihi + 1;
  }
  if (i_snap < i_lo || i_snap > i_hi)
    return -1;
  return i_snap;
}

amrex::Box extract_slab(const amrex::FArrayBox &fab, int normal_axis,
                        int slice_idx, int numvars,
                        std::vector<CCTK_REAL> &out_buf) {
  assert(normal_axis >= 0 && normal_axis < 3);
  assert(numvars >= 0 && numvars <= fab.nComp());

  const amrex::Box &box = fab.box();
  if (slice_idx < box.smallEnd(normal_axis) ||
      slice_idx > box.bigEnd(normal_axis))
    return amrex::Box();

  const int nx = box.length(0);
  const int ny = box.length(1);
  const int nz = box.length(2);
  const std::ptrdiff_t cell_count = std::ptrdiff_t(nx) * ny * nz;

  int a, b;
  switch (normal_axis) {
  case 0:
    a = 1;
    b = 2;
    break;
  case 1:
    a = 0;
    b = 2;
    break;
  case 2:
    a = 0;
    b = 1;
    break;
  default:
    assert(0);
    a = 0;
    b = 1;
  }
  const int na = box.length(a);
  const int nb = box.length(b);

  out_buf.resize(std::size_t(numvars) * std::size_t(na) * std::size_t(nb));

  const std::ptrdiff_t strides[3] = {1, std::ptrdiff_t(nx),
                                     std::ptrdiff_t(nx) * ny};
  const std::ptrdiff_t s_a = strides[a];
  const std::ptrdiff_t s_b = strides[b];
  const std::ptrdiff_t s_n = strides[normal_axis];
  const std::ptrdiff_t base = (slice_idx - box.smallEnd(normal_axis)) * s_n;

  const CCTK_REAL *const src = fab.dataPtr();
  CCTK_REAL *const dst = out_buf.data();

  for (int vi = 0; vi < numvars; ++vi) {
    const CCTK_REAL *const src_v = src + vi * cell_count + base;
    CCTK_REAL *const dst_v =
        dst + std::size_t(vi) * std::size_t(na) * std::size_t(nb);
    for (int j = 0; j < nb; ++j)
      for (int i = 0; i < na; ++i)
        dst_v[i + std::size_t(j) * na] = src_v[i * s_a + j * s_b];
  }

  amrex::IntVect lo = box.smallEnd();
  amrex::IntVect hi = box.bigEnd();
  lo.setVal(normal_axis, slice_idx);
  hi.setVal(normal_axis, slice_idx);
  return amrex::Box(lo, hi, box.ixType());
}

} // namespace PlanesX
