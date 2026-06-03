#ifndef CARPETX_PLANESX_PLANES_HXX
#define CARPETX_PLANESX_PLANES_HXX

#include <cctk.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>

#include <string>
#include <vector>

namespace PlanesX {

// normal_axis: 0 -> yz, 1 -> xz, 2 -> xy. elevation is already rounded
// to frac_precision; tag is the filename/meshname fragment.
struct plane_spec_t {
  int normal_axis;
  CCTK_REAL elevation;
  std::string tag;
};

// Rounds each elevation to frac_precision digits; warns once per unique
// spec on precision overflow.
std::vector<plane_spec_t> parse_planes(const std::string &spec,
                                       int int_precision, int frac_precision);

std::string format_plane_tag(int normal_axis, CCTK_REAL rounded_elevation,
                             int int_precision, int frac_precision);

// Returns -1 if elevation lies outside the level along normal_axis.
// VC indextype snaps to x0+i*dx; CC to x0+(i+0.5)*dx.
int snap_to_grid_index(const plane_spec_t &plane, const amrex::Geometry &geom,
                       const amrex::IndexType &indextype);

// Collapses fab along normal_axis at slice_idx into out_buf (size
// numvars * in-plane cells, var-major). Returns empty box if out of range.
amrex::Box extract_slab(const amrex::FArrayBox &fab, int normal_axis,
                        int slice_idx, int numvars,
                        std::vector<CCTK_REAL> &out_buf);

} // namespace PlanesX

#endif // #ifndef CARPETX_PLANESX_PLANES_HXX
