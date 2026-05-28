#ifndef CARPETX_CARPETX_IO_SILO_PLANES_HXX
#define CARPETX_CARPETX_IO_SILO_PLANES_HXX

#include "io_planes.hxx"

#include <cctk.h>

#ifdef HAVE_CAPABILITY_Silo

#include <string>
#include <vector>

namespace CarpetX {

void OutputSiloPlanes(const cGH *cctkGH, const std::vector<bool> &output_group,
                      const std::vector<plane_spec_t> &planes,
                      const std::string &output_dir,
                      const std::string &output_file);

} // namespace CarpetX

#endif

#endif // #ifndef CARPETX_CARPETX_IO_SILO_PLANES_HXX
