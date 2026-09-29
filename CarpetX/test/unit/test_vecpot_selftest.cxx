// Compiles and runs CarpetX/src/prolongate_3d_rf2_vecpot_test.cxx -- the
// self-test that CarpetX runs at startup -- against the AMReX/Cactus shim, so
// that it is known to build and pass before it reaches the toolkit.
#define CARPETX_VECPOT_STANDALONE_TEST 1
#include "shim/amrex_shim.hxx"

#include "../../src/prolongate_3d_rf2_vecpot.cxx"
#include "../../src/prolongate_3d_rf2_vecpot_test.cxx"

int main() {
  CarpetX::test_prolongate_3d_rf2_vecpot();
  std::printf("PASSED (self-test did not abort)\n");
  return 0;
}
