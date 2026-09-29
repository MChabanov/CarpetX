# Standalone unit tests for the vector-potential prolongation

    ./run.sh

Nothing is needed but a C++17 compiler: no Cactus, no AMReX build, no MPI.
The point is to make the operator verifiable on a laptop, since the Einstein
Toolkit testsuite can only be run where a full build exists.

| binary | what it covers |
|---|---|
| `test_vecpot_kernels` | the arithmetic in `src/prolongate_3d_rf2_vecpot_kernels.hxx`, against the identities the construction must satisfy |
| `test_vecpot_interp` | the real interpolater, `src/prolongate_3d_rf2_vecpot.cxx`, compiled against the shim: iteration boxes, index-type conversions, `CoarseBox`, the three launches, mask handling |
| `test_vecpot_selftest` | `src/prolongate_3d_rf2_vecpot_test.cxx`, the self-test CarpetX runs at startup, so that it is known to build and pass before it reaches the toolkit |

## What is checked

All at a deliberately non-unit cell aspect ratio (dX, dY, dZ = 0.7, 1.3, 0.35),
which is where the construction is least obviously correct:

* the four owning-type iteration boxes partition the fine edges exactly
  (36 + 24 + 24 + 16 per component per 2^3 coarse cells, disjoint and covering);
* restriction compatibility -- a coarse flux equals the mean of its four fine
  fluxes, which is what makes the operator compatible with
  `amrex::average_down_edges`;
* the fine fluxes agree with Toth & Roe's face-based formulas, both on the
  surface of a coarse cell and in its interior;
* a pure gauge coarse field, A = grad lambda, produces exactly zero fine
  magnetic field -- the property that distinguishes this operator from generic
  Lagrange interpolation of A_i;
* gauge covariance: the prolonged A is the discrete gradient of the *trilinear*
  interpolant of lambda;
* the fine field is exactly divergence free;
* masked (already refined) edges are never overwritten, and a coarse face that
  receives a mix of copied and prolonged edges still conserves flux;
* second-order convergence of B on a smooth field.

## The shim

`shim/` is a small functional stand-in for the slice of AMReX and Cactus that
the operator touches.  Its size is a statement of how much AMReX surface the
operator depends on.  It is test-only and is never compiled into CarpetX.

The operator's headers carry one guard, `CARPETX_VECPOT_STANDALONE_TEST`, whose
only effect is to skip `#include "driver.hxx"`; CarpetX itself never defines it.

The shim is not a substitute for building against real AMReX -- it cannot catch
a signature mismatch -- but it does catch everything that is our own: index
arithmetic, template parameters, box conversions, and macro-expansion hazards.
It has already earned its keep once, by catching that `surface_interp<D,P1,P2>`
cannot appear inside `AMREX_LOOP_3D`, whose top-level commas split the template
argument list.  That is why the kernels expose `surface_interp_n<D*4+P1*2+P2>`.
