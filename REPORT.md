# Implementation report: flux-preserving prolongation for the staggered vector potential

Branch `vecpot-prolongation`, six commits on top of `18d81f96`.

This implements the operator derived in `PROLONG_A.md` — the edge-staggered
reformulation of Tóth & Roe (2002) for a vector potential $A_i$ on cell edges.
Read `PROLONG_A.md` first for the construction; this document covers what was
built, how far it is verified, what is deliberately missing, and what to expect
when it is first compiled.

**The code has never been compiled against real AMReX or Cactus.** There is no
Cactus tree on the machine it was written on, so `agent_scripts/build.sh` could
not run. §3 explains what was done instead and what that does and does not
cover; §7 says what to expect on first build.

---

## 1. What it does

A group declares that it is one component of an edge-staggered vector potential
by carrying, on **every** member:

```
TAGS='prolongation_type="vecpot" vector_potential="Avec_x Avec_y Avec_z"'
```

with the centerings `{cvv}`, `{vcv}`, `{vvc}` in x, y, z order. CarpetX then
prolongs the three components in one coupled pass instead of group by group.

The operator

* reproduces the coarse magnetic fluxes exactly, so it is compatible with the
  restriction `amrex::average_down_edges` already performs;
* gives $\nabla\cdot B = 0$ identically on the fine grid;
* produces **no magnetic field at all** from a pure gauge coarse field — the
  property that distinguishes it from interpolating $A_i$ with a generic
  stencil, which is divergence-free but not flux-conserving;
* never overwrites a fine edge that the fine grid evolves.

---

## 2. Commits

| # | Commit | Step of `PROLONG_A.md` |
|---|---|---|
| 1 | `5742d8e6` | 0 — kernels, interpolater, standalone tests |
| 2 | `e8ab5700` | 1 — group-triple discovery |
| 3 | `06f98dd9` | 2 — `MakeNewLevelFromCoarse` |
| 4 | `820e8fee` | 3 — refined patch, pre-fill, solve masks; the ghost paths |
| 5 | `467e6b9c` | 4 — `RemakeLevel` |
| 6 | `208adaff` | 6 — `TestVecPotProlongation`, status |

4447 insertions, 13 deletions across 32 files.

### New files

```
CarpetX/src/prolongate_3d_rf2_vecpot_kernels.hxx   the arithmetic, dependency-free
CarpetX/src/prolongate_3d_rf2_vecpot.hxx/.cxx      the amrex::Interpolater subclass
CarpetX/src/prolongate_3d_rf2_vecpot_groups.cxx    the vector_potential tag
CarpetX/src/prolongate_3d_rf2_vecpot_test.cxx      the startup self-test
CarpetX/test/unit/                                 standalone tests + AMReX shim
TestVecPotProlongation/                            integration test thorn
PROLONG_A.md                                       design, derivation, verification
```

### Changed files

`driver.{hxx,cxx}` (GroupData fields, operator selection, two dispatch sites),
`fillpatch.{hxx,cxx}` (three new entry points), `sync_restrict.cxx` (three
dispatch sites, partial-triple SYNC), `param.ccl`, `make.code.defn`.

---

## 3. How it was verified without a toolkit build

The arithmetic lives in `prolongate_3d_rf2_vecpot_kernels.hxx`, free of AMReX
and Cactus dependencies and templated on the array accessor. A small
*functional* AMReX/Cactus shim (`CarpetX/test/unit/shim/`, ~290 lines) then lets
the **real interpolater source** be compiled and executed with nothing but a
C++17 compiler:

```
./CarpetX/test/unit/run.sh        # ~2 seconds, no Cactus, no AMReX build
```

Four binaries. The operator's own headers carry one guard,
`CARPETX_VECPOT_STANDALONE_TEST`, whose only effect is to skip
`#include "driver.hxx"`; CarpetX itself never defines it.

| binary | covers |
|---|---|
| `test_vecpot_kernels` | the arithmetic, against the identities the construction must satisfy |
| `test_vecpot_interp` | **the real `prolongate_3d_rf2_vecpot.cxx`**: iteration boxes, index-type conversions, `CoarseBox`, the three launches, mask handling |
| `test_vecpot_selftest` | the self-test CarpetX runs at startup, so it is known to build and pass before it reaches the toolkit |
| `test_vecpot_tags` | the `vector_potential` validator against a synthetic group registry — chiefly its error paths |

### Results

All at a deliberately non-unit cell aspect ratio $(\Delta X,\Delta Y,\Delta Z) =
(0.7,\,1.3,\,0.35)$, which is where the construction is least obviously correct:

```
[1] owning-type boxes partition the fine edges
  distinct fine edges written                          13872  ok
  edges written more than once                     0.000e+00  ok

[2] restriction compatibility (Eq. 7) and (Eq. 11), (Eq. 15)
  coarse flux equals the mean of its four fine fluxes  3.553e-15  ok
  fine surface fluxes match Toth Eq. (11)              3.553e-15  ok
  fine interior fluxes match Toth Eq. (15)             3.109e-15  ok

[3] divergence and gauge
  fine div B (random coarse data)                      8.882e-15  ok
  fine B from a pure gauge coarse field                3.997e-15  ok
  prolonged A equals grad of trilinear lambda (Eq. 31) 1.776e-15  ok

[4] copy-from-fine rule
  masked fine edges are not overwritten                0.000e+00  ok
  copied halves differ from the equal split by         9.335e-01  (nonzero: the
                                                                  test is real)
  mixed-provenance coarse faces still conserve flux    4.441e-16  ok
  fine div B with mixed provenance                     7.105e-15  ok

[5] accuracy
  observed convergence order                               2.589  ok
```

The discrete curls used for *verification* are written out longhand in the
tests, independently of the cyclic templates in the header, so a mistake in the
cyclic direction algebra cannot cancel against itself.

### What the shim does and does not cover

It catches everything that is our own — index arithmetic, template parameters,
box conversions, macro-expansion hazards. It **cannot** catch a signature
mismatch against the real AMReX. Its size is also a statement of how little
AMReX surface the operator touches.

It earned its keep twice:

1. **A real compile error.** `surface_interp<D, P1, P2>(...)` cannot appear
   inside `AMREX_LOOP_3D`, which is a fixed five-argument macro and splits its
   body on top-level commas; angle brackets do not group. This would have
   failed the toolkit build. Hence the packed wrapper
   `surface_interp_n<D*4 + P1*2 + P2>`.

2. **Two invariants that would otherwise be silent corruption.** A test
   allocates the coarse buffers at exactly `CoarseBox(target)` and the fine
   buffers at exactly the target, for four differently offset boxes including
   negative indices. The shim's array accessor aborts on out-of-range access,
   so a clean run pins both that `CoarseBox` is not under-sized and that the
   operator writes nothing outside its target. In a real build either would
   corrupt memory quietly rather than crash.

One further bug was found by review rather than by test: the tag validator's
symmetry check originally recursed into the full validator, so two members of a
triple would have recursed until the stack overflowed at startup. It is now
split into a raw parse plus validation.

---

## 4. Design notes worth knowing before reviewing

**Three launches, indexed by the coarse object that owns each fine edge.** Each
of the four cases iterates over `amrex::convert(c_fine_region, T)` with `T` the
index type of the owning coarse object — the coarse edge, one of the two coarse
face types, or the coarse cell. This is the draft's own decomposition (24 + 24 +
6 fine edges per coarse cell) and it is an exact partition, verified. The
owning-type boxes are self-clipping; a single nodal box for all four cases would
write out of bounds at the top nodal index.

**The interior closure reads only the destination array.** That is what makes
the copy-from-fine rule work without a branch: the closure sees whatever values
are on the coarse cell's faces, prolonged or copied. It also matches
`amrex::facediv_int`'s signature, which likewise takes no coarse array.

**No cell-aspect-ratio weights anywhere.** `amrex::facediv_int` — the same
operator for face-staggered $B$ — carries twelve five-term expressions with
weights like `dx*(2*dz*dz+dx*dx)/(8*dy*(dz*dz+dx*dx))`. That structure is
exactly what the edge formulation removes. **If an aspect-ratio weight ever
appears in `prolongate_3d_rf2_vecpot_kernels.hxx`, it is a bug.** Useful
invariant to review against.

### Two deliberate deviations from `amrex::FaceDivFree`

1. **Driven through `interp_vecpot`, not `interp_arr`.** Edge data has four
   classes of "known" object where face data has one, so seven solve masks are
   needed; AMReX's `interp_arr` signature carries three. Nothing inside AMReX
   can call us anyway — its Array FillPatch path asserts
   `ixType().nodeCentered(d)`, which edge data fails — so CarpetX owns the call
   site. Both virtual overrides fail loudly.

2. **Takes the coarse cell size, not a `Geometry`.** The cell size is all the
   operator uses; taking it directly keeps it a pure function of its inputs and
   testable without constructing a `Geometry`.

### Other choices

* The ghost path runs **synchronously** rather than through the task managers.
  The scalar path overlaps its coarse-patch copy with other groups' work; doing
  the same here means interleaving seven mask builds and three patch copies,
  which is worth doing only once the operator has been exercised in a real run.
* A `SYNC` naming only part of a triple is **quietly extended** to all three,
  with a one-time warning, rather than rejected. The components are coupled and
  the others would have had to be synced anyway.
* A triple is treated as fillable only when **every** member is: a coupled
  operator cannot prolong two thirds of a vector potential.

---

## 5. Not implemented

### Subcycling (Step 5)

Deliberate. It needs the time-blended coarse patch that
`FillPatch_Prolongate`'s `w_new`/`cmfab_old` provide, and **a blended coarse
edge is not the restriction of either fine state**. Flux conservation across the
refinement boundary would therefore degrade by the time-interpolation error
rather than hold exactly. That is a physics decision, not a plumbing one, and it
should be made against a measured number.

Until it is made, a vector-potential group prolonged from a time-misaligned
coarse level raises an error naming the limitation, instead of silently
producing a slightly non-conservative field. If you want subcycling, the useful
first step is to measure the seam residual (§7.8 of `PROLONG_A.md`) under a
blend and decide whether its size is acceptable.

### The high-order Stage-1 split (Eq. 32)

Also deliberate. It is pure gauge — it changes $A$ but not $B$ — but the draft
requires it to be added *together with* compensating terms in the face and cell
formulas, and those have not been derived or verified here. Implementing it
blind seemed worse than leaving it out.

### The gauge scalar $\Phi$

Guidance only, since there is no $\Phi$ group in this repository. Tag it
`prolongation_type="ddf"`, `prolongation_order=1`: trilinear on a `{vvv}` group
is already `prolongate_ddf_3d_rf2_c000_o1`, and §2.7 of `PROLONG_A.md` confirms
that our $P_\lambda$ is exactly trilinear, so the prolonged $(A_i, \Phi)$ pair
is then the exact gauge transform of the coarse pair.

---

## 6. Open questions to measure, not assume

Both are in §8 of `PROLONG_A.md` with more detail.

**Reflection symmetry.** The correction mixes components, and CarpetX applies
parities per component in `apply_boundary_conditions` with no knowledge of the
coupling. If the coarse patch's boundary fill is parity-consistent the discrete
curl comes out right automatically, but this has not been checked. A sign error
here is invisible in periodic tests, which is what the current parameter file
runs. Worth adding a reflection variant to `TestVecPotProlongation`.

**Seam conservation.** A coarse face straddling the edge of `ba_fine_patch`
receives a mix of copied and prolonged edges. This is benign — verified — but
only because the coarse data at and under the fine boundary is the restriction
of the fine data. That is ordinary AMR hygiene, but here it is load-bearing for
flux conservation rather than just accuracy, so it is worth asserting rather
than assuming.

---

## 7. Bring-up

Suggested order:

1. `./CarpetX/test/unit/run.sh` — should pass in seconds, and confirms the
   arithmetic and the operator's own logic before anything else is in play.
2. `./agent_scripts/build.sh` — **expect this to need a round or two.** The
   shim catches our errors, not API drift. The likeliest friction points, all
   read off `../amrex` but not compiled:
   * `make_mf_crse_mask` / `make_mf_crse_patch` visibility in `amrex::detail`
     for your AMReX version (`fillpatch.cxx` already relies on the same
     `using namespace amrex::detail` guarded by `AMREX_RELEASE_NUMBER >= 240500`);
   * the `FabArrayBase::CPC` constructor arguments;
   * `ParallelCopyToGhost`'s exact signature.
3. `./agent_scripts/test.sh` — the existing suite should be unaffected; nothing
   selects the new operator unless a group asks for it.
4. `TestVecPotProlongation/test/puregauge.par`, watching `max |Bvec_*|`.

**The pass criterion that matters is the fourth one.** `max |Bvec_x|`,
`|Bvec_y|`, `|Bvec_z|` must stay at round-off on both levels. Anything larger
means the refinement boundary is generating magnetic field from a pure gauge
field, which is the single failure the whole construction exists to prevent.
`divB` is identically zero for any edge-staggered $A$, prolonged or not, so it
checks the test thorn's bookkeeping rather than the operator.

The initial data is the *exact* discrete gradient of a smooth scalar —
differenced between precisely the two vertices bounding each edge — so $B$
vanishes identically rather than to truncation order. A failure is therefore
unambiguous.
