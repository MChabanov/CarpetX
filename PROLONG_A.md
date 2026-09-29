# Flux-preserving prolongation for the staggered vector potential

Implementation plan for the `interp_arr` route.

Reference: *"Flux-preserving prolongation for vector potential evolution with
mesh refinement"* (draft, `apssamp.tex`), Sec. V — the edge-staggered
reformulation of Tóth & Roe (2002).

Terminology: the draft's construction has three **stages** (coarse edges,
coarse faces, cell interior). The implementation plan in §6 is divided into
**steps**. The two are unrelated; see §6.0.

---

## 0. Decisions taken

| question | decision |
|---|---|
| which stages of the construction | **all three.** Stages 1–3 are needed in every prolongation path anyway (§2.4); nothing is scoped out. |
| copy from an already refined neighbour (draft Eq. 33) | **in scope, and it is the governing rule.** Where fine data exists, it wins: such an edge is not part of the prolongation at all. This makes AMReX's `solve_mask` load-bearing (§4, §5.2) and forces the moment form of Stage 3 (§2.5). |
| limiters | **unlimited centred slopes** for now. Limiting stays deferred (§9). |
| standalone thorn | **deferred, not rejected.** The kernel and its `Interpolater` subclass could live outside CarpetX — thorns here already `REQUIRES AMReX` (`CarpetXRegrid/configuration.ccl`) and `driver.hxx` is exported (`CarpetX/interface.ccl:17`) — but the `interp_arr` path it needs cannot; that is Steps 2–5. Extract after Step 3, when the registry interface can be shaped around a real consumer. See §5.7. |

The one place the "all stages" decision is left narrow: the **high-order
Stage-1 split** (draft Eq. 32) is pure gauge — it changes $A$ but not $B$ — so
it ships as a parameter defaulting to off (Step 5). Say so if you want it on by
default.

---

## 1. Goal and scope

Add a prolongation operator for a vector potential $A_i$ staggered on cell
edges that

* reproduces the coarse magnetic fluxes exactly (conservative across a
  refinement boundary),
* yields $\nabla\cdot B = 0$ identically on the fine grid (automatic for any
  edge representation),
* generates **no** magnetic field from a pure gauge coarse field,
* is local: a fine edge on a coarse edge uses only that coarse edge; a fine
  edge in a coarse face uses only that face,
* never overwrites a fine edge that the fine level already evolves.

**Out of scope:** limited slopes, and genuinely higher-order *field*
prolongation. See §9.

---

## 2. The operator in CarpetX variables

### 2.1 Storage convention

CarpetX stores **components** $\tilde A_i$, not circulations. With
$\mathcal A_x = \Delta X\,\tilde A_x$, $W = \Delta X\Delta Y\,B^z$,
$V = \Delta Z\Delta X\,B^y$, and the slope definition
$W_y = \tfrac18\big(W^{0,4,k} - W^{0,-4,k}\big)$, Eq. (30) of the draft
collapses to a single rule:

$$\tilde a_x^{\rm fine} \;=\; \big\langle \tilde A_x\big\rangle
\;+\; \frac{\Delta Y^2}{8}\big\langle \partial_y B^z\big\rangle
\;-\; \frac{\Delta Z^2}{8}\big\langle \partial_z B^y\big\rangle$$

where $\langle\cdot\rangle$ is the arithmetic mean over the coarse edges (or
faces) the fine edge sits between, and **each correction term appears only if
averaging actually happened in that direction**. Cyclic for $A_y$, $A_z$.

Mnemonic: *halve, average, and add one quarter of the relevant transverse flux
slope once for each direction in which averaging was performed.*

### 2.2 Centering and index types

| quantity | CENTERING | CarpetX `indextype` | AMReX `IndexType` |
|---|---|---|---|
| $A_x$ (x-edge) | `{cvv}` | `(1,0,0)` | CELL, NODE, NODE |
| $A_y$ (y-edge) | `{vcv}` | `(0,1,0)` | NODE, CELL, NODE |
| $A_z$ (z-edge) | `{vvc}` | `(0,0,1)` | NODE, NODE, CELL |
| $B^x$ (x-face) | `{vcc}` | `(0,1,1)` | NODE, CELL, CELL |
| $B^y$ (y-face) | `{cvc}` | `(1,0,1)` | CELL, NODE, CELL |
| $B^z$ (z-face) | `{ccv}` | `(1,1,0)` | CELL, CELL, NODE |
| $\Phi$ (vertex) | `{vvv}` | `(0,0,0)` | NODE, NODE, NODE |

> **Gotcha.** CarpetX and AMReX use *opposite* constants: CarpetX
> `indextype[d] == 0` means vertex/node, AMReX `IndexType::NODE == 1`. The
> conversion lives at `CarpetX/src/driver.cxx:991` with a comment. Every new
> index-type manipulation must go through the same convention.

`StaggeredWaveToyX/interface.ccl` is the working precedent for declaring such
groups, including the `parities={-1 +1 +1}` tags each component needs.

### 2.3 Discrete curls

$B$ is never stored; it is reconstructed on the fly by discrete Stokes. With
$(i,j,k)$ the index of the face in its own index space and $(h_X,h_Y,h_Z)$ the
spacings of the level in question:

```
Bx(i,j,k) = [ Ay(i,j,k) - Ay(i,j,k+1) ] / hZ
          + [ Az(i,j+1,k) - Az(i,j,k) ] / hY          // x-face {vcc}

By(i,j,k) = [ Az(i,j,k) - Az(i+1,j,k) ] / hX
          + [ Ax(i,j,k+1) - Ax(i,j,k) ] / hZ          // y-face {cvc}

Bz(i,j,k) = [ Ax(i,j,k) - Ax(i,j+1,k) ] / hY
          + [ Ay(i+1,j,k) - Ay(i,j,k) ] / hX          // z-face {ccv}
```

These are exactly Eqs. (5a)–(5c) divided by the face areas. They are the only
place where the three components meet. They are used on the coarse grid
(§2.4) and, with $h = \Delta/2$, on the fine grid (§2.5).

### 2.4 Stages 1 and 2, and the naive Stage 3

In the prolongation kernel CarpetX computes
(`prolongate_3d_rf2_impl.hxx:1423`)

```cpp
const vect<int,dim> icrse = ifine >> 1;
const vect<int,dim> off   = ifine & 0x1;
```

For $A_x$ with centering `(CC, VC, VC)`, write `off = (ox, oy, oz)` and
`icrse = (ix, iy, iz)`. **`ox` never enters** — the equal-split gauge of
Stage 1 gives both halves of a coarse edge the same value. The four cases of
Table I are just the four parities of `(oy, oz)`:

```
// (oy,oz) = (0,0)  --  fine edge lies on a coarse edge (Stage 1)
a = Ax(ix, iy, iz);

// (oy,oz) = (1,0)  --  mid-line in y of the coarse z-face (ix,iy,iz)  (Stage 2)
a = 0.5 * ( Ax(ix,iy,iz) + Ax(ix,iy+1,iz) )
  + (dY/16) * ( Bz(ix,iy+1,iz) - Bz(ix,iy-1,iz) );

// (oy,oz) = (0,1)  --  mid-line in z of the coarse y-face (ix,iy,iz)  (Stage 2)
a = 0.5 * ( Ax(ix,iy,iz) + Ax(ix,iy,iz+1) )
  - (dZ/16) * ( By(ix,iy,iz+1) - By(ix,iy,iz-1) );

// (oy,oz) = (1,1)  --  interior of the coarse cell (ix,iy,iz)  (Stage 3, naive)
a = 0.25 * sum_{dy,dz in {0,1}} Ax(ix, iy+dy, iz+dz)
  + (dY/32) * sum_{dz in {0,1}} ( Bz(ix,iy+1,iz+dz) - Bz(ix,iy-1,iz+dz) )
  - (dZ/32) * sum_{dy in {0,1}} ( By(ix,iy+dy,iz+1) - By(ix,iy+dy,iz-1) );
```

For $A_y$ the relevant parities are `(oz, ox)`, the "+" term is $B^x$ averaged
in $z$ and the "−" term is $B^z$ averaged in $x$; for $A_z$ the parities are
`(ox, oy)`, "+" is $B^y$ averaged in $x$, "−" is $B^x$ averaged in $y$:

```
A_y, (oz,ox) = (1,0):  0.5*( Ay(ix,iy,iz) + Ay(ix,iy,iz+1) )
                       + (dZ/16)*( Bx(ix,iy,iz+1) - Bx(ix,iy,iz-1) );
A_y, (oz,ox) = (0,1):  0.5*( Ay(ix,iy,iz) + Ay(ix+1,iy,iz) )
                       - (dX/16)*( Bz(ix+1,iy,iz) - Bz(ix-1,iy,iz) );
A_z, (ox,oy) = (1,0):  0.5*( Az(ix,iy,iz) + Az(ix+1,iy,iz) )
                       + (dX/16)*( By(ix+1,iy,iz) - By(ix-1,iy,iz) );
A_z, (ox,oy) = (0,1):  0.5*( Az(ix,iy,iz) + Az(ix,iy+1,iz) )
                       - (dY/16)*( Bx(ix,iy+1,iz) - Bx(ix,iy-1,iz) );
```

Every index above type-checks: each $B$ lookup lands on the index type of the
face it names. That is itself a useful consistency check when transcribing.

> **All four parities occur in every prolongation path.** They are parities of
> the fine index, nothing more. Ghost fills at a refinement boundary hit the
> cell-interior case exactly as often as a regrid does. There is no path that
> needs "only Stages 1 and 2".

### 2.5 Stage 3, provenance-agnostic (the form to implement)

The Stage-3 expression above is written in terms of **coarse** slopes, which is
only valid when the coarse cell's six faces were themselves prolonged by
Eqs. (11)/(28). Under the decision of §0 some of those faces are instead
**copied from an already refined neighbour**, and the draft's Eq. (33) says the
slope must then be recovered as a **moment of the cell's own fine surface
fluxes**:

$$W_y^{0,0,2s} = \sum_{i,j=\pm1} j\,w^{i,j,2s}.$$

This is not an alternative formula for the copied case — it is the *same*
formula in both cases. Substituting Eq. (11) into the moment reproduces the
coarse slope identically (verified in §2.7). So **implement Stage 3 in moment
form only**, and the copy case needs no branch at all:

Going one step further, the $\langle A\rangle$ term can be read from the fine
array too: by (C2) the coarse edge value *is* the mean of its two fine halves,
so `<Ax> = ½(a_h0 + a_h1)` whether those halves were prolonged (equal, both
$=A$) or copied (mean $=A$ by restriction). Doing so makes Stage 3 depend on
**nothing but the destination array** — matching `amrex::facediv_int`'s
signature exactly (§4.1) and removing a precondition (§8.8). Verified in §2.7.

```
// A_x interior of coarse cell (ix,iy,iz).  Reads ONLY the fine array.
// ax_f/by_f/bz_f: fine edges, and fine fluxes on this cell's own y- and
// z-faces by discrete curl (2.3, h = d/2) of the edges already written by
// Stages 1-2 -- whether those were prolonged or copied from the fine level.
a = 0.125 * sum_{h,dy,dz in {0,1}} ax_f(2ix+h, 2(iy+dy), 2(iz+dz))
  + (dY/16) * sum_{s,a,b in {0,1}} (2b-1) * bz_f(2ix+a, 2iy+b,   2iz+2s)
  - (dZ/16) * sum_{s,a,c in {0,1}} (2c-1) * by_f(2ix+a, 2iy+2s, 2iz+c );

// A_y interior
a = 0.125 * sum_{h,dx,dz} ay_f(2(ix+dx), 2iy+h, 2(iz+dz))
  + (dZ/16) * sum_{s,b,c} (2c-1) * bx_f(2ix+2s, 2iy+b,   2iz+c )
  - (dX/16) * sum_{s,a,b} (2a-1) * bz_f(2ix+a,  2iy+b,   2iz+2s);

// A_z interior
a = 0.125 * sum_{h,dx,dy} az_f(2(ix+dx), 2(iy+dy), 2iz+h)
  + (dX/16) * sum_{s,a,c} (2a-1) * by_f(2ix+a,  2iy+2s, 2iz+c )
  - (dY/16) * sum_{s,b,c} (2b-1) * bx_f(2ix+2s, 2iy+b,  2iz+c );
```

The `s` index runs over the cell's two opposing faces in the relevant
direction. This is why the kernel must run in **two phases** (§5.2): Stage 3
reads the destination array written by Stages 1–2, including at the neighbouring
coarse index `(ix,iy,iz+1)`, so it cannot be fused into the same launch.

### 2.6 What the averaging part already is

The $\langle\tilde A_x\rangle$ part alone — injection along the edge, linear
transverse — is **exactly the existing operator**

```
prolongate_ddf_3d_rf2_c100_o1
  = prolongate_3d_rf2<CC,VC,VC, CONS,POLY,POLY, 0,1,1, FB_NONE>
```

Verified: `coeffs1d<CC,CONS,0>::coeffs == {1}` and `interp1d<CC,CONS,0>`
returns `crse(0)` for both `off` values (a copy into both halves);
`interp1d<VC,POLY,1>` returns `crse(0)` for `off==0` and
`½(crse(0)+crse(1))` for `off==1`. The tensor product of the two reproduces
all four cases of §2.4 with the correction terms dropped. This is the
regression baseline of §7.6.

Why the correction cannot be made single-component:

$$\partial_y B^z - \partial_z B^y
= \partial_x\big(\partial_y A_y + \partial_z A_z\big)
- \big(\partial_y^2 + \partial_z^2\big) A_x .$$

The first term needs $A_y$ and $A_z$. Dropping it and keeping only
$-(\partial_y^2+\partial_z^2)A_x$ would merely turn the average into a
higher-order interpolation of $A_x$ alone — the "generic Lagrange stencil on
$A_i$" that Sec. V C rejects: divergence-free but *not* flux-conserving.
**The coupling is essential.** This is what forces the `interp_arr` route.

### 2.7 Status: the formulas are verified

A standalone numerical check (numpy, random coarse edge data on a
$10^3$-coarse-cell block, spacings $\Delta X,\Delta Y,\Delta Z = 0.7, 1.3,
0.35$ — deliberately non-unit aspect ratio) confirms, for **all three
components and all four parity cases**, to machine precision
($\le 3\times10^{-15}$):

| check | result |
|---|---|
| restriction compatibility, Eq. (7): the four fine fluxes on each coarse face sum to the coarse flux | ✅ |
| Tóth surface fluxes, Eq. (11): $u^{\pm2,j,k} = \tfrac14(U + jU_y + kU_z)$ | ✅ |
| Tóth interior fluxes, Eq. (15) | ✅ |
| third-order moments $U_{xyz},V_{xyz},W_{xyz}$ vanish, as Sec. V A predicts when the surfaces come from Eq. (11) | ✅ |
| pure gauge $A = \nabla\lambda$ ⟹ fine $B \equiv 0$ | ✅ |
| gauge covariance, Eq. (31): prolonged $A$ equals the discrete gradient of the **trilinear** interpolant of $\lambda$ | ✅ |
| **moment form of Stage 3 (§2.5) ≡ coarse-slope form (§2.4)** | ✅ (4.4e-16) |
| **mixed coarse face** (some bounding edges copied from a fine field, the rest prolonged): the four fine fluxes still sum to the coarse flux | ✅ (≤5.6e-17), see §8.3 |
| **Stage 3 from the fine array only** — no coarse $A$ at all (§2.5, §4.1) | ✅ (4.4e-16) |

In component form the Eq. (11) reference reads, for a quadrant $(j,k)=(\pm1)$
of a coarse x-face,

$$b^x_{\rm fine} = B^x + \tfrac{j}{8}\big(B^x_{j+1}-B^x_{j-1}\big)
+ \tfrac{k}{8}\big(B^x_{k+1}-B^x_{k-1}\big)$$

and the Eq. (15) reference (third-order terms vanishing)

$$b^{x,\rm mid} = \tfrac12\big(b^{x,L} + b^{x,R}\big)
+ \frac{\Delta X}{8\Delta Y}\!\!\sum_{i,j,k=\pm1}\!\! ij\,b^{y,\,i,2j,k}
+ \frac{\Delta X}{8\Delta Z}\!\!\sum_{i,j,k=\pm1}\!\! ik\,b^{z,\,i,j,2k}.$$

**Consequence for §6:** Step 0 is a port of an already-validated algorithm, not
a derivation. The risk in this project is plumbing. The verification scripts
(`verify_prolong_a.py`, `moment_form.py`) live in this session's scratchpad and
should be committed alongside the C++ self-test; they are the reference
implementation for §7.1–§7.3 and the oracle to bisect against.

---

## 3. Architectural constraint

CarpetX's prolongation is built as a **tensor product of independent 1-D
stencils**: `coeffs1d<CENT,INTP,ORDER,T>` → `interp1d<...>` → `call_stencil_3d`
(`prolongate_3d_rf2_impl.hxx:1220`). The new operator violates that contract in
three ways: the correction is irreducibly 3-D, it reads other components, and
Stage 3 reads the destination array. It therefore **cannot** be a new
`interpolation_t` enum value. It must be a separate `amrex::Interpolater`
subclass, a sibling of `prolongate_3d_rf2`, overriding `interp_arr`.

Since $A_x$, $A_y$, $A_z$ necessarily have different index types, they are
different Cactus groups → different `amrex::MultiFab`s. `Interpolater::interp()`
sees one `FArrayBox`; `interp_arr()` sees three.

---

## 4. What AMReX gives us

The §0 decision — fine data wins, and the prolongation consumes it — is
**already solved structurally in AMReX**, and `FaceDivFree` (the face-centred
$B$ analogue of this operator) is the worked precedent. The pattern, from
`AMReX_FillPatchUtil_I.H:638-780`:

1. `mf_crse_patch[d] = make_mf_crse_patch<MF>(fpc, ncomp, ixtype)` — the coarse
   source, as CarpetX builds today.
2. `mf_refined_patch[d] = make_mf_refined_patch<MF>(fpc, ncomp, ixtype, ratio)`
   — the destination, `convert(refine(coarsen(fpc.ba_fine_patch, ratio), ratio),
   ixtype)`, i.e. **snapped out to whole coarse cells**. That snapping is what
   makes a per-coarse-cell interior closure well defined.
3. It is then **pre-filled from the fine level**:
   `FillPatchSingleLevel(mf_refined_patch[d], time, fmf_time, ...)`. So when
   `interp_arr` is entered, the destination already holds the true fine values
   wherever the fine level has them.
4. `solve_mask[d] = make_mf_crse_mask<iMultiFab>(fpc, ncomp, ixtype, ratio)` —
   one entry per **coarse** index of that index type, on
   `convert(coarsen(fpc.ba_fine_patch, ratio), ixtype)`. Set to 1 everywhere,
   then to 0 through a `CPC` between `coarsen(refined_patch.ba)` and
   `coarsen(fine.ba)` — i.e. 0 exactly where the fine level already has data.
5. The kernel honours it by early return. `facediv_face_interp`
   (`AMReX_Interp_3D_C.H:132`) opens with

   ```cpp
   if (mask) { if (!mask(ci, cj, ck, n)) { return; } }
   ```

   leaving the pre-filled fine value untouched. Its interior phase then reads
   the destination array's face values regardless of provenance — mechanically
   the same thing our §2.5 moment form does.
6. Afterwards `mf[d]->ParallelCopyToGhost(mf_refined_patch[d], ...)` (aliased
   case, which is ours) and a final `FillPatchSingleLevel(*mf[d], nghost, ...)`,
   so the fine level's own data wins again at the end.

Other useful pieces:

* `virtual void Interpolater::interp_arr(...)` — `AMReX_Interpolater.H:89`.
  Default implementation aborts; we override. `interp()` is pure virtual in the
  base, so it must be implemented too (§5.2).
* The **mixed-index-type recipe**, `AMReX_FillPatchUtil_I.H:660-682`: build one
  *cell-centred dummy* `FPinfo` so the box decomposition is index-type
  independent, then convert per direction — same `fpc.dm_patch`, same box count.
* `fine_region` is passed **cell-centred** (`:735-751`); the operator derives
  each component's nodal target box via `amrex::convert`, which grows the high
  end by one per nodal direction.

**What AMReX does not give us:** `FillPatchTwoLevels(Array<MF*,3>...)` and
`InterpFromCoarseLevel(Array<MF*,3>...)` both assert
`mf[d]->ixType().nodeCentered(d)` (`:638`, `:1161`) — hardwired for **face**
arrays. Edge arrays are nodal in the two *transverse* directions and fail it.
There is also no `FillPatchInterp` overload for `Array<MF*,3>`.

That is not a blocker: **CarpetX does not call `FillPatchTwoLevels` or
`InterpFromCoarseLevel` at all** (they appear only in two comments,
`fillpatch.hxx:89` and `:97`). `fillpatch.cxx` reimplements the machinery by
hand in a coroutine style for task overlap, so we control the call site and can
copy AMReX's recipe without inheriting its face-only assertion. **We do have to
copy items 2–4 above, which CarpetX currently does not do** — today it uses
`make_mf_fine_patch` (unsnapped, no pre-fill, no mask). See §5.5.

### 4.1 Follow `FaceDivFree`'s shape

`amrex::FaceDivFree` (`AMReX_Interpolater.H:546`, `AMReX_Interpolater.cpp:1303-1545`,
kernels at `AMReX_Interp_3D_C.H:132` and `:219`) is the face-centred ($B$) Tóth
operator — the same construction this one replaces with edges. It is the closest
thing in the tree to what we are building, and the implementation should mirror
it wherever there is no reason to diverge.

| aspect | `FaceDivFree` (faces, $B$) | this operator (edges, $A$) |
|---|---|---|
| entry point | `interp_arr`; `interp()` is `amrex::Abort("FaceDivFree does not work on a single FArrayBox. Call 'interp_arr' instead.")` (`:1341`) | **same**, with `CCTK_ERROR` |
| loop index | **coarse**: `amrex::convert(c_fine_region, types[d])`, one thread per coarse face writing its $r^2$ fine faces | **same**: one thread per coarse index, writing the 2 fine edges it owns per case |
| phases | 2 launches — one fused 3-direction face phase, then the interior phase over `c_fine_region` | **3** — coarse edges, coarse faces, coarse cells (§5.2). The extra one is difference 3 below |
| mask | `if (mask) { if (!mask(ci,cj,ck,n)) { return; } }` at the top of the coarse-index kernel (`Interp_3D_C.H:140`) | **identical** — per component, null-tolerant, early return |
| interior closure inputs | `facediv_int(ci,cj,ck,n, finearr, ratio, cellSize)` — **no coarse array at all** | **same**, given §2.5's fine-array form |
| interior closure form | 12 expressions, 24 named fine loads, aspect-ratio weights like `dx*(2*dz*dz+dx*dx)/(8*dy*(dz*dz+dx*dx))` | ¼Σ + two moment sums. **No aspect-ratio weights** |
| ratio ≠ 2 | per-thread cached LU solvers, `Vector<std::map<Box,Solver_t>> m_lusolver` | **N/A** — CarpetX is rf2, so only the closed form matters |
| object state | per-thread solver cache, sized in the ctor from `OpenMP::get_max_threads()` | **stateless** |
| `CoarseBox` | `coarsen(fine,ratio).grow(1)` (`:1311`) | `.grow(2)` — see below |
| code layout | kernels are free `AMREX_GPU_HOST_DEVICE` templates in `AMReX_Interp_3D_C.H`; the class only orchestrates | mirror it: kernels in a `..._vecpot_kernels.hxx`, thin class |

**The headline.** `facediv_int` (`AMReX_Interp_3D_C.H:219`) is Tóth
Eqs. (15)–(17) written out in full: 24 fine loads and twelve five-term
expressions carrying weights $\frac{\Delta x(2\Delta z^2+\Delta x^2)}
{8\Delta y(\Delta z^2+\Delta x^2)}$. That is precisely the structure the draft
says the edge formulation removes —

> "Neither the cell aspect ratios … nor the third-order quantities … appear
> anywhere. The elaborate weighting of Eqs. (15)–(17) is exactly what is
> required for the interior fluxes to be the curl of an edge field; once the
> potential is made the primary variable, that structure is automatic and what
> remains is plain averaging."

Our Stage 3 is a mean of eight fine edges plus two signed sums. Same result
(§2.7 checks it against Eq. 15 directly), a fraction of the arithmetic, and no
aspect-ratio algebra to get wrong. If the implementation ever grows an
aspect-ratio weight, that is a bug.

**Where we deliberately differ.**

1. **`CoarseBox` grows by 2, not 1.** `FaceDivFree` interpolates $B$ directly,
   so its tangential slopes reach one coarse cell. We reconstruct $B$ from $A$
   by a discrete curl, which is one derivative further out, so our Stage-2
   slopes reach two coarse cells in $A$ (§5.3).
2. **Three classes of owned object, not one.** `FaceDivFree` has exactly one
   kind of "known" object per direction — a coarse face — so one mask answers
   everything. Edge data has three: coarse edges (Stage 1), coarse faces
   (Stage 2), coarse cells (Stage 3). AMReX's `solve_mask` answers only the
   first directly; the other two must be derived (§8.2). This is the one place
   the correspondence genuinely breaks, and the main design task inside the
   kernel.
3. **An extra class of owning object, hence an extra launch.** Our Stage 1
   (fine edges lying *on* a coarse edge) has no counterpart in `FaceDivFree`,
   because a fine face never lies on a coarse face's boundary the way a fine
   edge lies on a coarse edge. Arithmetically it is free — pure injection — but
   it is a fourth owning index type, so we get three launches where
   `FaceDivFree` gets two, and seven masks where it gets three (§5.2, §8.2).

---

## 5. Design

### 5.1 New files

```
CarpetX/src/prolongate_3d_rf2_vecpot.hxx        // class declaration
CarpetX/src/prolongate_3d_rf2_vecpot.cxx        // kernel + instantiation
CarpetX/src/prolongate_3d_rf2_vecpot_test.cxx   // standalone self-test
```

plus a new test thorn `TestVecPotProlongation/` (§7.5). Add the `.cxx` files to
`CarpetX/src/make.code.defn`.

### 5.2 The Interpolater subclass

```cpp
namespace CarpetX {

// Flux-preserving prolongation for an edge-staggered vector potential.
// Operates on all three components at once; see PROLONG_A.md.
class prolongate_3d_rf2_vecpot final : public amrex::Interpolater {
public:
  ~prolongate_3d_rf2_vecpot() override;

  // NOTE: for this operator CoarseBox takes and returns CELL-CENTRED boxes,
  // matching the interp_arr convention. Grows by 2 coarse cells per direction.
  amrex::Box CoarseBox(const amrex::Box &fine, int ratio) override;
  amrex::Box CoarseBox(const amrex::Box &fine,
                       const amrex::IntVect &ratio) override;

  // Pure virtual in the base, but not meaningful for this operator.
  void interp(...) override;  // -> CCTK_ERROR

  void interp_arr(amrex::Array<amrex::FArrayBox *, dim> const &crse,
                  int crse_comp,
                  amrex::Array<amrex::FArrayBox *, dim> const &fine,
                  int fine_comp, int ncomp,
                  const amrex::Box &fine_region_cc,
                  const amrex::IntVect &ratio,
                  amrex::Array<amrex::IArrayBox *, dim> const &solve_mask,
                  const amrex::Geometry &crse_geom,
                  const amrex::Geometry &fine_geom,
                  amrex::Vector<amrex::Array<amrex::BCRec, dim> > const &bcr,
                  int actual_comp, int actual_state,
                  amrex::RunOn gpu_or_cpu) override;
};

extern prolongate_3d_rf2_vecpot prolongate_vecpot_3d_rf2;

} // namespace CarpetX
```

**Kernel structure — three launches, indexed by the coarse object that owns
each fine edge.** Threads carry a **coarse** index, not a fine one, as in
`FaceDivFree` (§4.1). That is a departure from CarpetX's existing prolongation
kernels, which loop over fine points and derive `icrse = ifine >> 1`
(`prolongate_3d_rf2_impl.hxx:1423`). Coarse indexing makes the mask check a
single early return and lets a slope be computed once and shared by the edges
that use it.

The correct iteration box for each case is `amrex::convert(c_fine_region, T)`
with `T` the index type of the **coarse object that owns** those fine edges.
That is the draft's own decomposition (Sec. V B: "the 54 fine edges split into
24 lying on the twelve coarse edges, 24 lying in the interior of the six coarse
faces, and 6 joining the cell centre to the six face centres"):

| launch | case (for $A_x$) | owning object | iterate over | reads |
|---|---|---|---|---|
| 1 | `(oy,oz)=(0,0)` | coarse x-edge `{cvv}` | `convert(c_fine_region, (C,N,N))` | coarse $A$ |
| 2 | `(1,0)` | coarse z-face `{ccv}` | `convert(c_fine_region, (C,C,N))` | coarse $A$, $B$ |
| 2 | `(0,1)` | coarse y-face `{cvc}` | `convert(c_fine_region, (C,N,C))` | coarse $A$, $B$ |
| 3 | `(1,1)` | coarse cell `{ccc}` | `c_fine_region` | **destination array only** (§2.5) |

Launch 1 fuses the three edge types, launch 2 the three face types — both with
`AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG`, which takes exactly three box/lambda
pairs. Launch 3 is a single `ParallelFor` over coarse cells; within it, form the
cell's surface fine fluxes by discrete curl (§2.3 with $h=\Delta/2$) of the
destination array, then apply §2.5. Launches 1 and 2 are independent of each
other; launch 3 must follow both, because Stage 3 at cell `(ix,iy,iz)` reads
Stage-2 output at `(ix,iy,iz+1)`.

> **Do not give every case the same iteration box.** On a patch of $2^3$ coarse
> cells the $A_x$ destination holds 100 fine edges, partitioned **exactly** as
> 36 (edge-owned) + 24 + 24 (face-owned) + 16 (cell-owned) — verified disjoint
> and covering for all three components. A single nodal box for all four cases
> writes out of bounds at the top nodal index, where only the `off == 0` cases
> exist. The owning-type boxes are self-clipping; nothing else is.

Each launch writes only where `solve_mask` says "solve"; where it says "known"
the destination already holds the fine value and is left alone. The mask for
each launch is indexed on that launch's owning type — see §8.2.

**One templated kernel, not twelve transcriptions.** §2.4 spells the cases out
longhand for readability, but do **not** implement them that way. The rule of
§2.1 is already parameterized by the edge direction $d$ and the two transverse
parities $p_b$, and in discrete form it is

$$a_d \;=\; \big\langle A_d\big\rangle
\;+\; \sum_{b\neq d} p_b\,\varepsilon_{dbc}\,\frac{\Delta_b}{16}\,
\Big\langle B^c_{+1_b} - B^c_{-1_b}\Big\rangle$$

where both $\langle\cdot\rangle$ are means over the directions that were
averaged — all $p=1$ directions for the $A$ term, and the *other* transverse
direction for the $B$ term (which is what turns the $1/16$ into the $1/32$ of
the cell-interior case). Write it once as

```cpp
template <int D, int P1, int P2> struct vecpot_stencil;  // 12 instantiations
```

and generate the instantiations with a `constexpr_for`, as AMReX does for its
refinement-ratio cases (`AMReX_Interpolater.cpp:1507`). The `(P1,P2) == (1,1)`
specialization uses §2.5's moment form instead; the other nine use the formula
above.

Two reasons this matters beyond brevity. In the six mixed instantiations — the
only ones with a $B$ term, since `(0,0)` has none and `(1,1)` uses the moment
form — the $B$ lookup index is just the owning object's index with the
$b$-direction shifted $\pm1$, so **the index types line up automatically**:
when $p_b = 1$ the owning type is CELL in $b$, which is exactly what $B^c$ is,
and it is NODE in the remaining transverse direction, again matching $B^c$. The
template therefore cannot emit a type-inconsistent lookup; a hand transcription
can. And it makes the §4.1 invariant — *an
aspect-ratio weight anywhere is a bug* — checkable by reading one function
rather than twelve cyclic permutations, which is precisely where a sign error
would hide.

Other notes:

* No template parameters. The centering triple is fixed; the field order is
  fixed at 2.
* `ratio` must be `(2,2,2)`; assert it.
* Keep `ncomp` general (a group may hold several independent vector potentials)
  but require the same `ncomp` on all three components. **The coupling is
  component-wise**: component `c` of $A_x$ reads component `c` of $A_y$, $A_z$,
  never a different one. Note `solve_mask` is also per-component
  (`mask(ci,cj,ck,n)`).
* `solve_mask[d]` may legitimately be `nullptr` (no fine data anywhere, e.g. the
  regrid path). Guard exactly as `facediv_face_interp` does; a null mask means
  "solve everything".
* The `prolongate_per_group` parameter (`param.ccl:463`) is meaningless here —
  this operator always processes the triple together. Ignore it, and say so in
  a comment so nobody wires it up.

### 5.3 Required ghosts and `CoarseBox`

Reach of the `(oy,oz) = (1,1)` case **when prolonging $A_x$**, relative to the
base coarse index `(ix,iy,iz)`:

| needs | x | y | z |
|---|---|---|---|
| `Ax` | `ix` | `iy-1 .. iy+2` | `iz-1 .. iz+2` |
| `Ay` | `ix .. ix+1` | `iy-1 .. iy+1` | `iz .. iz+1` |
| `Az` | `ix .. ix+1` | `iy .. iy+1` | `iz-1 .. iz+1` |

and cyclically for $A_y$, $A_z$. **`required_ghosts = 2` uniformly** is the
envelope over all three targets, and is what `CoarseBox` should implement —
CarpetX grows symmetrically (`prolongate_3d_rf2_impl.hxx:1262`:
`amrex::coarsen(fine, 2)` then `crse.grow(d, required_ghosts[d])`).

The moment form of §2.5 does not widen this: it reads the *destination* array
inside the current coarse cell, not more coarse data. Sweep 2's coarse-$B$
reach is what sets the envelope.

> **Why 2 and not 3.** In a nodal direction the fine node indices of a target
> region run one past `coarsen(fine_region_cc)` on the high side, so one might
> expect `2 + 1 = 3`. But that extra index is always reached with `off == 0`,
> whose stencil has no `+2` reach — the `+2` comes only from the `off == 1`
> cases, whose base index stays inside the coarsened box. Both ends are covered
> by 2. Worth a comment in the code, or someone will "fix" it.

For comparison, `ddf` o1 for `c100` needs `{0,1,1}`.

### 5.4 Group-triple discovery

`get_interpolator()` (`driver.cxx:671`) keys on
`(prolongation_type, prolongation_order, indextype)` and returns one
`Interpolater*` per group. The triple must be recognised *before* that.

Exact precedent: `get_group_fluxes()` (`driver.cxx:445`) parses
`TAGS='fluxes="fluxx fluxy fluxz"'` into three group indices, resolving bare
names against the group's own implementation and asserting there are exactly
`dim` of them. Copy it almost verbatim:

```
CCTK_REAL Avec_x TYPE=gf CENTERING={cvv}
  TAGS='parities={-1 +1 +1} prolongation_type="vecpot"
        vector_potential="Avec_x Avec_y Avec_z"'
```

Every member carries the *same* tag naming all three, so the triple is
discoverable from any one of them. The new tag is independent of the existing
`fluxes` (reflux) tag; both may be present.

Validation in `get_group_vector_potential()`:

* exactly `dim` names, all valid `CCTK_GF` groups;
* `indextype` of entry `d` is the edge type for direction `d`
  (`indextype[e] == (e == d)` in CarpetX convention);
* all three have the same `numvars` and the same `nghostzones`;
* all three carry the same tag value (so the relation is symmetric);
* `prolongation_type == "vecpot"` on all three.

Store on `GroupData` (`driver.hxx:421`):

```cpp
std::array<int, dim> vecpot_groups;   // group indices, -1 if not a triple
int vecpot_direction;                 // 0/1/2, or -1
```

`GroupData::GroupData` already computes `indextype` before calling
`get_interpolator` (`driver.cxx:943-946`), so the tag can be parsed there.
Add `"vecpot"` to `param.ccl`'s `prolongation_type` keyword list.

All three groups on a level share the same base `BoxArray` and
`DistributionMapping` — `GroupData` is constructed from the level's `ba`/`dm`
and only converts the index type (`driver.cxx:991-999`). This is what makes the
`interp_arr` route viable; assert it at the call site. Multipatch runs resolve
the triple within one `patchdata`, never across patches.

### 5.5 `FillPatchInterp_arr` and the buffer change

This is the largest single piece of work, because CarpetX's buffers must gain
the refined-patch/pre-fill/mask machinery of §4 items 2–4.

Add to `fillpatch.cxx`:

```cpp
void FillPatchInterp_arr(
    std::array<amrex::MultiFab *, dim> const &fine_refined_patch,
    std::array<const amrex::MultiFab *, dim> const &crse_patch,
    std::array<const amrex::iMultiFab *, dim> const &solve_mask,
    int ncomp, const amrex::Geometry &cgeom, const amrex::Geometry &fgeom,
    const amrex::Box &fdomain_g_cc, const amrex::IntVect &ratio,
    amrex::Interpolater *mapper,
    const std::array<amrex::Vector<amrex::BCRec>, dim> &bcrecs);
```

looping `MFIter` over the cell-centred dummy decomposition, assembling the
three `FArrayBox*` / `IArrayBox*` triples, converting the target box to
cell-centred, and calling `mapper->interp_arr(...)` once per box. Model on
`AMReX_FillPatchUtil_I.H:717-755`.

Changes to the surrounding code:

* Replace `make_mf_fine_patch` with `make_mf_refined_patch` (snapped to whole
  coarse cells) for the vecpot path. **Non-negotiable**: the interior closure
  of §2.5 is defined per coarse cell and needs all six faces of that cell in
  the destination FAB.
* Pre-fill the refined patch from the fine `MultiFab` before interpolating.
  In `FillPatch_Prolongate` the fine source is `mfab` itself, so this is a
  `ParallelCopy` from `mfab` into the refined patch — and it must complete
  before `interp_arr` runs.
* Build `solve_mask` per direction via the `CPC` recipe of §4 item 4.
* Sequence the existing `FillBoundary` on `mfab` **before** the pre-fill, so
  the pre-fill sees valid fine ghost data too.

> **BC ordering constraint.** `Prolongate_Start` currently does
> `coarsegroupdata.apply_boundary_conditions(crse_patch)` immediately before
> `FillPatchInterp` (`fillpatch.cxx:167-175`). In the `_arr` version, **all
> three** coarse patches must have their boundary conditions applied before
> **any** of them is read, because prolonging $A_x$ reads $A_y$ and $A_z$.
> Applying them lazily per component inside the loop silently reads unfilled
> outer-boundary data. Apply all three up front, then interpolate.

Likewise the three coarse-patch `ParallelCopy_nowait` calls should be issued
together and joined together, to keep the existing overlap behaviour.

### 5.6 Call sites

Six sites across five paths (all pass `groupdata.interpolator`):

| path | site | notes |
|---|---|---|
| `Sync` / prolongate ghosts | `sync_restrict.cxx:449`, `:635` (`FillPatch_ProlongateGhosts`) | highest frequency; has fine data, so the mask matters most here |
| `ProlongateRestrictedGFs` | `sync_restrict.cxx:225` (`FillPatch_ProlongateOnly`) | shares `Prolongate_Start`/`_Finish` with the above |
| `MakeNewLevelFromCoarse` | `driver.cxx:1608` (`FillPatch_NewLevel`) | no fine data exists ⇒ null mask, no pre-fill; simplest case |
| `RemakeLevel` | `driver.cxx:1776` (`FillPatch_RemakeLevel`) | old fine data exists and is copied over the result afterwards |
| subcycling RK boundary | `subcycling.cxx:495` (`Prolongate_Start`/`_Finish`) | persistent buffers, see `driver.hxx:452-467` |

Each needs an `_arr` sibling taking `std::array<...,dim>` of
`GroupData`/`MultiFab`. The per-group loops in `sync_restrict.cxx` and
`driver.cxx` must skip the two non-leading members of a triple and dispatch the
whole triple once, on the leading member (`vecpot_direction == 0`).

### 5.7 Alternatives considered and rejected

**A single fused launch over a unified lattice.** All staggered objects —
coarse and fine, vertices, edges, faces, cell centres — are points of *one*
integer lattice in units of half a fine cell, with the parity vector giving the
index type. That is the draft's own convention (Sec. II: coarse centre at 0,
coarse faces at $\pm2$, fine centres at $\pm1$), and it invites a single loop
over coarse cells writing all 24 fine edges each one owns, computing every
edge's displacement from the cell's base vertex.

Rejected, for one decisive reason and three supporting ones. **Decisive:**
Stage 3 in coarse cell $(i,j,k)$ reads the fine fluxes on that cell's *upper*
faces, which under half-open ownership belong to cells $(i{+}1,j,k)$ etc. —
concurrent threads. Fusing would force each thread to recompute those faces
rather than read them, which roughly doubles the Stage-2 arithmetic (a face's
flux needs all twelve of its fine edges) *and* would still have to consult that
face's mask and fall back to the pre-filled fine data, so the copy path is not
avoided, only complicated. **Supporting:** AMReX stores each index type in its
own `FArrayBox`, so a unified index costs a shift-and-divide per access; one
loop would consult four differently-indexed masks per thread instead of one per
launch; and 24 formulas in a single thread body is register pressure and warp
divergence, where the per-owning-type launches are uniform — which is exactly
why `FaceDivFree` puts one thread on one coarse face.

The idea is right about the *formula* and wrong only about the *loop*: see the
templated kernel in §5.2, which is the same parameterization applied at compile
time. It is also worth keeping the draft's doubled index convention as the
internal notation of the self-test, so the C++ can be diffed against the paper
line by line.

**A new `interpolation_t` value inside `prolongate_3d_rf2`.** Impossible: that
framework is a tensor product of independent 1-D stencils, and this operator is
irreducibly 3-D, reads other components, and reads its own destination. See §3.

**A post-prolongation correction pass.** A thorn-level routine scheduled after
`SYNC` that adds the cross-component correction to already-prolonged ghosts.
Rejected on three counts: nothing exposes which fine edges were prolonged
versus copied; Cactus hands a thorn one refinement level at a time, so the
coarse data is not reachable; and the correction would run after the ghost data
had already been consumed.

**Dropping the cross-component terms.** Keeping only
$-(\partial_y^2+\partial_z^2)A_x$ turns the average into a higher-order
interpolation of $A_x$ alone — divergence-free but not flux-conserving, the
approach Sec. V C explicitly rejects. See §2.6.

**A standalone thorn — deferred, not rejected.** See §0. The blocker is that
`fillpatch.cxx` never calls `interp_arr`; registration alone is already
reachable (`GroupData::interpolator` is a public member of an exported struct).
The CarpetX-side work is generic — a face-centred $B$ operator would reuse it
unchanged — so the natural split is the Array path plus a registry in CarpetX,
and the stencil in a thorn. Revisit after Step 3, when there is a real consumer
to shape the registry around.

---

## 6. Staged plan

### 6.0 Ordering principle

Ordered by **plumbing complexity and testability**, not by math coverage: §2.7
has validated the math, and §2.4 shows every path needs all of it. Each step
should build and pass the existing testsuite
(`./agent_scripts/build.sh && ./agent_scripts/test.sh`).

Note the deliberate inversion relative to the obvious order: **`FillPatch_NewLevel`
(regrid) is the simplest path**, because no fine data exists there — null mask,
no pre-fill. It is the natural first integration target even though it is not
the most frequently executed.

### Step 0 — kernel and offline self-test

* `prolongate_3d_rf2_vecpot.{hxx,cxx}`: the class, `CoarseBox`, and the
  three-launch `interp_arr` kernel of §5.2, with §2.5's moment form for launch 3.
* `prolongate_3d_rf2_vecpot_test.cxx`: a C++ port of the numpy verification of
  §2.7, run once behind a `static` guard like the `test_interp1d` instances at
  `prolongate_3d_rf2_impl.hxx:1324`. Include a **masked** variant: pre-fill some
  faces with independent "fine" data, mask them, and check that (a) they are not
  overwritten and (b) the interior edges come out as the moment form predicts.
* Not yet reachable from `get_interpolator`.

**Deliverable:** a kernel reproducing the numpy reference to machine precision,
with correct masked behaviour. A mismatch is a transcription bug; bisect against
the scripts.

### Step 1 — group-triple discovery

* `get_group_vector_potential()` in `driver.cxx`, per §5.4.
* `vecpot_groups` / `vecpot_direction` on `GroupData`.
* `"vecpot"` in `param.ccl` and in `get_interpolator`.
* Nothing calls `interp_arr` yet; a `"vecpot"` group errors out loudly if
  reached. Verify tag parsing with a dummy thorn.

### Step 2 — `FillPatch_NewLevel_arr` (regrid from coarse)

* `FillPatchInterp_arr` (§5.5) in its simplest form: null mask, no pre-fill.
* Dispatch from `driver.cxx:1608`.
* Validity bookkeeping: the loop starting at `driver.cxx:1551` sets
  `groupdata.valid` per group; mark the triple valid together and make the
  `error_if_invalid` pre-checks cover all three coarse groups.

**Deliverable:** creating a new fine level from coarse produces a
flux-conserving, gauge-clean fine $A$. Test: pure-gauge data, assert $B\equiv0$.

### Step 3 — refined patch, pre-fill and mask; the ghost paths

* The buffer change of §5.5: `make_mf_refined_patch`, pre-fill from the fine
  level, `solve_mask` via `CPC`.
* `FillPatch_Prolongate_arr` + `Prolongate_Start_arr` / `Prolongate_Finish_arr`,
  preserving the coroutine structure and the BC ordering constraint.
* Dispatch from `sync_restrict.cxx:449`, `:635`, `:225`.
* Resolve the partial-triple question (§8.7) in `sync_filter_groups`
  (`sync_restrict.cxx:127`-ish), the single funnel for all four `Sync` entry
  points.
* This is the step that creates mixed coarse faces (§8.3). Add the seam
  conservation assertion (§7.8) with it, not after it.

**Deliverable:** a static refinement hierarchy prolongs ghosts correctly, fine
data is demonstrably never overwritten, and flux conservation holds across the
seam.

### Step 4 — `FillPatch_RemakeLevel_arr`

* `driver.cxx:1776`. Reuses everything from Step 3; the wrinkle is that old
  fine data is copied over the result afterwards.

### Step 5 — subcycling

* `subcycling.cxx:495` uses `Prolongate_Start`/`_Finish` with *persistent*
  buffers (`driver.hxx:452-467`), sized from one `FPinfo` per group; a triple
  needs three, sharing one cell-centred `FPinfo`, plus the masks.
* Last. Most intricate path, orthogonal to the math.

### Step 6 — companions

* **Gauge scalar $\Phi$.** Sec. V C's third remark: prolonging $\Phi$ (or
  $\sqrt\gamma\,\Phi$) with the *same* multilinear stencil makes the prolonged
  pair the exact gauge transform of the coarse pair, so a regrid injects no
  gauge mismatch. §2.7 confirms our $P_\lambda$ is exactly trilinear, and
  trilinear on a `{vvv}` group is already `prolongate_ddf_3d_rf2_c000_o1` —
  a **documentation change only**: tag the $\Phi$ group with
  `prolongation_type="ddf"`, `prolongation_order=1`.
* **High-order Stage-1 split**, Eq. (32):
  `a = Ax(ix,iy,iz) + s*(1/8)*( Ax(ix+1,iy,iz) - Ax(ix-1,iy,iz) )` with
  `s = (ox == 0 ? -1 : +1)`. Pure gauge — changes $A$, not $B$ — so it must be
  added *together with* the compensating terms in the face and cell formulas
  (Sec. V E, last paragraph) or not at all. Parameter, default off (§0).
  Note the moment form of §2.5 absorbs the face-formula compensation
  automatically; only the cell formula's $\langle A\rangle$ term needs revisiting.

---

## 7. Testing

### 7.1 Gauge test (the headline)

Coarse `Ax(i,j,k) = (lam(i+1,j,k) - lam(i,j,k)) / dX` and cyclically, for a
random vertex scalar `lam`. Then every prolonged fine flux must be zero to
machine precision, and (stronger, Eq. 31) the prolonged $A$ must equal the
discrete gradient of the *trilinear* interpolant of `lam`. Both confirmed in
§2.7. This is what distinguishes the operator from generic Lagrange
interpolation of $A_i$.

### 7.2 Flux conservation

Random coarse $A$, non-unit aspect ratio. Assert Eq. (7) (restriction
compatibility) and Eqs. (11)/(15) (Tóth equivalence), using the component-form
references in §2.7.

### 7.3 Divergence

Fine $\nabla\cdot B = 0$ cell by cell. Identically true for any edge
representation (Eq. 6), so it is a cheap assertion that the test harness's own
curl bookkeeping is right.

### 7.4 Convergence

Smooth analytic $A$; $\|B_{\rm fine} - B_{\rm exact}\|$ → second order under
resolution doubling. $A$ itself converges at second order with an
$O(\Delta^2)$ *gauge* error (§2.6): the residual is
$(\Delta Y^2/8)\,\partial_x\partial_y A_y$, a pure gradient, invisible in $B$.

### 7.5 Integration test thorn

`TestVecPotProlongation/`, modelled on `StaggeredWaveToyX`:

* three groups `{cvv}`, `{vcv}`, `{vvc}` with the right `parities` and the
  `vector_potential` tag; a `{vvv}` gauge scalar;
* initial data 1: pure gauge — regrid, assert $B\equiv0$;
* initial data 2: a smooth force-free field — flux conservation across the
  refinement boundary over several regrids;
* **a fine-data-preservation test**: perturb the fine level's $A$ in its valid
  region, run a ghost fill, assert the valid region is bit-identical;
* **a symmetry variant** exercising a reflection boundary (§8.4);
* reference output per `agent_docs/generating-reference-output.md`.

### 7.6 Regression guard

With the correction terms forced to zero, the operator must reproduce
`prolongate_ddf_3d_rf2_c100_o1` (and cyclic) bit for bit (§2.6). Keep as a
compile-time-switchable debug mode; it isolates plumbing bugs from math bugs.

### 7.7 Restriction round-trip

`restrict ∘ prolong` on a coarse edge must be the identity. Cheap, and it
guards the interaction with `average_down_edges` (§8.5).

### 7.8 Seam conservation

On every coarse face that straddles the edge of `ba_fine_patch` — i.e. one that
receives a mix of copied and prolonged edges (§8.3) — assert that the four fine
fluxes sum to the coarse flux. This is the assertion that turns §8.8's
restriction-consistency precondition from an assumption into a checked
invariant, and it is the one test that would catch a `Restrict` ordering
mistake or a stale coarse level. Run it in the ghost-fill path (Step 3) and
again under subcycling (Step 5), where it is expected to pick up the
time-interpolation error.

---

## 8. Risks and open questions

1. **Coarse validity of all three components.** The correction reads $A_y$,
   $A_z$ while prolonging $A_x$. In `Sync`, groups are processed one at a time;
   the triple must be dispatched as a unit *and* the coarse level must be valid
   and synced for all three at that moment. CarpetX's `valid`/`why_valid_t`
   bookkeeping catches violations only at runtime and only with
   `poison_undefined_arrays` on. Add explicit `error_if_invalid` checks for all
   three coarse groups at the head of the `_arr` paths.

2. **Seven masks, one per owning object.** Provenance has to be answered at
   three granularities — is this coarse *edge* / *face* / *cell* already covered
   by fine data? This is less work than it first looks:
   `make_mf_crse_mask(fpc, ncomp, idx_type, ratio)`
   (`AMReX_FillPatchUtil_I.H:432`) takes an **arbitrary** `idx_type` and returns
   a mask on `convert(coarsen(fpc.ba_fine_patch, ratio), idx_type)`. So build
   three edge masks, three face masks and one cell mask with the same helper,
   each matching its launch's iteration box (§5.2). No derivation from the edge
   mask, and no new machinery.

   What is left is only the **coverage predicate**: AMReX fills the mask by a
   `CPC` between `coarsen(refined_patch.ba)` and `coarsen(fine.ba)`
   (`AMReX_FillPatchUtil_I.H:703-714`), and for the face and cell types the
   corresponding "fine" box array has to be built by conversion rather than
   taken from an existing `MultiFab`. Getting the boundary convention right
   there — a coarse face on the seam is covered only if all twelve of its fine
   edges are — is the real design task. Note Stage 3's *inputs* need no mask at
   all: §2.5 reads the destination array, which is correct either way.

3. **The seam: mixed coarse faces are real, and benign.** A coarse face
   straddling the edge of `ba_fine_patch` genuinely receives a *mix* — some
   bounding edges copied from the fine level, the rest prolonged from coarse.

   *Why the mix happens.* The three components have different index types, so
   "covered by fine data" is a different set for each. Take a fine level filling
   the half-space $x<0$ and a coarse cell $C=[0,1]^3$ entirely in the ghost
   region. Look at $C$'s z-face at $z=0$ and its four bounding coarse edges:

   | bounding edge | index type in x | covered by fine data? |
   |---|---|---|
   | $A_x$ at cell-x 0 (×2, at node-y 0 and 1) | CELL | no — cell 0 is ghost |
   | $A_y$ at node-x 0 | NODE | **yes** — a fine box over cells $[-4,-1]$ owns nodes $[-4,0]$ |
   | $A_y$ at node-x 1 | NODE | no |

   Nodal data sticks out one index further than cell data, so the fine valid
   region for $A_y$ reaches one node past the fine valid region for $A_x$. The
   coarse face straddles that offset. This is generic at every refinement
   boundary, not a corner case.

   *Why it does not break conservation.* The sum of the four fine fluxes
   through a coarse face equals the circulation of the **eight perimeter** fine
   edges — the four mid-line edges are interior to the face and cancel between
   adjacent quadrants. Those eight are the two halves of each of the four
   bounding coarse edges, and conservation needs only that each pair sums to the
   coarse edge value. A prolonged pair satisfies this by construction (Stage 1
   splits equally). A **copied** pair satisfies it too, because
   `average_down_edges` *is* $\tfrac12(a_1+a_2)$ — the restriction that produced
   the coarse value is exactly the sum rule (C2). Verified numerically (§2.7):
   copying one coarse edge's halves changes individual fine fluxes by O(2) while
   leaving every affected coarse-face sum exact to $5.6\times10^{-17}$, and
   $\nabla\cdot B = 0$ throughout.

   *What this actually costs.* Nothing, given the precondition in item 8: the
   coarse data at and under the fine boundary must be the restriction of the
   fine data. That is the ordinary AMR requirement, not something new to this
   operator — but here it is load-bearing for flux conservation, so assert it
   (§7.8) rather than assume it.

4. **Symmetry boundaries and parities.** The correction mixes components, so
   parity bookkeeping must be consistent across it: under a reflection $A_x$
   carries `parities={-1 +1 +1}` and the induced parity of $B^z$ is a product of
   the $A_x$/$A_y$ parities. CarpetX applies parities per component in
   `apply_boundary_conditions`, with no knowledge of the coupling. If the coarse
   patch's boundary fill is parity-consistent the discrete curl comes out right
   automatically, but this has not been checked. **Add a reflection variant to
   the test thorn (§7.5)** — a sign error here is invisible in periodic tests.

5. **Restriction is already correct — verified, don't change.**
   `sync_restrict.cxx:1236` dispatches rank-1 groups to
   `amrex::average_down_edges`. Its kernel `amrex_avgdown_edges`
   (`AMReX_MultiFabUtil_3D_C.H:223`) computes, for an x-edge,
   `crse(i,j,k) = (1/facx) * sum_{iref} fine(i*facx+iref, j*facy, k*facz)` —
   the arithmetic mean of the two fine edges at the coincident transverse
   nodes, which in component form is exactly condition (C2). Keep §7.7 as a
   cheap guard. Note the `#warning "TODO: Allow different restriction
   operators, and ensure this is conservative"` at `sync_restrict.cxx:1224`.

6. **This operator enforces $D = 0$.** Unlike Tóth, which preserves an arbitrary
   coarse divergence, condition (C1) presupposes $D = 0$ — a vector potential
   cannot represent a nonzero coarse divergence. Correct and harmless for MHD,
   but `"vecpot"` is **not** a general-purpose staggered-field prolongation.
   Assert the edge index types and reject anything else.

7. **Partial-triple SYNC.** The `Sync` paths take the group list from the
   schedule (`sync_filter_groups`, `sync_restrict.cxx:127`-ish, the single
   funnel for all four entry points). Nothing stops a thorn from
   `SYNC: Avec_x` alone, which the coupled operator cannot honour. Decide:
   silently expand to the full triple (recommended, with a one-time warning) or
   error out. Expanding is safe but means a SYNC statement does more than it
   says.

8. **Restriction consistency.** One thing still depends on the coarse data at
   and under the fine boundary being the restriction of the fine data: flux
   conservation across a mixed face (item 3). (The $\langle A\rangle$ term no
   longer does — §2.5's fine-array form removed that dependence, following
   `facediv_int`; see §4.1.) It is correct only after `Restrict` has run. Two
   places to check:

   * the schedule ordering in the `Sync` paths — `Restrict` then
     `ProlongateRestrictedGFs` (`sync_restrict.cxx:307`, `:798`) already has
     the right shape, but confirm it for the vecpot triple;
   * the **subcycling** path (Step 5), where the fine level is deliberately
     misaligned in time with the coarse and the coarse patch is a time-blend of
     two snapshots (`FillPatch_Prolongate`'s `w_new`/`cmfab_old`). A blended
     coarse edge is not the restriction of either fine state. Expect a
     conservation defect of order the time-interpolation error here, and measure
     it before assuming it is acceptable.

9. **GPU.** The three-launch structure follows `FaceDivFree`'s: two fused
   3-direction launches then one over coarse cells (§4.1, §5.2). Use the same idioms —
   `AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG` for the fused phase,
   `AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG` for the interior — and honour the
   `runon` argument rather than assuming a target. Register pressure is higher than the existing
   kernels (three `Array4`s plus curl reconstructions plus the mask). Write it
   device-ready from the start; do not retrofit.

10. **`CoarseBox` index-type convention.** For the scalar operators `CoarseBox`
    takes a box in the group's own index type; for the Array path AMReX passes
    cell-centred boxes. Mixing these silently gives an off-by-one coarse patch.
    Document it on the method and assert `fine.ixType().cellCentered()`.

11. **Six call sites plus a buffer change is a large surface.** The staging in
    §6 exists to keep each step reviewable. Resist adding all `_arr` variants at
    once.

---

## 9. Deliberately deferred

* **Limiters** (§0). Eq. (9) mentions limiting "where the scheme requires it".
  Unlimited centred slopes first; limiting couples the three components'
  stencils further and interacts with the ENO/minmod machinery in a way that
  needs its own design. Revisit once §7.4 convergence is established.
* Higher-order prolongation of the *field* (Sec. V E, last paragraph): requires
  replacing the second-order surface prolongation Eq. (11) and hence the
  corrections in Eq. (30) by higher-order counterparts. Out of scope in the
  draft too.
* Curvilinear / GR densitisation. The draft notes (end of Sec. II) that the
  same manipulations apply verbatim to $\sqrt\gamma A_i$ on edges and
  $\sqrt\gamma B^i$ on faces. If the evolved variable is already densitised,
  nothing changes; if not, the operator is being applied to the wrong variable.
  Decide and document which convention the tagged groups carry.
