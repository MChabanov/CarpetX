# 2D Plane Output for Silo and openPMD

Axis-aligned 2D slabs (xy / xz / yz) at chosen world-coordinate elevations,
written through the Silo and openPMD code paths. File-per-plane, AMR-aware
(per-level snap), Cartesian patches only. Implemented in
`CarpetX/src/io_planes.{hxx,cxx}` (spec parsing, tag formatting, per-level
snap, 3D-FAB slab extraction), `io_silo_planes.{hxx,cxx}` (`OutputSiloPlanes`),
and `io_openpmd.{hxx,cxx}` (`OutputOpenPMDPlanes`); dispatched from
`io.cxx::OutputGH`.

## Parameters

| Parameter | Default | Purpose |
|---|---|---|
| `out_silo_planes` | `""` | Specs, e.g. `"xy:0.0, xy:10.0, xz:-5.0"` |
| `out_silo_plane_vars` | `""` | Variable regex |
| `out_silo_planes_every` | `-1` | Cadence (-1 ⇒ `out_silo_every`) |
| `out_openpmd_planes` / `_plane_vars` / `_planes_every` | — | Same, for openPMD |
| `out_planes_int_precision` | `4` | Min integer digits in tag |
| `out_planes_frac_precision` | `3` | Fractional digits in tag |

Spec syntax: `<axes>:<elevation>`, comma-separated; `<elevation>` is the
world coordinate along the normal axis.

## Snapping & tags

Elevations are rounded to `frac_precision` digits at parse (one-time warning on
overflow); the rounded value is both the tag and the snap target. Per level,
snap to the nearest vertex (VC: `x0 + i·dx`) or cell centre (CC: `x0 +
(i+½)·dx`). A plane outside every Cartesian (patch, level) along its normal —
or a non-Cartesian patch — is skipped with a one-time warning and **no file is
written**.

Tag: `<plane>_<axis>_<sign><int>p<frac>`, e.g. `xy_z_pos0012p500` (z=+12.5),
`xz_y_neg0003p000` (y=−3). Survives Silo `legalize_name` and HDF5/ADIOS2 rules.
Silo mesh names get a 2-char in-plane centering suffix (`cv`/`vc`) for rank-1
staggered groups (per-group mesh); rank-0/2 groups share a mesh (no suffix).

## Output files

```
Silo:
  <out_dir>/<sim>.<tag>.it<N>.silo_planes.dir/<sim>.<tag>.it<N>.p<NNNNNN>.silo
  <out_dir>/<sim>.<tag>.it<N>.silo            metafile (multimesh/var + DBmrgtree)
  <out_dir>/<sim>.<tag>.silo_planes.visit     VisIt index
openPMD:
  <out_dir>/<sim>.<tag>.it<N>.<ext>           one file per iteration (.bp5/.h5/…)
  <out_dir>/<sim>.<tag>.openpmd.visit         VisIt index
```

- **Silo** includes ghost cells (`DBOPT_LO/HI_OFFSET` markers); the metafile
  carries an AMR `DBmrgtree` (`amr_decomp/levels/patches` + `lvlRatios`/
  `ijkExts`/`xyzExts`/`rank` mrgvars) so VisIt shadows coarse data with finer
  levels where they overlap. The mrgtree (and the multimesh's
  `DBOPT_MRGTREE_NAME`) is written **only when at least one coarse→fine
  parent-child relationship exists** (`emit_amr` in `io_silo_planes.cxx`);
  single-level output, or a plane that crosses only one level, is a plain
  multimesh — see *mrgtree pitfalls* below.
- **openPMD** is interior-only (ghosts stripped), MPI-collective `storeChunk`,
  one 2D mesh per (group, patch, level); cell-centring via per-component
  `setPosition`. Each iteration also carries AMR-hierarchy attributes
  (`numPatches`, `patchSuffixes`, per-(patch,level) `chunkInfo`/`iteration_*`)
  and plane keys `planeTag`/`planeNormalAxis`/`planeElevation`. Empty
  `chunkInfo` is omitted (a zero-length attribute is unreadable by openpmd-api).
  A `(patch, level, group)` whose boxes don't reach the plane is skipped
  entirely (no mesh record): `snap_to_grid_index` only checks the level's full
  refined domain, so a finer level can pass it without actually covering the
  plane, and defining a mesh with no stored chunk makes openpmd-api emit "No
  extent found" read warnings. The skip uses the replicated BoxArray so all
  ranks agree (mesh creation is collective).

## Viewing in VisIt

VisIt reads **Silo** natively: open the `.silo_planes.visit` index and
pseudocolor — the mrgtree gives correct AMR level shadowing. openPMD output
can be written as HDF5 (`CarpetX::openpmd_format = "HDF5"` → `.h5`; nothing in
the writer prevents it, as long as the linked openPMD-api has the HDF5
backend), but direct VisIt plotting of openPMD requires VisIt's openPMD reader
plugin, which is not in every stock build — so Silo is the reliable VisIt path.

### mrgtree pitfalls (empty child map; 3D name mismatch)

The Silo AMR `DBmrgtree` has two sharp edges, both found while debugging a VisIt
`SIGSEGV` (`avtSiloFileFormat::GetMesh` → `DBFreeGroupelmap` → `cfree`) on
single-level plane output.

1. **All-empty child map crashes VisIt.** The child map (`_wmrgtree_chldMaps`,
   a `DBPutGroupelmap`) has one segment per component listing that component's
   finer-level children. With only one level — or a plane that intersects a
   single level — *every* segment is empty. Silo then stores no segment-data
   array; on read-back `segment_data == NULL` while `num_segments > 0`, so VisIt
   dereferences a NULL array when freeing the map and crashes. The numeric test
   never sees this: `plane_readers/silo_reader.py` reads the leaf per-box
   `.silo` files directly and never opens the metafile multimesh or walks the
   mrgtree, so CI stays green while VisIt dies. **Fix applied:** the plane writer
   emits the mrgtree only when `total_children > 0`, else a plain multimesh
   (`emit_amr`, `io_silo_planes.cxx`).

2. **The 3D writer (`io_silo.cxx`) hides the same hazard behind a bug.** It names
   the mrgtree object `"mrgTree"` (`DBPutMrgtree`) but sets the multimesh's
   `DBOPT_MRGTREE_NAME` to `"mrgtree"` — a case mismatch. Silo/HDF5 names are
   case-sensitive, so VisIt's `DBGetMrgtree` lookup fails, the tree is silently
   ignored, and the multimesh is read as plain. That is why single-level 3D Silo
   output opens fine while the (correctly-wired) plane output crashed. The cost
   is that **3D AMR Silo output has never actually received the intended mrgtree
   level-shadowing in VisIt** — the feature is dead code on the read side.

   *Future fix:* align the two names in `io_silo.cxx` so VisIt loads the 3D
   mrgtree — but this **must** land together with the same `emit_amr`/`nlevels >
   1` guard as the plane writer, otherwise single-level 3D output will start
   crashing exactly as the planes did. Consider factoring the shared mrgtree
   construction (group maps, region arrays, ratio/extent mrgvars, the
   `total_children > 0` guard) out of both writers so the two paths cannot drift
   again. A cheaper CI safety net: have the verifier also open the metafile
   multimesh (not just the leaf files) so a tree VisIt would reject fails the
   test.

## Testing (thorn `TestPlanes`)

`TestPlanes` is a rigorous, CI-integrated test suite. It declares one grid
function per index-type centering (the eight V/C combinations) and fills each,
on the interior, with the same analytic axis-distinct field

    f(x, y, z) = x + 100 y + 10000 z

via `WRITES interior` + `SYNC` at `initial` and `postregrid`. Because `f` is
linear and `CCTK_INITIAL` re-runs on every level, each level's interior holds
`f` exactly (independent of prolongation); interior+SYNC (not "everywhere")
also keeps validity tracking consistent so the poison check after
`MakeNewLevelFromCoarse` doesn't fire on not-yet-synced ghosts. The distinct
decade weights make any axis swap/transpose visible.

| Parfile (`TestPlanes/test/`) | Coverage |
|---|---|
| `planes-single-level.par` | uniform grid; 3 principal planes + 1 offset/axis; all 8 centerings; both writers |
| `planes-amr.par` | 3-level AMR (refined cube); every plane crosses every level (per-level snap, mrgtree shadowing, openPMD AMR attrs) |
| `planes-edge-cases.par` | half-integer / negative / over-precise / out-of-domain elevations |

The grid is small with a **distinct cell count and non-zero origin per axis**
(`16×20×24`, `x∈[-4,12]`, `y∈[-2,18]`, `z∈[-6,18]`, `dx=1`): distinct sizes
turn an in-plane transpose into a shape/value mismatch and the offsets catch
`ProbLo==0` bugs. `max_grid_size=8`/`blocking_factor=2` force several boxes,
and the single-level + AMR cases run on **2 ranks**, exercising the Silo MPI
gather / multi-file metafile and openPMD collective writes. The AMR cube uses
`ddf` prolongation with the poison guard left on.

`verify_planes.py` checks, for every emitted value, that it equals `f` at the
point's reconstructed coordinate — in-plane coords from the file's own mesh
metadata, the normal coord from an independent replay of `snap_to_grid_index`
(so neither can mask the other) — plus **coarse-level coverage** (the interior
points must tile the whole domain plane → catches silently dropped data) and
that out-of-domain planes produced **no** file. Reading is isolated in the
reusable, geometry-free `TestPlanes/test/plane_readers/` package (one kernel
per format, returning a uniform `Slab` list); `verify_planes.py` owns only the
parfile/geometry and comparison logic.

Both writers are numerically verified. openPMD via `openpmd_api` (pip); Silo
`.silo` (`DB_HDF5`) via the LLNL `Silo` module (apt `python3-silo`, not on
PyPI): its `DBfile` has only `GetToc`/`GetVar`/`GetVarInfo`, so `GetVarInfo`
gives each object's HDF5 component paths (`value0`, `coord0`/`coord1`),
centering and interior range, and `h5py` reads the arrays (reshape `(nb,na)` —
the buffer is axis_a-fastest, so trust the stride, not h5py's 2D shape). Only
interior cells are checked (AMR coarse-fine ghosts are prolongation-filled, not
exact `f`). Without the Silo module it degrades to a smoke check;
`--silo-mode inspect` dumps the layout/API.

A **golden-reference regression** check sits on top: committed binary reference
output under `TestPlanes/test/golden/<parfile>/`, generated by
`scripts/make-golden-planes.sh` (runs only the executable — no Python — so it
works on Frontier). Present golden is read back through the same reader and
compared data-level (point→value, keyed on centering, so it is independent of
MPI decomposition and of Silo's per-box variable names); absent golden is a
soft skip. Only the small openPMD `.bp5` golden is committed; the (large)
Silo golden is a stub — when committed it compares the same way.

CI runs `scripts/test-planes.sh` after the build (`.github/workflows/ci.yml`;
also via `agent_scripts/test.sh` locally): it pip-installs `openpmd-api`/`h5py`,
apt-installs `python3-silo`, runs the parfiles, and verifies both writers. To
regenerate golden after an intentional change, run `make-golden-planes.sh` and
commit `TestPlanes/test/golden/`.

## Why Cartesian-only

`CoordinatesX::vertex_coords` always stores world-space Cartesian (x,y,z), even
on spherical patches — so `z=12` is literal world-z. A constant-z surface
through a curvilinear patch is curved in index space, deferred to a future
interpolation path.
