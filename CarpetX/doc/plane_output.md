# 2D Plane Output for Silo and openPMD

Axis-aligned 2D slabs (xy / xz / yz) at chosen elevations, written
through the existing Silo and openPMD code paths. File-per-plane,
AMR-aware (per-level snap), Cartesian patches only.

## Status

- **Commit 1**: shared infrastructure — parameters, spec parsing,
  tag formatting, per-level grid snap, 3D-FAB slab extraction.
- **Commit 2**: Silo writer for non-staggered (rank-0 all-VC and
  rank-3 all-CC) groups; dispatch wired in `io.cxx::OutputGH`;
  smoke-test parfile under `TestOutput`.
- **Commit 3**: Silo per-group-mesh path for rank-1 / rank-2
  staggered groups. All 3D indextypes now handled; mesh sharing key
  extended with an in-plane centering tag so per-group meshes don't
  collide with the shared mesh. Avoids the `assert(0)` at
  `io_silo.cxx:1118`.
- **Commit 4**: openPMD writer (`OutputOpenPMDPlanes`). Per-plane
  local `openPMD::Series` (file-per-plane via the `it%08T` template),
  per-(group, patch, level) 2D mesh, per-component `storeChunk` with
  shared_ptr aliasing. Cell-centring encoded via `setPosition` per
  axis (no per-group-mesh dichotomy needed). Interior-only data,
  matching the existing 3D openPMD convention.
- **Commit 5**: Silo plane writer now emits a full AMR `DBmrgtree`
  alongside the multimesh (levelmaps + childmaps, `lvlRatios` /
  `ijkExts` / `xyzExts` / `rank` mrgvars, plus `DBOPT_EXTENTS` /
  `DBOPT_ZONECOUNTS` / `DBOPT_MRGTREE_NAME` on the multimesh). VisIt
  now shadows coarse-level data with finer levels where they
  overlap, instead of overlaying all levels. Slab enumeration is
  level-major to match the 3D writer's ordering convention; child
  relationships use 2D in-plane box refine-by-2 and overlap test.
- **Commit 6 (this commit)**: openPMD plane writer now writes
  per-iteration AMR-hierarchy attributes mirroring the 3D openPMD
  writer: `numDims` / `numPatches` / `patchSuffixes`, per-patch
  `numLevels` / `levelSuffixes`, per-(patch, level) `chunkInfo`
  (2D bounds of intersecting FABs, reversed for openPMD C order) and
  `iteration_num` / `iteration_den`. Plus plane-specific attributes
  `planeTag`, `planeNormalAxis`, `planeElevation` so readers can
  identify each plane slab. Intersection is computed in world
  coordinates (indextype-agnostic).

## Parameters

| Parameter | Default | Purpose |
|---|---|---|
| `out_silo_planes` | `""` | Specs, e.g. `"xy:0.0, xy:10.0, xz:-5.0"` |
| `out_silo_plane_vars` | `""` | Variable regex |
| `out_silo_planes_every` | `-1` | Cadence (-1 ⇒ `out_silo_every`) |
| `out_openpmd_planes` / `_plane_vars` / `_planes_every` | — | Same, for openPMD |
| `out_planes_int_precision` | `4` | Min integer digits in tag |
| `out_planes_frac_precision` | `3` | Fractional digits in tag |

Spec syntax: `<axes>:<elevation>`, comma-separated. `<elevation>` is
world-coordinate along the normal axis (z, y, or x).

## Tag format

`<plane>_<axis>_<sign><int>p<frac>` — e.g. `xy_z_pos0012p500` for
z=+12.500, `xz_y_neg0003p000` for y=−3.000. Survives Silo's
`legalize_name` and HDF5 / ADIOS2 identifier rules. Mesh names get
an extra 2-char in-plane centering suffix (`cv` / `vc`) when the
per-group mesh path is used for in-plane rank-1 staggered groups;
rank-0 / rank-2 groups omit the suffix and share a mesh.

## Snapping

Elevations are rounded to `frac_precision` digits at parse (warning
fires once per unique spec on precision overflow). The rounded value
is both the tag content and the snap target. Per level, snap to the
nearest vertex (VC: `x0 + i·dx`) or cell center (CC: `x0 + (i+½)·dx`).
Out-of-domain elevations and non-Cartesian patches produce a one-time
warning and are skipped.

## Output files

Per plane: one data file per ioproc plus a metafile (Silo) or one
file per iteration (openPMD), under per-iteration subdirectories:

```
Silo:
  <out_dir>/<sim>.<plane_tag>.it<...>.silo_planes.dir/
      <sim>.<plane_tag>.it<...>.p<6-digit-ioproc>.silo
  <out_dir>/<sim>.<plane_tag>.it<...>.silo               metafile
  <out_dir>/<sim>.<plane_tag>.silo_planes.visit          VisIt index

openPMD:
  <out_dir>/<sim>.<plane_tag>.it<...>.<ext>              one per iteration
  <out_dir>/<sim>.<plane_tag>.openpmd.visit              VisIt index
```

Silo metafile holds `DBPutMultimesh` / `DBPutMultivar` references to
all per-block 2D quadmeshes/quadvars, plus a `DBmrgtree` with
`amr_decomp/levels/patches` regions and `lvlRatios` / `ijkExts` /
`xyzExts` / `rank` mrgvars. VisIt uses the mrgtree to shadow
coarse-level data with finer levels where they overlap.

openPMD output uses MPI-collective `storeChunk`; cell-centring is
encoded per-component via `setPosition({0.5, 0.5})` (etc.) and the
mesh extent is the level's vertex grid along the two in-plane axes.
Data is interior-only (ghosts stripped), matching the existing 3D
openPMD output convention; Silo plane data includes ghost cells with
`DBOPT_LO_OFFSET` / `HI_OFFSET` markers, matching its 3D counterpart.

Each openPMD iteration also carries AMR-hierarchy attributes
(`numPatches`, `patchSuffixes`, per-patch `numLevels<suffix>` and
`levelSuffixes<suffix>`, per-(patch, level) `chunkInfo<suffix>` /
`iteration_num<suffix>` / `iteration_den<suffix>`), plus
plane-specific keys `planeTag`, `planeNormalAxis`, `planeElevation`.
Downstream tools can reconstruct which FABs contributed to which
slab without re-running the snap logic.

## File layout

```
CarpetX/src/io_planes.{hxx,cxx}        spec parsing, tag formatting,
                                       per-level snap, slab extraction
CarpetX/src/io_silo_planes.{hxx,cxx}   OutputSiloPlanes; shared mesh
                                       for in-plane rank-0/2, per-group
                                       mesh for rank-1 staggered groups
CarpetX/src/io_openpmd.{hxx,cxx}       OutputOpenPMDPlanes; per-plane
                                       Series, per-(group,patch,level)
                                       2D Mesh, setPosition for centring
CarpetX/src/io.cxx                     dispatch (silo & openpmd planes)
CarpetX/param.ccl                      8 new parameters
CarpetX/src/make.code.defn             io_planes.cxx, io_silo_planes.cxx
TestOutput/test/output-silo-planes.par     smoke-test (Silo)
TestOutput/test/output-openpmd-planes.par  smoke-test (openPMD)
```

## Testing

Thorn `TestPlanes` provides a rigorous, CI-integrated test suite for the plane
writers. It declares one grid function per index-type centering (the eight V/C
combinations: vertex `gf000`, cell `gf111`, the three face- and three
edge-centred variants) and fills each with the same analytic, axis-distinct
field

```
f(x, y, z) = x + 100 y + 10000 z
```

evaluated at the grid function's own centering-dependent world coordinate. The
interior is set directly (writes interior, then SYNC) at `initial` and
`postregrid`, and `CCTK_INITIAL` re-runs on every level, so each level's
interior holds `f` exactly -- independent of prolongation -- which is the data
the openPMD writer emits and the numeric check uses. Writing interior + SYNC
(not "everywhere") also keeps validity tracking consistent with a freshly
regridded level, so the poison check after `MakeNewLevelFromCoarse` does not
fire on not-yet-synced ghosts. The distinct decade weights make an axis swap or
transpose immediately visible.

Parfiles (`TestPlanes/test/`):

| Parfile | Coverage |
|---|---|
| `planes-single-level.par` | uniform grid; xy/xz/yz principal planes + one offset per direction; all 8 centerings; Silo + openPMD |
| `planes-amr.par` | 3-level AMR (refined cube about the centre); every plane crosses every level (per-level snap, Silo mrgtree shadowing, openPMD per-(patch,level) attributes) |
| `planes-edge-cases.par` | half-integer / negative / over-precise elevations, and an out-of-domain elevation that must be skipped |

The grid is deliberately small and uses a distinct cell count and a distinct,
non-zero origin per axis (`16×20×24` cells, `x∈[-4,12]`, `y∈[-2,18]`,
`z∈[-6,18]`, `dx=1`): the distinct sizes make a transposed in-plane axis a
shape/value mismatch, and the non-zero origins catch `ProbLo == 0` assumptions.
The AMR run refines a half-width-4 cube about the centre `(4,8,6)` over three
levels (so the cell-centred normal coordinate differs per level, e.g. `6.0` on
the coarse level vs `6.125` on the finest).

Verification (`TestPlanes/test/verify_planes.py`) reads the written plane files
back and asserts that every stored value equals `f` at the point's
reconstructed world coordinate. The two coordinate sources are deliberately
independent: in-plane coordinates come from each file's own mesh metadata
(catching wrong `gridSpacing`/offset/`position`/centering), while the
normal-axis coordinate is obtained by independently replaying
`snap_to_grid_index` from the parfile geometry (catching a wrong slice). A
compact per-slab summary is also diffed against committed golden files
(`TestPlanes/test/golden/`); a missing golden file is a soft skip (the analytic
check is the gate), while a present one is enforced.

openPMD is read via `openpmd_api` (pip-installable). Silo `.silo` files are
read with the LLNL `Silo` Python module when available (cleanest API; e.g.
conda-forge `silo-nompi`, and present in the ET container); otherwise the Silo
path falls back to a structural smoke check (valid HDF5 + datasets present,
since the writer uses the `DB_HDF5` driver). `verify_planes.py --silo-mode
inspect` dumps the on-disk Silo layout and module API to ease finalizing the
reader in a given container.

Beyond the per-point value check, the verifier also checks **coarse-level
coverage**: the union of interior in-plane points on level 0 must tile the
whole domain plane. This catches data that is silently *missing* rather than
wrong -- the typical multi-rank failure mode (a rank's slab dropped, or a
gather mismatch). The grids are forced into several boxes (`max_grid_size=8`,
`blocking_factor=2`), and the single-level and AMR cases run on **two ranks**,
so the Silo `MPI_Send`/`MPI_Recv` gather, the multi-file metafile, and the
openPMD per-rank collective `storeChunk` are all exercised.

CI runs `scripts/test-planes.sh` after the build (see
`.github/workflows/ci.yml`), which executes the parfiles (single-level and AMR
on two ranks) and runs the verifier.
The same script runs in the local container loop via `agent_scripts/test.sh`.
To refresh golden files after an intentional change, run the script with
`UPDATE_GOLDEN=1`.

## Why Cartesian-only

`CoordinatesX::vertex_coords` always stores world-space Cartesian
(x, y, z) — even on spherical patches. So `z=12` is literal world-z,
not "third logical coordinate = 12". Constant-z through a curvilinear
patch is a curved sheet in index space, deferred to a future
interpolation path.
