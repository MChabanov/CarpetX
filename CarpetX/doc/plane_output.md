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

## Why Cartesian-only

`CoordinatesX::vertex_coords` always stores world-space Cartesian
(x, y, z) — even on spherical patches. So `z=12` is literal world-z,
not "third logical coordinate = 12". Constant-z through a curvilinear
patch is a curved sheet in index space, deferred to a future
interpolation path.
