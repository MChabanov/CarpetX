# 2D Plane Output for Silo and openPMD

Axis-aligned 2D slabs (xy / xz / yz) at chosen elevations, written
through the existing Silo and openPMD code paths. File-per-plane,
AMR-aware (per-level snap), Cartesian patches only.

## Status

- **Commit 1**: shared infrastructure — parameters, spec parsing,
  tag formatting, per-level grid snap, 3D-FAB slab extraction.
- **Commit 2 (this commit)**: Silo writer for non-staggered (rank-0
  all-VC and rank-3 all-CC) groups; dispatch wired in
  `io.cxx::OutputGH`; smoke-test parfile under `TestOutput`. Skips
  rank-1/2 staggered groups and non-Cartesian patches with one-time
  warnings.
- **Commit 3**: Silo per-group-mesh fallback for rank-1 / rank-2
  staggered groups (avoids the `assert(0)` at `io_silo.cxx:1118`).
- **Commit 4**: openPMD writer, mirroring the Silo design.

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
an extra 2-char in-plane centering suffix (`vv`/`cc`/`cv`/`vc`) when
the per-group fallback path is used (commit 3).

## Snapping

Elevations are rounded to `frac_precision` digits at parse (warning
fires once per unique spec on precision overflow). The rounded value
is both the tag content and the snap target. Per level, snap to the
nearest vertex (VC: `x0 + i·dx`) or cell center (CC: `x0 + (i+½)·dx`).
Out-of-domain elevations and non-Cartesian patches produce a one-time
warning and are skipped.

## Output files

Per plane: one data file per ioproc plus a metafile, all under a
per-iteration subdirectory:

```
<out_dir>/<sim>.<plane_tag>.it<8-digit-iter>.silo_planes.dir/
    <sim>.<plane_tag>.it<...>.p<6-digit-ioproc>.silo
<out_dir>/<sim>.<plane_tag>.it<...>.silo                metafile
<out_dir>/<sim>.<plane_tag>.silo_planes.visit          VisIt index
```

The metafile holds `DBPutMultimesh`/`DBPutMultivar` references to all
per-block 2D quadmeshes/quadvars. (Full AMR `DBmrgtree` deferred —
data still displays in VisIt, but level shadowing is manual.)

## File layout

```
CarpetX/src/io_planes.{hxx,cxx}        spec parsing, tag formatting,
                                       per-level snap, slab extraction
CarpetX/src/io_silo_planes.{hxx,cxx}   OutputSiloPlanes — rank-0/3 path
CarpetX/src/io.cxx                     dispatch (out_silo_planes block)
CarpetX/param.ccl                      8 new parameters
CarpetX/src/make.code.defn             two new sources registered
TestOutput/test/output-silo-planes.par smoke-test parfile
```

## Why Cartesian-only

`CoordinatesX::vertex_coords` always stores world-space Cartesian
(x, y, z) — even on spherical patches. So `z=12` is literal world-z,
not "third logical coordinate = 12". Constant-z through a curvilinear
patch is a curved sheet in index space, deferred to a future
interpolation path.
