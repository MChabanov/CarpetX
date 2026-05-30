#!/usr/bin/env python3
"""Verify CarpetX 2D plane output (Silo and openPMD) for thorn TestPlanes.

Every grid function written by TestPlanes holds the analytic field

    f(x, y, z) = x + 100*y + 10000*z

evaluated at its own centering-dependent world coordinate. This script reads
the plane files back and checks, for every value written, that it equals f at
the point's reconstructed world coordinate -- a reference-free, analytic check
of the snap, slab extraction, axis ordering and per-centering coordinates.

In-plane coordinates are taken from the file's own mesh metadata (so a wrong
gridSpacing/offset/position is caught); the normal-axis coordinate is obtained
by independently replaying CarpetX's snap_to_grid_index from the parfile
geometry (so a wrong slice is caught). The two are deliberately sourced
differently so neither can mask a bug in the other.

If --golden-dir is given and holds committed reference output for this parfile
(<golden-dir>/<sim>/, the binary .bp5/.silo files produced by
scripts/make-golden-planes.sh), the golden files are read back through the same
readers and the data is compared against the fresh output -- a regression
snapshot on top of the analytic check. The comparison is data-level (point
coordinate -> value), so it is robust to MPI-decomposition differences.

Usage:
    verify_planes.py --parfile P.par --out-dir DIR [--golden-dir D]
                     [--silo-mode {auto,verify,smoke,skip,inspect}]
                     [--require-openpmd] [--require-silo]

Exit status is non-zero on any verification failure.
"""

import argparse
import glob
import math
import os
import re
import sys

# f(x,y,z) = X_WEIGHTS . (x,y,z)
X_WEIGHTS = (1.0, 100.0, 10000.0)
AXIS_NAMES = ("x", "y", "z")

# Default comparison tolerances (real64). expected magnitudes reach ~3e5.
DEFAULT_RTOL = 1.0e-9
DEFAULT_ATOL = 1.0e-6


def f_analytic(x, y, z):
    return X_WEIGHTS[0] * x + X_WEIGHTS[1] * y + X_WEIGHTS[2] * z


def lround_half_away(x):
    """Match C++ std::lround: round half away from zero."""
    return math.floor(x + 0.5) if x >= 0 else math.ceil(x - 0.5)


# ---------------------------------------------------------------------------
# Parfile parsing (geometry + plane specs + tag precision)
# ---------------------------------------------------------------------------

class Geometry:
    def __init__(self, par):
        self.ncells = [int(par["carpetx::ncells_%s" % a]) for a in AXIS_NAMES]
        self.xmin = [float(par["carpetx::%smin" % a]) for a in AXIS_NAMES]
        self.xmax = [float(par["carpetx::%smax" % a]) for a in AXIS_NAMES]
        self.max_num_levels = int(par.get("carpetx::max_num_levels", "1"))
        self.int_precision = int(par.get("carpetx::out_planes_int_precision", "4"))
        self.frac_precision = int(par.get("carpetx::out_planes_frac_precision", "3"))

    def dx0(self, axis):
        return (self.xmax[axis] - self.xmin[axis]) / self.ncells[axis]

    def dx(self, axis, level):
        return self.dx0(axis) / (2 ** level)

    def ncells_at(self, axis, level):
        return self.ncells[axis] * (2 ** level)

    def in_domain(self, axis, elevation):
        return self.xmin[axis] <= elevation <= self.xmax[axis]

    def snap_coord(self, axis, level, cell_centered, elevation):
        """Replay CarpetX::snap_to_grid_index and return the world coordinate
        of the snapped slice, or None if the elevation is out of range."""
        x0 = self.xmin[axis]
        dx = self.dx(axis, level)
        ncells = self.ncells_at(axis, level)
        r = (elevation - x0) / dx
        if cell_centered:
            i_snap = lround_half_away(r - 0.5)
            i_lo, i_hi = 0, ncells - 1
            if i_snap < i_lo or i_snap > i_hi:
                return None
            return x0 + (i_snap + 0.5) * dx
        else:
            i_snap = lround_half_away(r)
            i_lo, i_hi = 0, ncells
            if i_snap < i_lo or i_snap > i_hi:
                return None
            return x0 + i_snap * dx


def round_elevation(elev, frac_precision):
    """Replay parse_planes elevation rounding."""
    scale = 10 ** frac_precision
    sign = -1.0 if elev < 0 else 1.0
    total = lround_half_away(abs(elev) * scale)
    return sign * (total / scale)


def parse_parfile(path):
    """Return a dict of lowercased 'thorn::param' / '$var' -> string value.

    Handles double-quoted (possibly multi-line) and bare values, strips full
    and trailing '#' comments, and resolves '$var' references.
    """
    with open(path) as fh:
        text = fh.read()

    # Strip '#' comments that are not inside double quotes, line by line.
    lines = []
    in_str = False
    for raw in text.splitlines():
        buf = []
        for ch in raw:
            if ch == '"':
                in_str = not in_str
                buf.append(ch)
            elif ch == "#" and not in_str:
                break
            else:
                buf.append(ch)
        lines.append("".join(buf))
    clean = "\n".join(lines)

    params = {}
    # NAME = "value"  (DOTALL for multi-line quoted strings)
    pat = re.compile(r'([\w:$]+)\s*=\s*("(?:[^"]*)"|\S+)', re.DOTALL)
    for m in pat.finditer(clean):
        name = m.group(1).strip().lower()
        val = m.group(2).strip()
        if val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
        params[name] = val

    # Resolve "$var" references.
    for name, val in list(params.items()):
        v = val.strip()
        if v.startswith("$") and v.lower() in params:
            params[name] = params[v.lower()]
    return params


def parse_plane_specs(spec_string):
    """Replay parse_planes: return list of (normal_axis, raw_elevation)."""
    specs = []
    for item in spec_string.split(","):
        item = item.strip()
        if not item or ":" not in item:
            continue
        axes, elev = item.split(":", 1)
        axes = axes.strip()
        normal = {"yz": 0, "xz": 1, "xy": 2}.get(axes)
        if normal is None:
            continue
        try:
            specs.append((normal, float(elev.strip())))
        except ValueError:
            continue
    return specs


# ---------------------------------------------------------------------------
# Slab record produced by the readers
# ---------------------------------------------------------------------------

class Slab:
    """One 2D slab of a single variable at one (patch, level)."""
    def __init__(self, fmt, var, centering, normal_axis, level, patch,
                 axis_a, axis_b, coords_a, coords_b, values):
        self.fmt = fmt                  # "openpmd" or "silo"
        self.var = var                  # e.g. "gf000"
        self.centering = centering      # (cx, cy, cz), 1 == cell-centred
        self.normal_axis = normal_axis  # 0/1/2
        self.level = level
        self.patch = patch
        self.axis_a = axis_a
        self.axis_b = axis_b
        self.coords_a = coords_a        # 1D world coords along axis_a (len na)
        self.coords_b = coords_b        # 1D world coords along axis_b (len nb)
        self.values = values            # 2D [b][a]


def centering_from_varname(var):
    m = re.search(r"gf([012])([012])([012])", var.lower())
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2)), int(m.group(3)))


# ---------------------------------------------------------------------------
# openPMD reader
# ---------------------------------------------------------------------------

def read_openpmd(out_dir, sim, geom):
    """Yield (tag, normal_axis, elevation, [Slab, ...]) per plane file group.

    Returns (records, ok, msg). records is a list of tuples; ok is False only
    if openPMD support is unavailable.
    """
    try:
        import openpmd_api as io
    except Exception as exc:  # noqa: BLE001
        return [], False, "openpmd_api unavailable: %s" % exc

    import numpy as np

    # Discover plane file groups: <sim>.<tag>.it<NNNNNNNN>.<ext>
    exts = ("bp5", "bp", "bp4", "h5")
    files = []
    for ext in exts:
        files += glob.glob(os.path.join(out_dir, "%s.*.it*.%s" % (sim, ext)))
    groups = {}
    for path in files:
        base = os.path.basename(path)
        m = re.match(r"^%s\.(.+)\.it(\d+)\.(\w+)$" % re.escape(sim), base)
        if not m:
            continue
        tag, _it, ext = m.group(1), m.group(2), m.group(3)
        groups.setdefault((tag, ext), True)

    records = []
    for (tag, ext) in sorted(groups):
        # The writer emits no file for an out-of-domain elevation; that is
        # asserted separately by find_out_of_domain_files. Defensively, if such
        # a file does exist, do not open it here -- its empty chunkInfo
        # attributes would make openpmd-api throw at Series open -- so that the
        # clean "writer produced a file" error is reported instead of a crash.
        parsed = parse_plane_tag(tag)
        if parsed is not None and not geom.in_domain(parsed[0], parsed[1]):
            records.append((tag, parsed[0], parsed[1], []))
            continue
        pattern = os.path.join(out_dir, "%s.%s.it%%08T.%s" % (sim, tag, ext))
        series = io.Series(pattern, io.Access.read_only)
        for it_index in series.iterations:
            it = series.iterations[it_index]
            normal_axis = int(it.get_attribute("planeNormalAxis"))
            elevation = float(it.get_attribute("planeElevation"))

            slabs = []
            for mesh_name in it.meshes:
                mesh = it.meshes[mesh_name]
                lev_m = re.search(r"_lev(\d+)", mesh_name)
                patch_m = re.search(r"_patch(\d+)", mesh_name)
                level = int(lev_m.group(1)) if lev_m else 0
                patch = int(patch_m.group(1)) if patch_m else 0

                ggo = list(mesh.grid_global_offset)   # [off_b, off_a]
                gsp = list(mesh.grid_spacing)          # [sp_b, sp_a]

                axis_a = 1 if normal_axis == 0 else 0
                axis_b = 1 if normal_axis == 2 else 2

                for comp_name in mesh:
                    var = comp_name
                    centering = centering_from_varname(var)
                    if centering is None:
                        continue
                    rc = mesh[comp_name]
                    pos = list(rc.position)            # [pos_b, pos_a]

                    chunks = rc.available_chunks()
                    for ch in chunks:
                        off = list(ch.offset)          # [o_b, o_a]
                        ext_ = list(ch.extent)         # [c_b, c_a]
                        data = rc.load_chunk(off, ext_)
                        series.flush()
                        arr = np.asarray(data)         # shape [c_b, c_a]

                        nb, na = ext_[0], ext_[1]
                        coords_b = [ggo[0] + (off[0] + j + pos[0]) * gsp[0]
                                    for j in range(nb)]
                        coords_a = [ggo[1] + (off[1] + i + pos[1]) * gsp[1]
                                    for i in range(na)]
                        slabs.append(Slab(
                            "openpmd", var, centering, normal_axis, level,
                            patch, axis_a, axis_b, coords_a, coords_b,
                            arr.reshape(nb, na)))
            records.append((tag, normal_axis, elevation, slabs))
        del series
    return records, True, "ok"


# ---------------------------------------------------------------------------
# Silo reader (best effort: full verify if a reader is available, else smoke)
# ---------------------------------------------------------------------------

def silo_meta_files(out_dir, sim):
    """Per-plane Silo metafiles (no .p<proc> suffix)."""
    return sorted(p for p in glob.glob(os.path.join(out_dir, "%s.*.it*.silo" % sim))
                  if ".silo_planes.dir" not in p)


def silo_data_files(out_dir, sim):
    """Per-process Silo data files inside the .silo_planes.dir subdirectories."""
    return sorted(glob.glob(os.path.join(
        out_dir, "%s.*.it*.silo_planes.dir" % sim, "*.silo")))


# Inverse of CarpetX::format_plane_tag, e.g. "xy_z_pos0006p000" -> (2, 6.0).
def parse_plane_tag(tag):
    parts = tag.split("_")
    if len(parts) < 3:
        return None
    normal = {"yz": 0, "xz": 1, "xy": 2}.get(parts[0])
    if normal is None:
        return None
    mag = parts[2]
    sign = -1.0 if mag.startswith("neg") else 1.0
    mag = mag[3:]
    if "p" in mag:
        ip, fp = mag.split("p", 1)
        elev = float(ip) + (float(fp) / (10 ** len(fp)) if fp else 0.0)
    else:
        elev = float(mag)
    return normal, sign * elev


# Silo centering constants (silo.h)
DB_NODECENT = 110
DB_ZONECENT = 111


def read_silo(out_dir, sim, geom):
    """Return (records, mode, msg). mode is 'verify' (numeric) or 'smoke'.

    Full numeric verification uses the LLNL Silo Python module if it is
    importable (the cleanest, most reliable reader; available e.g. from
    conda-forge `silo-nompi` and expected in the einsteintoolkit/carpetx
    image). Reading the DB_HDF5-backed files directly with h5py is possible but
    the internal layout is version-dependent, so absent the Silo module we fall
    back to a structural smoke check (valid HDF5 + datasets present). Use
    --silo-mode inspect to dump the real container layout/API.
    """
    try:
        import Silo  # noqa: F401
    except Exception:  # noqa: BLE001
        return _silo_smoke(out_dir, sim)
    try:
        records = _read_silo_module(out_dir, sim, geom)
    except Exception as exc:  # noqa: BLE001
        recs, mode, msg = _silo_smoke(out_dir, sim)
        return recs, mode, "Silo module present but read failed (%s); %s" % (
            exc, msg)
    if not records:
        recs, mode, msg = _silo_smoke(out_dir, sim)
        return recs, mode, "Silo module present but produced no slabs; " + msg
    return records, "verify", "read via Silo Python module"


def _silo_attr(obj, *names):
    for n in names:
        if hasattr(obj, n):
            return getattr(obj, n)
    raise AttributeError("none of %s on %r" % (names, type(obj)))


def _read_silo_module(out_dir, sim, geom):
    """Read per-process Silo plane data files via the LLNL `Silo` module.

    Defensive about exact SWIG attribute names; any mismatch raises and the
    caller falls back to the smoke check. Confirm the real API with
    --silo-mode inspect.
    """
    import Silo
    import numpy as np

    records_by_tag = {}
    for path in silo_data_files(out_dir, sim):
        base = os.path.basename(os.path.dirname(path))  # ...<tag>.it<N>.silo_planes.dir
        m = re.match(r"^%s\.(.+)\.it(\d+)\.silo_planes\.dir$" % re.escape(sim), base)
        if not m:
            continue
        tag = m.group(1)
        parsed = parse_plane_tag(tag)
        if parsed is None:
            continue
        normal_axis, elevation = parsed
        axis_a, axis_b = ([1, 2], [0, 2], [0, 1])[normal_axis]

        db = Silo.Open(path, Silo.DB_READ)
        try:
            toc = db.GetToc()
            qvar_names = list(_silo_attr(toc, "qvar_names"))
            slabs = []
            for qvname in qvar_names:
                qv = db.GetQuadvar(qvname)
                centering = centering_from_varname(qvname)
                if centering is None:
                    continue
                vals = np.asarray(_silo_attr(qv, "vals")[0], dtype=float)
                dims = list(_silo_attr(qv, "dims"))           # [na, nb]
                cent = int(_silo_attr(qv, "centering"))
                meshname = _silo_attr(qv, "meshname", "meshid")
                if isinstance(meshname, bytes):
                    meshname = meshname.decode()
                lev_m = re.search(r"\.rl(\d+)", qvname)
                patch_m = re.search(r"\.m(\d+)", qvname)
                level = int(lev_m.group(1)) if lev_m else 0
                patch = int(patch_m.group(1)) if patch_m else 0

                qm = db.GetQuadmesh(meshname)
                mcoords = _silo_attr(qm, "coords")
                coord_a = np.asarray(mcoords[0], dtype=float)
                coord_b = np.asarray(mcoords[1], dtype=float)

                na, nb = int(dims[0]), int(dims[1])
                values = vals.reshape(nb, na)               # row-major, a fastest
                if cent == DB_ZONECENT:
                    coords_a = [(coord_a[i] + coord_a[i + 1]) / 2 for i in range(na)]
                    coords_b = [(coord_b[j] + coord_b[j + 1]) / 2 for j in range(nb)]
                else:
                    coords_a = list(coord_a[:na])
                    coords_b = list(coord_b[:nb])
                slabs.append(Slab("silo", qvname, centering, normal_axis,
                                  level, patch, axis_a, axis_b,
                                  coords_a, coords_b, values))
        finally:
            db.Close()

        rec = records_by_tag.setdefault(tag, [tag, normal_axis, elevation, []])
        rec[3].extend(slabs)
    return [tuple(r) for r in records_by_tag.values()]


def _silo_smoke(out_dir, sim):
    files = silo_meta_files(out_dir, sim) + silo_data_files(out_dir, sim)
    if not files:
        return [], "smoke", "no Silo plane files found"
    try:
        import h5py
    except Exception:  # noqa: BLE001
        for path in files:
            if os.path.getsize(path) == 0:
                raise SystemExit("Silo file is empty: %s" % path)
        return [], "smoke", "%d Silo files present (size>0; no h5py)" % len(files)

    n_datasets = 0
    for path in files:
        try:
            with h5py.File(path, "r") as h5:
                def _count(name, obj):
                    nonlocal n_datasets
                    if isinstance(obj, h5py.Dataset):
                        n_datasets += 1
                h5.visititems(_count)
        except Exception as exc:  # noqa: BLE001
            raise SystemExit("Silo file is not valid HDF5: %s (%s)" % (path, exc))
    if n_datasets == 0:
        raise SystemExit("Silo files contain no datasets (writer produced "
                         "nothing?)")
    return [], "smoke", "%d Silo files valid HDF5, %d datasets" % (
        len(files), n_datasets)


def _inspect_one_hdf5(path, label, max_show_values=32):
    import h5py
    import numpy as np
    print("--- %s: %s ---" % (label, path))
    with h5py.File(path, "r") as h5:
        # Root attributes
        if len(h5.attrs):
            print("  / attrs: %s" % {k: _short_attr(v) for k, v in h5.attrs.items()})

        def show(name, obj):
            if isinstance(obj, h5py.Group):
                a = (" attrs=%s" % {k: _short_attr(v) for k, v in obj.attrs.items()}
                     if len(obj.attrs) else "")
                print("  [G] %s%s" % (name, a))
                return
            # Dataset
            line = "  [D] %s shape=%s dtype=%s" % (name, obj.shape, obj.dtype)
            if len(obj.attrs):
                line += " attrs=%s" % {k: _short_attr(v) for k, v in obj.attrs.items()}
            print(line)
            try:
                n = int(np.prod(obj.shape)) if obj.shape else 1
                if n and n <= max_show_values:
                    print("        values=%r" % (np.asarray(obj[()]).ravel().tolist()))
            except Exception as exc:  # noqa: BLE001
                print("        (could not read values: %s)" % exc)

        h5.visititems(show)


def _short_attr(v):
    try:
        import numpy as np
        if isinstance(v, bytes):
            return v.decode(errors="replace")
        arr = np.asarray(v)
        if arr.size > 12:
            return "<%s array, size %d>" % (arr.dtype, arr.size)
        return arr.tolist()
    except Exception:  # noqa: BLE001
        return repr(v)


def inspect_silo(out_dir, sim):
    """Dump the real Silo layout/API to help finalize the h5py reader."""
    print("=== Silo inspection ===")
    metas = silo_meta_files(out_dir, sim)
    datas = silo_data_files(out_dir, sim)
    print("metafiles: %d, data files: %d" % (len(metas), len(datas)))

    try:
        import h5py  # noqa: F401
        if datas:
            _inspect_one_hdf5(datas[0], "per-proc DATA file")
        if metas:
            _inspect_one_hdf5(metas[0], "METAFILE")
        if not datas and not metas:
            print("no Silo files found in %s" % out_dir)
    except Exception as exc:  # noqa: BLE001
        print("h5py inspection failed: %s" % exc)

    sample = (datas or metas)
    if sample:
        try:
            import Silo
            import numpy as np
            print("--- Silo module: API exploration ---")
            db = Silo.Open(sample[0], Silo.DB_READ)
            try:
                print("  dir(db): %s"
                      % [a for a in dir(db) if not a.startswith("_")])
            except Exception as e:  # noqa: BLE001
                print("  dir(db) failed: %s" % e)
            # NB: dir(toc) raises on this build; access fields by name only.
            toc = db.GetToc()
            qvn = list(getattr(toc, "qvar_names", []) or [])
            qmn = list(getattr(toc, "qmesh_names", []) or [])
            vn = list(getattr(toc, "var_names", []) or [])
            print("  #qvar=%d #qmesh=%d #simplevar=%d" % (len(qvn), len(qmn), len(vn)))
            print("  qvar[0]  = %r" % (qvn[0] if qvn else None))
            print("  qmesh[0] = %r" % (qmn[0] if qmn else None))
            print("  simplevars = %r" % (vn[:8],))

            def _show(label, val):
                if isinstance(val, dict):
                    print("    %s -> dict keys=%s" % (label, list(val.keys())))
                    for k, v in val.items():
                        try:
                            a = np.asarray(v)
                            sv = (a.shape if a.size > 16 else a.ravel().tolist())
                        except Exception:  # noqa: BLE001
                            sv = repr(v)[:100]
                        print("        [%r] %s = %s" % (k, type(v).__name__, sv))
                else:
                    try:
                        a = np.asarray(val)
                        sv = (a.shape if a.size > 16 else a.ravel().tolist())
                    except Exception:  # noqa: BLE001
                        sv = repr(val)[:200]
                    print("    %s -> %s %s" % (label, type(val).__name__, sv))

            for name, lab in ((qvn[0] if qvn else None, "qvar"),
                              (qmn[0] if qmn else None, "qmesh")):
                if name is None:
                    continue
                for meth in ("GetVarInfo", "GetVar"):
                    if not hasattr(db, meth):
                        continue
                    try:
                        _show("db.%s(%s)" % (meth, lab),
                              getattr(db, meth)(name))
                    except Exception as e:  # noqa: BLE001
                        print("    db.%s(%s) FAILED: %s" % (meth, lab, e))
            db.Close()
        except Exception as exc:  # noqa: BLE001
            print("Silo module exploration error: %s" % exc)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def verify_slabs(records, geom, rtol, atol):
    """Check every value == f at its reconstructed coordinate, and that the
    coarse level is fully covered (catches data dropped by a rank on a
    multi-process run). Returns (n_checked, errors)."""
    n_checked = 0
    errors = []
    eps = 1e-9
    # coverage[(tag, var)] -> set of distinct interior (a,b) coords on level 0
    coverage = {}
    cover_meta = {}
    for (tag, normal_axis, elevation, slabs) in records:
        rounded = round_elevation(elevation, geom.frac_precision)
        for s in slabs:
            cc_normal = bool(s.centering[normal_axis])
            ncoord = geom.snap_coord(normal_axis, s.level, cc_normal, rounded)
            if ncoord is None:
                errors.append("%s %s L%d: produced a slab but the elevation "
                              "%.6g snaps out of range" %
                              (tag, s.var, s.level, rounded))
                continue
            nb = len(s.coords_b)
            na = len(s.coords_a)
            local_errs = 0
            for j in range(nb):
                for i in range(na):
                    coord = [0.0, 0.0, 0.0]
                    coord[s.axis_a] = s.coords_a[i]
                    coord[s.axis_b] = s.coords_b[j]
                    coord[normal_axis] = ncoord
                    expected = f_analytic(*coord)
                    actual = float(s.values[j][i])
                    if abs(actual - expected) > atol + rtol * abs(expected):
                        if local_errs < 3:
                            errors.append(
                                "%s %s L%d patch%d at (a=%.6g,b=%.6g,n=%.6g): "
                                "got %.10g expected %.10g" %
                                (tag, s.var, s.level, s.patch, s.coords_a[i],
                                 s.coords_b[j], ncoord, actual, expected))
                        local_errs += 1
                    n_checked += 1
                    # Coverage: record interior in-plane points on level 0.
                    if s.level == 0:
                        a, b = s.coords_a[i], s.coords_b[j]
                        if (geom.xmin[s.axis_a] - eps <= a <= geom.xmax[s.axis_a] + eps
                                and geom.xmin[s.axis_b] - eps <= b <= geom.xmax[s.axis_b] + eps):
                            key = (tag, s.var)
                            coverage.setdefault(key, set()).add(
                                (round(a, 6), round(b, 6)))
                            cover_meta[key] = (s.centering, s.axis_a, s.axis_b)

    # Coarse-level completeness: the union of interior points must tile the
    # whole domain plane. A missing rank's contribution shows up as a shortfall.
    for key, pts in coverage.items():
        tag, var = key
        centering, axis_a, axis_b = cover_meta[key]
        na_exp = geom.ncells[axis_a] + (0 if centering[axis_a] else 1)
        nb_exp = geom.ncells[axis_b] + (0 if centering[axis_b] else 1)
        expected = na_exp * nb_exp
        if len(pts) != expected:
            errors.append(
                "%s %s L0 coverage: %d distinct interior points, expected %d "
                "(missing/extra data -- e.g. a rank's slab dropped)" %
                (tag, var, len(pts), expected))
    return n_checked, errors


# ---------------------------------------------------------------------------
# Golden files
# ---------------------------------------------------------------------------

def find_out_of_domain_files(out_dir, sim, geom):
    """The writers must emit no file for a plane whose elevation is out of
    domain. Flag any plane output file (openPMD or Silo) whose tag decodes to
    an out-of-domain elevation."""
    errors = []
    seen = set()
    patterns = ["%s.*.it*.h5", "%s.*.it*.bp", "%s.*.it*.bp4", "%s.*.it*.bp5",
                "%s.*.it*.silo"]
    for pat in patterns:
        for path in glob.glob(os.path.join(out_dir, pat % sim)):
            base = os.path.basename(path)
            m = re.match(r"^%s\.(.+)\.it(\d+)\.(\w+)$" % re.escape(sim), base)
            if not m:
                continue
            tag = m.group(1)
            if tag in seen:
                continue
            seen.add(tag)
            parsed = parse_plane_tag(tag)
            if parsed is not None and not geom.in_domain(parsed[0], parsed[1]):
                errors.append(
                    "writer produced an output file for out-of-domain plane "
                    "'%s' (elevation %.6g on axis %d); expected no file"
                    % (tag, parsed[1], parsed[0]))
    return errors


# Golden reference is the committed binary plane output (.bp5 / .silo) under
# <golden-dir>/<sim>/, generated by scripts/make-golden-planes.sh. The golden
# branch reads it back through the same reader and compares data values against
# the freshly produced output -- a regression snapshot on top of the analytic
# check. The comparison is purely data-level (point coordinate -> value), so it
# is robust to differences in MPI decomposition or chunk layout between the
# machine that produced the golden files and CI.

def slabs_to_value_map(records):
    m = {}
    for (tag, normal_axis, elevation, slabs) in records:
        for s in slabs:
            na = len(s.coords_a)
            nb = len(s.coords_b)
            for j in range(nb):
                for i in range(na):
                    key = (tag, s.var, s.level,
                           round(s.coords_a[i], 6), round(s.coords_b[j], 6))
                    m[key] = float(s.values[j][i])
    return m


def compare_to_golden(fresh_records, golden_records, rtol, atol, label):
    """Compare freshly produced data against the golden reference data.
    Returns (n_compared, errors)."""
    fm = slabs_to_value_map(fresh_records)
    gm = slabs_to_value_map(golden_records)
    fk, gk = set(fm), set(gm)
    errors = []
    missing = gk - fk
    extra = fk - gk
    if missing:
        errors.append("%s golden: %d reference point(s) missing from current "
                      "output (e.g. %s)" % (label, len(missing),
                                            sorted(missing)[0]))
    if extra:
        errors.append("%s golden: %d current point(s) absent from reference "
                      "(e.g. %s)" % (label, len(extra), sorted(extra)[0]))
    nmismatch = 0
    for k in fk & gk:
        if abs(fm[k] - gm[k]) > atol + rtol * abs(gm[k]):
            if nmismatch < 3:
                errors.append("%s golden mismatch at %s: current %.10g vs "
                              "reference %.10g" % (label, k, fm[k], gm[k]))
            nmismatch += 1
    if nmismatch:
        errors.append("%s golden: %d value(s) differ from reference"
                      % (label, nmismatch))
    return len(fk & gk), errors


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parfile", required=True)
    ap.add_argument("--out-dir", required=True,
                    help="directory containing the plane output")
    ap.add_argument("--golden-dir", default=None,
                    help="directory holding per-parfile golden subdirs "
                         "(<golden-dir>/<sim>/ with .bp5/.silo reference output)")
    ap.add_argument("--silo-mode", default="auto",
                    choices=["auto", "verify", "smoke", "skip", "inspect"])
    ap.add_argument("--require-openpmd", action="store_true")
    ap.add_argument("--require-silo", action="store_true")
    ap.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    ap.add_argument("--atol", type=float, default=DEFAULT_ATOL)
    args = ap.parse_args(argv)

    par = parse_parfile(args.parfile)
    geom = Geometry(par)
    sim = os.path.splitext(os.path.basename(args.parfile))[0]

    all_errors = []
    golden_subdir = (os.path.join(args.golden_dir, sim)
                     if args.golden_dir is not None else None)

    # --- Out-of-domain planes must produce no file at all ------------------
    ood = find_out_of_domain_files(args.out_dir, sim, geom)
    if ood:
        all_errors += ood
    else:
        print("out-of-domain planes: no files produced (as expected)")

    # --- openPMD -----------------------------------------------------------
    op_records, op_ok, op_msg = read_openpmd(args.out_dir, sim, geom)
    if not op_ok:
        msg = "openPMD: %s" % op_msg
        if args.require_openpmd:
            all_errors.append(msg)
        else:
            print("SKIP %s" % msg)
    else:
        n, errs = verify_slabs(op_records, geom, args.rtol, args.atol)
        all_errors += errs
        print("openPMD: %d plane file groups, %d values checked, %d errors"
              % (len(op_records), n, len(errs)))
        if golden_subdir is not None:
            if os.path.isdir(golden_subdir):
                g_records, g_ok, g_msg = read_openpmd(golden_subdir, sim, geom)
                if not g_ok:
                    all_errors.append("openPMD golden: %s" % g_msg)
                else:
                    ncmp, gerrs = compare_to_golden(op_records, g_records,
                                                    args.rtol, args.atol,
                                                    "openPMD")
                    all_errors += gerrs
                    print("openPMD golden: %d points compared vs reference, "
                          "%d issue(s)" % (ncmp, len(gerrs)))
            else:
                print("openPMD golden: %s absent (not enforced yet) -- "
                      "generate with scripts/make-golden-planes.sh"
                      % golden_subdir)

    # --- Silo --------------------------------------------------------------
    if args.silo_mode == "skip":
        print("Silo: skipped by request")
    elif args.silo_mode == "inspect":
        inspect_silo(args.out_dir, sim)
    else:
        si_records, si_mode, si_msg = read_silo(args.out_dir, sim, geom)
        print("Silo: mode=%s (%s)" % (si_mode, si_msg))
        if args.silo_mode == "verify" and si_mode != "verify":
            all_errors.append("Silo: full verification requested but reader "
                              "unavailable (%s)" % si_msg)
        if si_mode == "verify":
            n, errs = verify_slabs(si_records, geom, args.rtol, args.atol)
            all_errors += errs
            print("Silo: %d values checked, %d errors" % (n, len(errs)))
            if golden_subdir is not None:
                if os.path.isdir(golden_subdir):
                    g_records, g_mode, g_msg = read_silo(golden_subdir, sim,
                                                         geom)
                    if g_mode == "verify":
                        ncmp, gerrs = compare_to_golden(si_records, g_records,
                                                        args.rtol, args.atol,
                                                        "Silo")
                        all_errors += gerrs
                        print("Silo golden: %d points compared vs reference, "
                              "%d issue(s)" % (ncmp, len(gerrs)))
                    else:
                        print("Silo golden: reference present but not readable "
                              "here (%s)" % g_msg)
                else:
                    print("Silo golden: %s absent (not enforced yet)"
                          % golden_subdir)
        elif args.require_silo:
            all_errors.append("Silo: full verification required but only "
                              "smoke check available")

    if all_errors:
        print("\nFAILED with %d error(s):" % len(all_errors))
        for e in all_errors[:50]:
            print("  - %s" % e)
        return 1
    print("\nOK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
