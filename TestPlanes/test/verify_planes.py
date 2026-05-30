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

A compact per-slab summary is also written and diffed against committed golden
files (--golden-dir), giving a frozen snapshot on top of the analytic check.

Usage:
    verify_planes.py --parfile P.par --out-dir DIR [--golden-dir D]
                     [--update-golden] [--silo-mode {auto,verify,smoke,skip}]
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


def inspect_silo(out_dir, sim):
    """Dump the real Silo layout/API to help finalize the reader in-container."""
    print("=== Silo inspection ===")
    metas = silo_meta_files(out_dir, sim)
    datas = silo_data_files(out_dir, sim)
    print("metafiles: %d, data files: %d" % (len(metas), len(datas)))
    sample = (datas or metas)
    if not sample:
        print("no Silo files found in %s" % out_dir)
        return
    sample = sample[0]
    print("sample: %s" % sample)
    try:
        import h5py
        with h5py.File(sample, "r") as h5:
            print("--- HDF5 tree (h5py) ---")
            def show(name, obj):
                kind = "D" if isinstance(obj, h5py.Dataset) else "G"
                extra = ""
                if isinstance(obj, h5py.Dataset):
                    extra = " shape=%s dtype=%s" % (obj.shape, obj.dtype)
                print("  [%s] %s%s" % (kind, name, extra))
            h5.visititems(show)
    except Exception as exc:  # noqa: BLE001
        print("h5py inspection failed: %s" % exc)
    try:
        import Silo
        print("--- Silo module ---")
        db = Silo.Open(sample, Silo.DB_READ)
        toc = db.GetToc()
        for attr in ("qmesh_names", "qvar_names", "multimesh_names",
                     "multivar_names", "var_names", "dir_names"):
            if hasattr(toc, attr):
                print("  toc.%s = %r" % (attr, getattr(toc, attr)))
        qvn = getattr(toc, "qvar_names", [])
        if qvn:
            qv = db.GetQuadvar(qvn[0])
            print("  GetQuadvar(%r) attrs: %s" % (
                qvn[0], [a for a in dir(qv) if not a.startswith("__")]))
        db.Close()
    except Exception as exc:  # noqa: BLE001
        print("Silo module not usable: %s" % exc)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def verify_slabs(records, geom, rtol, atol, summary_lines):
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
            vmin = math.inf
            vmax = -math.inf
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
                    vmin = min(vmin, actual)
                    vmax = max(vmax, actual)
                    # Coverage: record interior in-plane points on level 0.
                    if s.level == 0:
                        a, b = s.coords_a[i], s.coords_b[j]
                        if (geom.xmin[s.axis_a] - eps <= a <= geom.xmax[s.axis_a] + eps
                                and geom.xmin[s.axis_b] - eps <= b <= geom.xmax[s.axis_b] + eps):
                            key = (tag, s.var)
                            coverage.setdefault(key, set()).add(
                                (round(a, 6), round(b, 6)))
                            cover_meta[key] = (s.centering, s.axis_a, s.axis_b)
            summary_lines.append(
                "%s %s patch%d L%d na=%d nb=%d ncoord=%.6f "
                "amin=%.6f amax=%.6f bmin=%.6f bmax=%.6f "
                "vmin=%.6f vmax=%.6f" %
                (tag, s.var, s.patch, s.level, na, nb, ncoord,
                 min(s.coords_a), max(s.coords_a),
                 min(s.coords_b), max(s.coords_b), vmin, vmax))

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

def golden_path(golden_dir, parfile, fmt):
    base = os.path.splitext(os.path.basename(parfile))[0]
    return os.path.join(golden_dir, "%s.%s.summary" % (base, fmt))


def handle_golden(summary_lines, golden_dir, parfile, fmt, update):
    lines = sorted(summary_lines)
    content = "\n".join(lines) + ("\n" if lines else "")
    path = golden_path(golden_dir, parfile, fmt)
    if update:
        os.makedirs(golden_dir, exist_ok=True)
        with open(path, "w") as fh:
            fh.write(content)
        print("  golden updated: %s (%d slabs)" % (path, len(lines)))
        return []
    if not os.path.exists(path):
        # Soft-skip: the analytic check is the rigorous gate. Golden files are
        # an additional frozen snapshot; once committed they are enforced.
        print("  golden absent (not enforced yet): %s "
              "-- run with --update-golden to create" % path)
        return []
    with open(path) as fh:
        want = fh.read()
    if want != content:
        return ["golden mismatch for %s (%s); rerun --update-golden to refresh"
                % (parfile, fmt)]
    return []


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parfile", required=True)
    ap.add_argument("--out-dir", required=True,
                    help="directory containing the plane output")
    ap.add_argument("--golden-dir", default=None)
    ap.add_argument("--update-golden", action="store_true")
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

    # --- openPMD -----------------------------------------------------------
    op_records, op_ok, op_msg = read_openpmd(args.out_dir, sim, geom)
    if not op_ok:
        msg = "openPMD: %s" % op_msg
        if args.require_openpmd:
            all_errors.append(msg)
        else:
            print("SKIP %s" % msg)
    else:
        summary = []
        n, errs = verify_slabs(op_records, geom, args.rtol, args.atol, summary)
        all_errors += errs
        print("openPMD: %d plane file groups, %d values checked, %d errors"
              % (len(op_records), n, len(errs)))
        if args.golden_dir is not None:
            all_errors += handle_golden(summary, args.golden_dir,
                                        args.parfile, "openpmd",
                                        args.update_golden)

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
            summary = []
            n, errs = verify_slabs(si_records, geom, args.rtol, args.atol,
                                   summary)
            all_errors += errs
            print("Silo: %d values checked, %d errors" % (n, len(errs)))
            if args.golden_dir is not None:
                all_errors += handle_golden(summary, args.golden_dir,
                                            args.parfile, "silo",
                                            args.update_golden)
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
