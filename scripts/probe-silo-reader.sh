#!/bin/bash
#
# One-shot diagnostic (intended to run in CI) to decide how to numerically read
# CarpetX Silo plane output. It produces a small Silo plane file and then
# checks, in order of preference:
#   (1) a Silo Python module (import Silo)
#   (2) an apt-installable Silo Python package (python3-silo / silo-python)
#   (3) Silo CLI tools (browser / silodiff / silex)
#   (4) the HDF5 layout, for an h5py-based reader (the .silo files use DB_HDF5)
#
# This script is a temporary aid; remove it (and the CI step that calls it)
# once the Silo reader is wired into verify_planes.py.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CARPETX_REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
CACTUS_DIR="${CACTUS_DIR:-$CARPETX_REPO/../workspace/Cactus}"
EXE="${EXE:-$CACTUS_DIR/exe/cactus_sim}"
TESTDIR="$CARPETX_REPO/TestPlanes/test"
W="$(mktemp -d)"

export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export LD_LIBRARY_PATH="/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial/lib:${LD_LIBRARY_PATH:-}"

echo "########## producing a Silo plane file ##########"
( cd "$W" && "$EXE" "$TESTDIR/planes-single-level.par" ) >/dev/null 2>&1 \
  && echo "  ok" || echo "  (run returned nonzero -- continuing)"

echo "########## (1) Silo Python module ##########"
python3 -c "import Silo, sys; print('  import Silo OK:', Silo.__file__)" 2>&1 | head -2

echo "########## (2) apt packages ##########"
if command -v apt-get >/dev/null 2>&1; then
  apt-get update -qq 2>/dev/null
  echo "  -- apt-cache search silo --"
  apt-cache search silo 2>/dev/null | sed 's/^/    /'
  for pkg in python3-silo silo-python python3-pydbsilo; do
    echo "  -- apt-get install -y $pkg --"
    apt-get install -y "$pkg" 2>&1 | tail -2 | sed 's/^/    /'
  done
  python3 -c "import Silo; print('  import Silo after apt: OK')" 2>&1 | head -1 | sed 's/^/  /'
else
  echo "  (no apt-get)"
fi

echo "########## (3) Silo CLI tools ##########"
for t in browser silex silodiff silofile silock silo2silo; do
  echo "  $t: $(command -v "$t" 2>/dev/null || echo none)"
done
ls /usr/local/bin 2>/dev/null | grep -i silo | sed 's,^,  bin: ,'

echo "########## (4) HDF5 layout via h5py ##########"
python3 "$TESTDIR/verify_planes.py" \
  --parfile "$TESTDIR/planes-single-level.par" \
  --out-dir "$W/planes-single-level" \
  --silo-mode inspect 2>&1 | head -120

echo "########## (5) per-centering orientation dump ##########"
python3 - "$W/planes-single-level" <<'PY'
import sys, glob, os, re, numpy as np
import Silo, h5py
out = sys.argv[1]
sc = lambda v: (v[0] if isinstance(v, (list, tuple)) else v)
files = sorted(glob.glob(os.path.join(out, "*.silo_planes.dir", "*.silo")))
db = Silo.Open(files[0], Silo.DB_READ)
toc = db.GetToc()
seen = set()
with h5py.File(files[0], "r") as h:
    for qv in list(getattr(toc, "qvar_names", [])):
        m = re.search(r"gf(\d\d\d)", qv)
        if not m or m.group(1) in seen:
            continue
        seen.add(m.group(1))
        info = db.GetVarInfo(qv)
        mi = db.GetVarInfo(sc(info["meshid"]))
        V = np.asarray(h[sc(info["value0"])][()])
        c0 = np.asarray(h[sc(mi["coord0"])][()])
        c1 = np.asarray(h[sc(mi["coord1"])][()])
        print("gf%s cent=%s dims1,2=(%s,%s) V.shape=%s len(c0)=%d len(c1)=%d"
              % (m.group(1), sc(info["centering"]),
                 sc(info["dims1"]), sc(info["dims2"]),
                 V.shape, len(c0), len(c1)))
        print("    c0[:2]=%s c1[:2]=%s | V[0,0]=%g dV[0,1]=%g dV[1,0]=%g"
              % (c0[:2].tolist(), c1[:2].tolist(), V[0, 0],
                 V[0, 1] - V[0, 0], V[1, 0] - V[0, 0]))
        print("    mesh min/max idx: a=[%s,%s] b=[%s,%s]"
              % (sc(mi["min_index1"]), sc(mi["max_index1"]),
                 sc(mi["min_index2"]), sc(mi["max_index2"])))
db.Close()
PY

echo "########## (6) numeric Silo read test (--silo-mode verify) ##########"
python3 "$TESTDIR/verify_planes.py" \
  --parfile "$TESTDIR/planes-single-level.par" \
  --out-dir "$W/planes-single-level" \
  --silo-mode verify 2>&1 | tail -30

echo "########## done ##########"
rm -rf "$W"
