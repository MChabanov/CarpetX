#!/bin/bash
#
# Run the TestPlanes parfiles with the freshly built executable and verify the
# 2D plane output (Silo + openPMD) analytically. Intended to run after a build,
# both in CI (.github/workflows/ci.yml) and locally.
#
# Environment:
#   CACTUS_DIR   Cactus source/run directory (default: $PWD/../workspace/Cactus)
#   EXE          Path to the cactus executable (default: $CACTUS_DIR/exe/cactus_sim)
#   CARPETX_DIR  CarpetX repo inside Cactus (default: $CACTUS_DIR/repos/CarpetX)
#   MPIRUN       MPI launcher (default: mpirun); set to empty to run serially.

set -eu

CACTUS_DIR="${CACTUS_DIR:-$PWD/../workspace/Cactus}"
EXE="${EXE:-$CACTUS_DIR/exe/cactus_sim}"
CARPETX_DIR="${CARPETX_DIR:-$CACTUS_DIR/repos/CarpetX}"
MPIRUN="${MPIRUN-mpiexec}"

TESTDIR="$CARPETX_DIR/TestPlanes/test"
GOLDEN="$TESTDIR/golden"
WORKDIR="${WORKDIR:-$CACTUS_DIR/test-planes-run}"

# OpenMPI refuses to run as root unless told otherwise (containers run as root).
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export LD_LIBRARY_PATH="/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial/lib:${LD_LIBRARY_PATH:-}"

if [ ! -x "$EXE" ]; then
  echo "✗ executable not found or not executable: $EXE" >&2
  exit 1
fi

# Ensure the Python readers are present. openpmd_api reads the ADIOS2/BP5 plane
# files; this needs network access (fine on GitHub runners, not on HPC login
# nodes). Errors are shown (not hidden), and PEP 668 "externally-managed"
# environments are handled with --break-system-packages.
echo "python: $(command -v python3)"; python3 --version || true
python3 -m pip --version 2>/dev/null || echo "  (no pip module)"

ensure_py() {
  mod="$1"; pkg="$2"
  if python3 -c "import $mod" 2>/dev/null; then
    echo "  $mod: already importable"
    return 0
  fi
  echo "  installing $pkg ..."
  python3 -m pip install "$pkg" \
    || python3 -m pip install --break-system-packages "$pkg" \
    || python3 -m pip install --user --break-system-packages "$pkg" \
    || true
  python3 -c "import $mod" 2>/dev/null
}

ensure_py openpmd_api openpmd-api && HAVE_OPENPMD=1 || HAVE_OPENPMD=0
ensure_py h5py h5py || true

if [ "$HAVE_OPENPMD" -ne 1 ]; then
  echo "✗ openpmd_api is required for rigorous plane verification but could" >&2
  echo "  not be imported or installed (see pip output above)." >&2
  exit 1
fi

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# parfile : nprocs
run_and_verify() {
  par_base="$1"
  nprocs="$2"
  par="$TESTDIR/$par_base.par"
  echo "================================================================"
  echo "Running $par_base on $nprocs proc(s)"
  rm -rf "$par_base"
  if [ -n "$MPIRUN" ] && [ "$nprocs" -gt 1 ]; then
    "$MPIRUN" -np "$nprocs" "$EXE" "$par"
  else
    "$EXE" "$par"
  fi
  echo "Verifying $par_base"
  python3 "$TESTDIR/verify_planes.py" \
    --parfile "$par" \
    --out-dir "$par_base" \
    --golden-dir "$GOLDEN" \
    --require-openpmd \
    --silo-mode "${SILO_MODE:-auto}"
}

# Run on >1 rank where it matters: the single-level and AMR cases exercise the
# Silo MPI gather / multi-file metafile and the openPMD per-rank collective
# chunks (the grids are forced into several boxes via max_grid_size). The
# edge-case run only needs to exercise parsing/snapping, so it runs serially.
run_and_verify planes-single-level 2
run_and_verify planes-amr 2
run_and_verify planes-edge-cases 1

echo "================================================================"
echo "✓ plane verification passed"
