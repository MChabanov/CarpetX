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
MPIRUN="${MPIRUN-mpirun}"

TESTDIR="$CARPETX_DIR/TestPlanes/test"
GOLDEN="$TESTDIR/golden"
WORKDIR="${WORKDIR:-$CACTUS_DIR/test-planes-run}"

# OpenMPI refuses to run as root unless told otherwise (containers run as root).
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH:-}"

if [ ! -x "$EXE" ]; then
  echo "✗ executable not found or not executable: $EXE" >&2
  exit 1
fi

# Best-effort: ensure the Python readers are present.
ensure_py() {
  python3 -c "import $1" 2>/dev/null && return 0
  echo "  (installing python module $2)"
  pip3 install --quiet "$2" 2>/dev/null || true
  python3 -c "import $1" 2>/dev/null
}
HAVE_OPENPMD=0
ensure_py openpmd_api openpmd-api && HAVE_OPENPMD=1
ensure_py h5py h5py || true

if [ "$HAVE_OPENPMD" -ne 1 ]; then
  echo "✗ openpmd_api is required for rigorous plane verification but is not" >&2
  echo "  importable and could not be installed." >&2
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
    --silo-mode "${SILO_MODE:-auto}" \
    ${UPDATE_GOLDEN:+--update-golden}
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
