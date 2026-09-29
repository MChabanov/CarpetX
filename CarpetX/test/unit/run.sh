#!/bin/bash
# Build and run the standalone unit tests for the vector-potential
# prolongation.  Needs nothing but a C++17 compiler; see README.md.
set -e
cd "$(dirname "$0")"
CXX=${CXX:-c++}
status=0
for t in test_vecpot_kernels test_vecpot_interp test_vecpot_selftest; do
  echo "=== $t ==="
  $CXX -std=c++17 -O2 -Wall -Wextra -I shim -o "$t" "$t.cxx"
  ./"$t" || status=1
  echo
done
if [ $status -eq 0 ]; then echo "all unit tests passed"; else echo "UNIT TESTS FAILED" >&2; fi
exit $status
