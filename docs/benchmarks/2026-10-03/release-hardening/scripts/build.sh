#!/bin/bash
# The four Release builds of _core this study measured, in the order and under
# the directory names the other scripts expect. Usage:
#
#   build.sh W clang   # before batch 1: clang FAST parked, clang plain in place
#   SEQ="plain fast fast plain plain fast" PREFIX=h pairs.sh W
#   build.sh W gcc     # before batch 2: clang pair renamed aside, GCC builds
#   SEQ="gplain gassert gassert gplain gplain gassert" PREFIX=g pairs.sh W
#   probes.sh W > raw/probes.txt
#
# W is the worktree measured (a checkout of 390b516 here). Every build is
# configured at W/build-bench with bench.py's own configure line (tools/bench.py,
# build()), plus a cached CMAKE_CXX_FLAGS for the hardened ones: bench.py
# reconfigures without CMAKE_CXX_FLAGS, so the cached value survives. Each build
# is then moved aside to W/build-park-<variant>; pairs.sh moves them back in
# turn. W/build-bench/.variant names the build in place. Python is the main
# checkout's venv, as in the runs.
set -eu
W=$1
STAGE=$2
PY=${PY:-/Users/skavhaug/projects/rasputin/.venv/bin/python}

configure_and_build() {  # variant, then extra cmake args
  local variant=$1
  shift
  cmake -S "$W" -B "$W/build-bench" -DCMAKE_BUILD_TYPE=Release \
    -DRASPUTIN_BUILD_PYTHON=ON -DRASPUTIN_BUILD_TESTS=OFF \
    "-DPython_EXECUTABLE=$PY" -DPYBIND11_FINDPYTHON=ON "$@"
  cmake --build "$W/build-bench" -j --target _core
  echo "$variant" > "$W/build-bench/.variant"
}

case $STAGE in
  clang)
    configure_and_build fast "-DCMAKE_CXX_FLAGS=-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST"
    mv "$W/build-bench" "$W/build-park-fast"
    configure_and_build plain
    ;;
  gcc)
    # The clang pair goes aside under the names probes.sh reads. After batch 1
    # one clang build is in build-bench and the other parked; .variant says which.
    cur=$(cat "$W/build-bench/.variant")
    other=$([ "$cur" = plain ] && echo fast || echo plain)
    mv "$W/build-bench" "$W/build-clang-$cur"
    mv "$W/build-park-$other" "$W/build-clang-$other"
    configure_and_build gassert -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16 \
      "-DCMAKE_CXX_FLAGS=-D_GLIBCXX_ASSERTIONS"
    mv "$W/build-bench" "$W/build-park-gassert"
    configure_and_build gplain -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16
    ;;
  *)
    echo "usage: build.sh W clang|gcc" >&2
    exit 2
    ;;
esac
