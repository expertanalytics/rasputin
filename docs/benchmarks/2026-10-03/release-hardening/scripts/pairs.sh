#!/bin/bash
# Release vs Release+hardening, back to back, one bench.py run each, in the
# order SEQ (default ABBA-balanced against drift). Each configured build dir is
# swapped in and out of <tree>/build-bench by mv; every one was configured at
# that path, so each is consistent when in place. The variant in place is named
# by build-bench/.variant; a parked one lives at build-park-<variant>.
# Usage: pairs.sh W, with the builds made by build.sh (see there for the order).
# Batch 1 (clang): SEQ="plain fast fast plain plain fast" PREFIX=h.
# Batch 2 (GCC 16): SEQ="gplain gassert gassert gplain gplain gassert" PREFIX=g.
set -u
W=${1:?usage: $0 W (the worktree measured) ...}
PY=${PY:-/Users/skavhaug/projects/rasputin/.venv/bin/python}
OUT=${OUT:-/Users/skavhaug/projects/rasputin_scratch/hardening-2026-10-03}
mkdir -p "$OUT/runs" "$OUT/meshes" "$OUT/logs"
use() {
  cur=$(cat "$W/build-bench/.variant")
  [ "$cur" = "$1" ] && return
  mv "$W/build-bench" "$W/build-park-$cur"
  mv "$W/build-park-$1" "$W/build-bench"
}
n=0
for v in ${SEQ:-plain fast fast plain plain fast}; do
  n=$((n + 1))
  use "$v"
  label="${PREFIX:-h}${n}-$v"
  echo "== $label $(date +%T)" ; pmset -g batt | head -1
  "$PY" "$W/tools/bench.py" run --label "$label" --tree "$W" --out-root "$OUT/runs" \
    --mesh-dir "$OUT/meshes/$label" --repeats "${REPEATS:-5}" > "$OUT/logs/$label.log" 2>&1
  echo "exit $? $(date +%T)"; tail -2 "$OUT/logs/$label.log"
done
echo DONE
