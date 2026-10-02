#!/bin/sh
# 15c-1 acceptance: back-to-back bench.py pairs in reverse order, branch first,
# so that a drift over the session (battery level, heat) does not always land
# on the branch. Same environment as pairs.sh:
#   PY=... BASE=... MESH=... sh .../scripts/pairs_reversed.sh 4 5
# Each branch run is judged afterwards with `bench.py compare`, since its
# baseline does not exist yet when it runs (its own verdict is NO BASELINE).
set -u
: "${PY:?}" "${BASE:?}" "${MESH:?}"
L=15c-1-acceptance
D="docs/benchmarks/$(date +%F)/$L"
for n in "$@"; do
  mkdir -p "$MESH/base-r$n/$L" "$MESH/new-r$n/$L"
  echo "== pair $n new $(date +%T)"; pmset -g batt | head -2
  "$PY" tools/bench.py run --label "$L/15c-1-r$n" --tree . --mesh-dir "$MESH/new-r$n"
  echo "exit $?"
  echo "== pair $n base $(date +%T)"; pmset -g batt | head -2
  "$PY" tools/bench.py run --label "$L/base-a130f7c-r$n" --tree "$BASE" \
    --mesh-dir "$MESH/base-r$n"
  echo "exit $?"
  echo "== pair $n compare"
  "$PY" tools/bench.py compare "$D/15c-1-r$n" --baseline "$D/base-a130f7c-r$n"
  echo "exit $?"
done
pmset -g batt | head -2
