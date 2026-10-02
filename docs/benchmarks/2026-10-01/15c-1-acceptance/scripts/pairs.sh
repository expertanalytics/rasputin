#!/bin/sh
# 15c-1 acceptance: back-to-back bench.py pairs, base (a130f7c) then branch.
# Run from the 15c-1 worktree root, under `caffeinate -ims`:
#   PY=<venv python> BASE=<detached worktree of a130f7c> MESH=<scratch dir> \
#     sh docs/benchmarks/2026-10-01/15c-1-acceptance/scripts/pairs.sh 1 2 3
# Each argument is a pair number. Evidence lands in
# docs/benchmarks/<date>/15c-1-acceptance/{base-a130f7c,15c-1}-r<N>/.
set -u
: "${PY:?}" "${BASE:?}" "${MESH:?}"
L=15c-1-acceptance
for n in "$@"; do
  mkdir -p "$MESH/base-r$n/$L" "$MESH/new-r$n/$L"
  echo "== pair $n base $(date +%T)"; pmset -g batt | head -2
  "$PY" tools/bench.py run --label "$L/base-a130f7c-r$n" --tree "$BASE" \
    --mesh-dir "$MESH/base-r$n"
  echo "exit $?"
  echo "== pair $n new $(date +%T)"; pmset -g batt | head -2
  "$PY" tools/bench.py run --label "$L/15c-1-r$n" --tree . \
    --mesh-dir "$MESH/new-r$n" \
    --baseline "docs/benchmarks/$(date +%F)/$L/base-a130f7c-r$n"
  echo "exit $?"
done
pmset -g batt | head -2
