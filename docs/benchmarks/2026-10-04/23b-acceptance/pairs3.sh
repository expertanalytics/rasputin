#!/bin/zsh
# 23b acceptance, experiment: E = c4fb2bf + exp-unswitch.diff (scan inner loop unswitched on frozen). Usage: pairs3.sh W B S E
# Usage: pairs.sh W B S
#   W  the branch worktree (its tools/bench.py drives every run)
#   B  a detached worktree of the merge base b4bcdc3
#   S  a scratch directory for the quality meshes
# Each tree's build-bench/ is built once beforehand with bench.build(); the
# runs rebuild it (a no-op). Order B N N B: two pairs, balanced against drift.
# A run is skipped (and the batch stops) when the machine is not on AC.
set -u
W=$1 B=$2 S=$3 E=$4
PY=$W/.venv/bin/python
D=$W/docs/benchmarks/2026-10-04
cd $W
run() { # label tree baseline-args
  local lab=$1 tree=$2 base=$3
  if ! pmset -g batt | grep -q "AC Power"; then echo "=== NOT ON AC, stop before $lab"; pmset -g batt; exit 3; fi
  echo "=== $lab start $(date -u +%FT%TZ)"; pmset -g batt
  $PY tools/bench.py run --label $lab --tree $tree --hardening on --repeats 5 \
     --mesh-dir $S/meshes/$lab ${=base}
  echo "=== $lab exit $? $(date -u +%FT%TZ)"; pmset -g batt
}
run 23b-exp-unswitch    $E "--baseline $D/23b-r4"
run 23b-base-b4bcdc3-r5 $B "--baseline $D/23b-exp-unswitch"
run 23b-r5              $W "--baseline $D/23b-base-b4bcdc3-r5"
run 23b-exp-unswitch-r2 $E "--baseline $D/23b-base-b4bcdc3-r5"
echo ALLDONE; pmset -g batt
