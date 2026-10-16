#!/bin/zsh
# 15f-3 acceptance (1): bench.py back to back against the merge base c193cb1.
# Usage: pairs.sh W B S   (W the branch worktree, B a worktree of c193cb1,
# S scratch for quality meshes). Each tree's build-bench/ is built beforehand
# with bench.build(); the runs rebuild it (a no-op). Two balanced batches in
# opposite order: B N N B, then N B B N. Stops if the machine is not on AC.
set -u
W=$1 B=$2 S=$3
PY=$W/.venv/bin/python
D=$W/docs/benchmarks/2026-10-04
cd $W
run() { # label tree baseline-args
  local lab=$1 tree=$2 base=$3
  if ! pmset -g batt | grep -q "AC Power"; then echo "=== NOT ON AC, stop before $lab"; pmset -g batt; exit 3; fi
  echo "=== $lab start $(date -u +%FT%TZ)"; pmset -g batt; sysctl vm.swapusage
  $PY tools/bench.py run --label $lab --tree $tree --hardening on --repeats 5 \
     --mesh-dir $S/meshes/$lab ${=base}
  echo "=== $lab exit $? $(date -u +%FT%TZ)"; pmset -g batt; sysctl vm.swapusage
}
run 15f-3-base-r1 $B "--baseline $W/docs/benchmarks/2026-10-03/24-on-r2"
run 15f-3-r1      $W "--baseline $D/15f-3-base-r1"
run 15f-3-r2      $W "--baseline $D/15f-3-base-r1"
run 15f-3-base-r2 $B "--baseline $D/15f-3-r2"
run 15f-3-r3      $W "--baseline $D/15f-3-base-r2"
run 15f-3-base-r3 $B "--baseline $D/15f-3-r3"
run 15f-3-base-r4 $B "--baseline $D/15f-3-r3"
run 15f-3-r4      $W "--baseline $D/15f-3-base-r4"
echo ALLDONE; pmset -g batt
