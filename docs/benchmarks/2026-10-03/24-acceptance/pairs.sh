#!/bin/zsh
# Increment 24's paired acceptance (docs/increments/24-release-hardening.md §7).
# Usage: pairs.sh W M OFF S
#   W   the branch worktree (its tools/bench.py drives every run; built ON)
#   M   a detached worktree of the previous merge (origin/master's tip)
#   OFF a detached worktree of the branch's commit, built OFF
#   S   a scratch directory for the quality meshes
# Each tree's build-bench/ is built once beforehand with bench.build(); the
# runs below rebuild it (a no-op) with the same --hardening value.
# Order M OFF ON ON OFF M: pairs (M,OFF) and (OFF,M) for "nothing else moved",
# (OFF,ON) and (ON,OFF) for the hardening cost.
set -u
W=$1 M=$2 OFF=$3 S=$4
PY=/Users/skavhaug/projects/rasputin/.venv/bin/python
D=$W/docs/benchmarks/2026-10-03
cd $W
run() { # label tree mode baseline-args
  local lab=$1 tree=$2 mode=$3 base=$4
  echo "=== $lab start $(date -u +%FT%TZ)"; pmset -g batt
  $PY tools/bench.py run --label $lab --tree $tree --hardening $mode --repeats 5 \
     --mesh-dir $S/meshes/$lab ${=base}
  echo "=== $lab exit $? $(date -u +%FT%TZ)"
}
run 24-base-6cdc8cc    $M   off ""
run 24-off             $OFF off "--baseline $D/24-base-6cdc8cc"
run 24-on              $W   on  "--baseline $D/24-off"
run 24-on-r2           $W   on  "--baseline $D/24-on"
run 24-off-r2          $OFF off "--baseline $D/24-base-6cdc8cc"
run 24-base-6cdc8cc-r2 $M   off "--baseline $D/24-off-r2"
echo ALLDONE; pmset -g batt
