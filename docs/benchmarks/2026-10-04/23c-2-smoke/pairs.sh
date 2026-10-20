#!/bin/zsh
# 23c-2: bench.py back to back against its base aacf37c (B, a detached worktree)
# and the branch (N). One balanced batch B N N B (time); full sweep, 5 repeats.
# Stops off AC, or if swap grew by more than 100 MB since the start.
# Usage: pairs.sh W B S
set -u
W=$1 B=$2 S=$3
PY=$W/.venv/bin/python
D=$W/docs/benchmarks/2026-10-04
cd $W
swap() { sysctl -n vm.swapusage | sed -E 's/.*used = ([0-9.]+)M.*/\1/'; }
S0=$(swap)
run() {
  pmset -g batt | grep -q "AC Power" || { echo "NOT ON AC before $1"; exit 3; }
  [ $(echo "$(swap) - $S0 > 100" | bc) = 1 ] && { echo "SWAP GREW before $1: $(swap) MB vs $S0"; exit 4; }
  echo "=== $1 start $(date -u +%FT%TZ) swap $(swap)"; pmset -g batt
  $PY tools/bench.py run --label $1 --tree $2 --hardening on --repeats 5 --mesh-dir $S/m23c/$1 ${=3}
  echo "=== $1 exit $? $(date -u +%FT%TZ) swap $(swap)"; pmset -g batt
}
run 23c2-base-r1 $B ""
run 23c2-r1      $W "--baseline $D/23c2-base-r1"
run 23c2-r2      $W "--baseline $D/23c2-base-r1"
run 23c2-base-r2 $B "--baseline $D/23c2-r2"
echo ALLDONE
