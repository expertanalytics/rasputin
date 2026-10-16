#!/bin/zsh
# Increment 20 fix: bench.py back to back against master 2060f14 (B, a detached worktree)
# and the branch (N). Batch 1 B N N B, batch 2 N B B N; full sweep, 5 repeats.
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
  $PY tools/bench.py run --label $1 --tree $2 --hardening on --repeats 5 --mesh-dir $S/m20/$1 ${=3}
  echo "=== $1 exit $? $(date -u +%FT%TZ) swap $(swap)"; pmset -g batt
}
run 20nd-base-r1 $B ""
run 20nd-r1      $W "--baseline $D/20nd-base-r1"
run 20nd-r2      $W "--baseline $D/20nd-base-r1"
run 20nd-base-r2 $B "--baseline $D/20nd-r2"
run 20nd-r3      $W "--baseline $D/20nd-base-r2"
run 20nd-base-r3 $B "--baseline $D/20nd-r3"
run 20nd-base-r4 $B "--baseline $D/20nd-r3"
run 20nd-r4      $W "--baseline $D/20nd-base-r4"
echo ALLDONE
