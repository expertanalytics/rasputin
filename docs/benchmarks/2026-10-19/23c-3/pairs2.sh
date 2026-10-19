#!/bin/zsh
# 23c-3 at 09c8aad4: bench.py back to back against 23c-2 (B, 321bac56) and master
# (M, cef9f3f0), detached scratch worktrees; W is the branch. Second batch, the order reversed:
# W B B W; full sweep, 5 repeats, hardening on. Stops off AC, or if swap grew
# by more than 300 MB since the start.
# Usage: pairs.sh W B M S   (S: scratch directory for the quality meshes)
set -u
W=$1 B=$2 M=$3 S=$4
PY=$W/.venv/bin/python
O=$W/docs/benchmarks/2026-10-19/23c-3
D=$O/$(date +%F)
cd $W
swap() { sysctl -n vm.swapusage | sed -E 's/.*used = ([0-9.]+)M.*/\1/'; }
S0=$(swap)
run() {
  pmset -g batt | grep -q "AC Power" || { echo "NOT ON AC before $1"; exit 3; }
  [ $(echo "$(swap) - $S0 > 300" | bc) = 1 ] && { echo "SWAP GREW before $1: $(swap) MB vs $S0"; exit 4; }
  echo "=== $1 start $(date -u +%FT%TZ) swap $(swap) load $(sysctl -n vm.loadavg)"; pmset -g batt
  caffeinate -i $PY tools/bench.py run --label $1 --tree $2 --hardening on --repeats 5 --out-root $O --mesh-dir $S/$1 ${=3}
  echo "=== $1 exit $? $(date -u +%FT%TZ) swap $(swap) load $(sysctl -n vm.loadavg)"; pmset -g batt
}
# the CPU's other users, every 15 s, while the batch runs
( while true; do echo "--- $(date -u +%T)"; ps -Ao pcpu,comm -r | sed -n 2,6p; sleep 15; done ) > $O/cpu2.log 2>&1 &
SAMPLER=$!
run head-r3 $W "--baseline $D/base-r2"
run base-r3 $B "--baseline $D/head-r3"
run base-r4 $B "--baseline $D/head-r3"
run head-r4 $W "--baseline $D/base-r4"
kill $SAMPLER
echo ALLDONE
