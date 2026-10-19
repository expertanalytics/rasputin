#!/bin/zsh
# 23c-4 at e970bde4 (W) against 23c-3 at bb5b2fcc (B, a detached scratch worktree):
# bench.py W B W, 5 repeats, full sweep, hardening on. head-r1 is judged against
# 23c-3's accepted head-r4 run. Stops off AC, or if swap grew by more than 300 MB.
# Usage: pairs.sh W B S   (S: scratch directory for the quality meshes)
set -u
W=$1 B=$2 S=$3
PY=$W/.venv/bin/python
O=$W/docs/benchmarks/2026-10-19/23c-4
D=$O/$(date +%F)
P3=$W/docs/benchmarks/2026-10-19/23c-3/2026-10-19/head-r4
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
run head-r1 $W "--baseline $P3"
run base-r1 $B "--baseline $D/head-r1"
run head-r2 $W "--baseline $D/base-r1"
echo ALLDONE
