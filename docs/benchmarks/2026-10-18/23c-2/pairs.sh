#!/bin/zsh
# 23c-2 at 4e140864: bench.py back to back against master cb713c87 (B, a detached
# worktree) and the branch (W). One balanced batch B W W B; full sweep, 5 repeats,
# hardening on. Stops off AC, or if swap grew by more than 100 MB since the start.
# Usage: pairs.sh W B S   (S: scratch directory for the quality meshes)
set -u
W=$1 B=$2 S=$3
PY=$W/.venv/bin/python
O=$W/docs/benchmarks/2026-10-18/23c-2
D=$O/$(date +%F)
cd $W
swap() { sysctl -n vm.swapusage | sed -E 's/.*used = ([0-9.]+)M.*/\1/'; }
S0=$(swap)
run() {
  pmset -g batt | grep -q "AC Power" || { echo "NOT ON AC before $1"; exit 3; }
  [ $(echo "$(swap) - $S0 > 100" | bc) = 1 ] && { echo "SWAP GREW before $1: $(swap) MB vs $S0"; exit 4; }
  echo "=== $1 start $(date -u +%FT%TZ) swap $(swap)"; pmset -g batt
  $PY tools/bench.py run --label $1 --tree $2 --hardening on --repeats 5 --out-root $O --mesh-dir $S/m23c2/$1 ${=3}
  echo "=== $1 exit $? $(date -u +%FT%TZ) swap $(swap)"; pmset -g batt
}
run base-r1 $B ""
run head-r1 $W "--baseline $D/base-r1"
run head-r2 $W "--baseline $D/base-r1"
run base-r2 $B "--baseline $D/head-r2"
echo ALLDONE
