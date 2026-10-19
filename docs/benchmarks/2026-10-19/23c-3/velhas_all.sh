#!/bin/zsh
# Every Velhas run of this acceptance, through velhas.sh, small to large.
# 10 m first, 23c-2 (B) and 23c-3 (W) back to back at --pieces 16, 64, 256
# (K6 and memory per piece); then 23c-3's table, cut at 64 and uncut, 50 to
# 1 m; then 23c-2 at 1 m, 64 pieces (the largest memory against 23c-2's).
# Stops if swap grew by more than 1 GB since the start.
# Usage: velhas_all.sh W B S
set -u
W=$1 B=$2 S=$3
V=${0:A:h}/velhas.sh
swap() { sysctl -n vm.swapusage | sed -E 's/.*used = ([0-9.]+)M.*/\1/'; }
S0=$(swap)
go() {
  [ $(echo "$(swap) - $S0 > 1024" | bc) = 1 ] && { echo "SWAP GREW before $*: $(swap) MB vs $S0"; exit 4; }
  zsh $V "$@" || exit $?
}
for p in 64 16 256; do go $B base $S 10 $p; go $W head $S 10 $p; done
for t in 50 20 10 5 2 1; do
  [ $t != 10 ] && go $W head $S $t 64
  go $W head $S $t 0
done
go $B base $S 1 64
echo VELHASDONE
