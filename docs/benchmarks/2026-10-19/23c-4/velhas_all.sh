#!/bin/zsh
# 23c-4's Velhas runs, through velhas.sh: 10 m, 23c-3 (B, bb5b2fcc) and 23c-4 (W)
# back to back at --pieces 64, 16, 256 (K6: piece files byte-identical).
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
echo VELHASDONE
