#!/bin/zsh
# 15f-3 acceptance (3): the Velhas piece, base c193cb1 and 15f-3 back to back,
# B N N B, through velhas.py run (20/10/5 m, 2 timed runs + 1 ascii each).
# Usage: velhas.sh W B S
set -u
W=$1 B=$2 S=$3
PY=$W/.venv/bin/python
V=$W/docs/benchmarks/2026-10-04/15f-3-acceptance
geo() { # label tree
  if ! pmset -g batt | grep -q "AC Power"; then echo "NOT ON AC before $1"; pmset -g batt; exit 3; fi
  echo "=== $1 start $(date -u +%FT%TZ) $(sysctl -n vm.swapusage)"
  $PY $V/velhas.py run $1 $2 $V/velhas/$1 $S/velhas/$1 --tolerances 20,10,5 --repeats 2
  echo "=== $1 exit $? $(date -u +%FT%TZ) $(sysctl -n vm.swapusage)"
}
geo base-r1 $B; geo 15f3-r1 $W; geo 15f3-r2 $W; geo base-r2 $B
echo ALLDONE
