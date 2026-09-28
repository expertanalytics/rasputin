#!/bin/zsh
cd /Users/skavhaug/projects/rasputin
SP=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad
NEW=/Users/skavhaug/projects/rasputin/build-bench/pkg
PY=.venv/bin/python; D=$SP/p16b0/drv.py
run() { echo "# $(date +%T) $*"; pmset -g batt | tail -1; /usr/bin/time -l $PY $D "$@" 2>&1 | grep '^RESULT\|Error\|error\|maximum resident' ; }
for tol in 1 10; do
  run --pkg $NEW mesh --side 48 --tol $tol --features none --min-angle 25
  run --pkg $NEW mesh --side 48 --tol $tol --features A --min-angle 25
  run --pkg $NEW mesh --side 48 --tol $tol --features A --min-angle 0
  run --pkg $NEW mesh --side 48 --tol $tol --features none --min-angle 0
done
echo "# done $(date +%T)"; pmset -g batt | tail -1
