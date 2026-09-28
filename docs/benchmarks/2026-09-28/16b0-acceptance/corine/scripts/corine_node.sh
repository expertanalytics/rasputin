#!/bin/zsh
cd /Users/skavhaug/projects/rasputin
SP=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad
NEW=/Users/skavhaug/projects/rasputin/build-bench/pkg
OLD=$SP/wt-master-dc372a5/build-bench/pkg
PY=.venv/bin/python; D=$SP/p16b0/drv.py
run() { echo "# $(date +%T) $*"; pmset -g batt | tail -1; $PY $D "$@" 2>&1 | grep '^RESULT\|Error\|error' ; }
for s in 12 24 36 48 96 144; do for L in A B; do run --pkg $NEW node --layout $L --side $s --repeats 3; done; done
for s in 12 24 36 48; do run --pkg $OLD node --layout A --side $s --repeats 1; done
for s in 12 24 48; do run --pkg $OLD node --layout B --side $s --repeats 1; done
for n in 300 600 1200; do run --pkg $NEW ladder --n $n --repeats 3; done
for n in 300 600; do run --pkg $OLD ladder --n $n --repeats 1; done
for n in 1000 4000 16000 64000; do run --pkg $NEW comb --n $n --repeats 3; done
for n in 1000 4000; do run --pkg $OLD comb --n $n --repeats 1; done
echo "# done $(date +%T)"; pmset -g batt | tail -1
