#!/bin/bash
# Section 6 "At the merge": the probe (docs/increments/30c-probes/dem_bytes.py) on the merge commit 26a5d839,
# a non-editable install in a scratch venv (git archive 26a5d839, then uv pip install ".[codecs]" pytest
# pytest-asyncio hypothesis; Python 3.14.7; base_6c729e97.txt records no Python version). Each mode twice; run 1 is the record,
# run 2 must equal it. Run from the worktree root (its tests/python and probe equal 26a5d839's:
# git diff --stat 26a5d839 HEAD -- tests dem_bytes.py prints nothing).
WT=/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/30c-merge
P=$S/venv/bin/python
R=$WT/docs/benchmarks/2026-10-06/30c-dem-read/raw
B=$WT/docs/increments/30c-probes/base_6c729e97.txt
cd $WT
for n in 1 2; do
  for mode in fixtures mesh; do
    echo "merge run $n $mode before $(date '+%H:%M:%S'): $(pmset -g batt | tr '\n\t' '  ')" >> $S/power.txt
    if [ $mode = fixtures ]; then
      PYTHONPATH=tests/python $P docs/increments/30c-probes/dem_bytes.py fixtures > $S/$mode$n.txt 2> $S/${mode}${n}_stderr.txt
    else
      $P docs/increments/30c-probes/dem_bytes.py mesh numedalslagen skiensvassdraget > $S/$mode$n.txt 2> $S/${mode}${n}_stderr.txt
    fi
    echo "exit $?" >> $S/${mode}${n}_stderr.txt
    echo "merge run $n $mode after $(date '+%H:%M:%S'): $(pmset -g batt | tr '\n\t' '  ')" >> $S/power.txt
  done
  cat $S/fixtures$n.txt $S/mesh$n.txt > $S/run$n.txt
done
