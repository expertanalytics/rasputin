#!/bin/bash
# The design's byte-identical probe (docs/increments/30c-probes/dem_bytes.py), run after the timing runs and
# never beside them, on each side's install (each venv's own python; pytest, pytest-asyncio and hypothesis
# installed into both scratch venvs), then compared with base_6c729e97.txt by section 6's two `comm` commands.
# The base side is run too, to show this base install reproduces the committed base file.
# Run from the worktree root, so tests/python is the branch's (section 6).
WT=/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/30c
R=$WT/docs/benchmarks/2026-10-06/30c-dem-read/raw
B=$WT/docs/increments/30c-probes/base_6c729e97.txt
cd $WT
: > $R/probe_power.txt
: > $R/probe_compare.txt
for side in base branch; do
  P=$S/venv-$side/bin/python
  for mode in fixtures mesh; do
    echo "$side $mode before $(date '+%H:%M:%S'): $(pmset -g batt | tr '\n\t' '  ')" >> $R/probe_power.txt
    if [ $mode = fixtures ]; then
      PYTHONPATH=tests/python $P docs/increments/30c-probes/dem_bytes.py fixtures > $R/probe_fixtures_$side.txt 2> $R/probe_fixtures_${side}_stderr.txt
    else
      $P docs/increments/30c-probes/dem_bytes.py mesh numedalslagen skiensvassdraget > $R/probe_mesh_$side.txt 2> $R/probe_mesh_${side}_stderr.txt
    fi
    echo "exit $?" >> $R/probe_${mode}_${side}_stderr.txt
    echo "$side $mode after $(date '+%H:%M:%S'): $(pmset -g batt | tr '\n\t' '  ')" >> $R/probe_power.txt
  done
  cat $R/probe_fixtures_$side.txt $R/probe_mesh_$side.txt > $S/run_$side.txt
  {
    echo "== $side: $(head -1 $R/probe_fixtures_$side.txt | sed "s|$S/||") / $(head -1 $R/probe_mesh_$side.txt | sed "s|$S/||")"
    echo "== $side: base lines missing or changed (must be empty):"
    comm -23 <(grep -E '^(fixture|mesh) ' $B | sort) <(grep -E '^(fixture|mesh) ' $S/run_$side.txt | sort)
    echo "== $side: lines the base file lacks (allowed on the branch: the red suite's new tests, section 6):"
    comm -13 <(grep -E '^(fixture|mesh) ' $B | sort) <(grep -E '^(fixture|mesh) ' $S/run_$side.txt | sort)
    echo "== $side: counts: base file $(grep -cE '^(fixture|mesh) ' $B), this run $(grep -cE '^(fixture|mesh) ' $S/run_$side.txt); $(grep '^fixtures:' $R/probe_fixtures_$side.txt); pytest summary: $(grep -E '^[0-9]+ (passed|failed)' $R/probe_fixtures_$side.txt | tail -1)"
    echo
  } >> $R/probe_compare.txt
done
cat $R/probe_compare.txt
