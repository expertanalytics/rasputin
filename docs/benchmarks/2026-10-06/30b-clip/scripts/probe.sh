#!/bin/bash
# The design's byte-identical probe (docs/increments/30b-probes/clip_bytes.py) on the branch install,
# run after the timing runs and never beside them, then compared with base_5e2fbe0.txt (section 6's check).
# The branch venv needs pytest, pytest-asyncio and hypothesis for `fixtures` (installed into the scratch venv).
WT=/Users/skavhaug/projects/rasputin/.claude/worktrees/clip-speed
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/30b
R=$WT/docs/benchmarks/2026-10-06/30b-clip/raw
P=$S/venv-branch/bin/python
cd $WT
: > $R/probe_power.txt
pmset -g batt | head -1 | sed 's/^/before fixtures: /' >> $R/probe_power.txt
PYTHONPATH=tests/python $P docs/increments/30b-probes/clip_bytes.py fixtures > $R/probe_fixtures_branch.txt 2> $R/probe_fixtures_branch_stderr.txt
echo "exit $?" >> $R/probe_fixtures_branch_stderr.txt
pmset -g batt | head -1 | sed 's/^/after fixtures: /' >> $R/probe_power.txt
pmset -g batt | head -1 | sed 's/^/before mesh: /' >> $R/probe_power.txt
$P docs/increments/30b-probes/clip_bytes.py mesh numedalslagen skiensvassdraget > $R/probe_mesh_branch.txt 2> $R/probe_mesh_branch_stderr.txt
echo "exit $?" >> $R/probe_mesh_branch_stderr.txt
pmset -g batt | head -1 | sed 's/^/after mesh: /' >> $R/probe_power.txt
B=docs/increments/30b-probes/base_5e2fbe0.txt
cat $R/probe_fixtures_branch.txt $R/probe_mesh_branch.txt > $S/run.txt
{
  echo "== base lines missing or changed on the branch (must be empty):"
  comm -23 <(grep -E '^(fixture|mesh) ' $B | sort) <(grep -E '^(fixture|mesh) ' $S/run.txt | sort)
  echo "== lines on the branch the base lacks (allowed: the red suite's 11, section 6):"
  comm -13 <(grep -E '^(fixture|mesh) ' $B | sort) <(grep -E '^(fixture|mesh) ' $S/run.txt | sort)
  echo "== counts: base $(grep -cE '^(fixture|mesh) ' $B), branch $(grep -cE '^(fixture|mesh) ' $S/run.txt)"
} > $R/probe_compare.txt
cat $R/probe_compare.txt
