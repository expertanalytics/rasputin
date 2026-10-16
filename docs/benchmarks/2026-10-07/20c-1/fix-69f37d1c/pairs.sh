#!/bin/bash
# 1 m benchmark and thread sweep (tools/bench.py run, defaults: 7908_3 tile and quarter,
# tolerance 1, threads default and 1 to 20, 5 repeats, hardening on), order ABBAAB:
# base = master ed125121 (scratch worktree), head = 20c-1 69f37d1c. Stops off AC.
# Copied from ../pairs.sh; scratch perf-20c1b, three pairs instead of two.
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c1b
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality
cd $W
for spec in base:r1 head:r1 head:r2 base:r2 base:r3 head:r3; do
  who=${spec%%:*}; r=${spec##*:}
  if [ $who = base ]; then tree=$S/base; else tree=$W; fi
  if ! pmset -g batt | grep -q "AC Power"; then echo "$(date -u +%FT%TZ) STOP: off AC"; exit 3; fi
  echo "$(date -u +%FT%TZ) start b1m-$who-$r $(pmset -g batt | tail -1)"
  .venv/bin/python tools/bench.py run --label b1m-$who-$r --tree $tree --out-root $S/out --mesh-dir $S/meshes1m 2>&1 | tail -2
  echo "$(date -u +%FT%TZ) end b1m-$who-$r exit ${PIPESTATUS[0]} $(pmset -g batt | tail -1)"
done
echo DONE
