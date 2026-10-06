#!/bin/bash
# ABBA pairs: base = 618328b (audit-crs), head = 1a84ea9 (audit-lattice). Stops off AC.
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad
W=/Users/skavhaug/projects/rasputin/.claude/worktrees
cd $W/audit-lattice
for spec in base:r1 head:r1 head:r2 base:r2 head:r3 base:r3 base:r4 head:r4; do
  who=${spec%%:*}; r=${spec##*:}
  if [ $who = base ]; then tree=$W/audit-crs; else tree=$W/audit-lattice; fi
  if ! pmset -g batt | grep -q "AC Power"; then echo "$(date -u +%FT%TZ) STOP: off AC"; exit 3; fi
  echo "$(date -u +%FT%TZ) start pra-$who-$r $(pmset -g batt | tail -1)"
  .venv/bin/python tools/bench.py run --label pra-$who-$r --tree $tree --out-root $S/out --mesh-dir $S/meshes 2>&1 | tail -2
  echo "$(date -u +%FT%TZ) end pra-$who-$r exit ${PIPESTATUS[0]}"
done
echo DONE
