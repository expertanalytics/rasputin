#!/bin/zsh
# Three bench.py pairs, base (origin/master 14f5fe3 in a detached worktree) then 16b-1/2.
set -u
cd /Users/skavhaug/projects/rasputin
SP=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad
WT=$SP/wt-master-14f5fe3
D=docs/benchmarks/2026-09-28/16b12-acceptance
COMMON=(--mesh-dir $SP/p16b12/meshes --domain docs/benchmarks/2026-09-26/quarter.geojson --domain tile --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --repeats 5 --threads 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
for r in ${PAIRS:-1 2 3}; do
  sfx=""; [ $r != 1 ] && sfx="-r$r"
  echo "== base r$r $(date +%T)"; pmset -g batt
  .venv/bin/python tools/bench.py run --label 16b12-acceptance/base-14f5fe3$sfx --tree $WT $COMMON; echo "exit $?"
  echo "== 16b12 r$r $(date +%T)"; pmset -g batt
  .venv/bin/python tools/bench.py run --label 16b12-acceptance/16b12$sfx --baseline $D/base-14f5fe3$sfx $COMMON; echo "exit $?"
done
echo "== done $(date +%T)"; pmset -g batt
