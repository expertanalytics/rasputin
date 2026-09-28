#!/bin/zsh
set -u
cd /Users/skavhaug/projects/rasputin
SP=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad
WT=$SP/wt-master-dc372a5
D=docs/benchmarks/2026-09-28/16b0-acceptance
COMMON=(--mesh-dir $SP/p16b0/meshes --domain docs/benchmarks/2026-09-26/quarter.geojson --domain tile --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --repeats 5 --threads 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
for r in 1 2; do
  sfx=""; [ $r = 2 ] && sfx="-r2"
  echo "== base r$r $(date +%T)"; pmset -g batt
  .venv/bin/python tools/bench.py run --label 16b0-acceptance/base-dc372a5$sfx --tree $WT $COMMON; echo "exit $?"
  echo "== 16b0 r$r $(date +%T)"; pmset -g batt
  .venv/bin/python tools/bench.py run --label 16b0-acceptance/16b0$sfx --baseline $D/base-dc372a5$sfx $COMMON; echo "exit $?"
done
echo "== done $(date +%T)"; pmset -g batt
