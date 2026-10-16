#!/bin/bash
# The reprojected Velhas run, base (618328b, audit-crs) then head (1a84ea9, audit-lattice), back to back. Stops off AC.
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad
W=/Users/skavhaug/projects/rasputin/.claude/worktrees
D=/Users/skavhaug/projects/rasputin_data/sao_francisco_piece
mkdir -p $S/vmeshes
cd $W/audit-lattice
for who in ${WHO:-base head}; do
  if [ $who = base ]; then tree=$W/audit-crs; else tree=$W/audit-lattice; fi
  if ! pmset -g batt | grep -q "AC Power"; then echo "$(date -u +%FT%TZ) STOP: off AC"; exit 3; fi
  echo "$(date -u +%FT%TZ) start vel-$who-r1 $(pmset -g batt | tail -1)"
  .venv/bin/python tools/bench.py run --label vel-$who-r1 --tree $tree --out-root $S/out --mesh-dir $S/vmeshes \
    --dem $D/bho2017_5k_76949_anadem_window_epsg4674.tif \
    --domain $D/bho2017_5k_76949_outline_epsg4674.geojson \
    --tolerance 10 -- --out-crs EPSG:31983 2>&1 | tail -4
  echo "$(date -u +%FT%TZ) end vel-$who-r1 exit ${PIPESTATUS[0]}"
done
echo DONE
