#!/bin/zsh
# 26b acceptance: Velhas piece at 10 and 5 m, with and without --land-cover, 3 runs each.
set -e
WT=/Users/skavhaug/projects/rasputin/.claude/worktrees/mapbiomas-26b
DATA=/Users/skavhaug/projects/rasputin_data
EV=$WT/docs/benchmarks/2026-10-12/26b-fractions/raw
OUT=/Users/skavhaug/projects/rasputin_scratch/sao_francisco_piece_meshes/26b-fractions
DEM=$DATA/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif
DOMAIN=$DATA/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson
RASPUTIN=$WT/.venv/bin/rasputin
echo "start $(date -u +%FT%TZ)" >> $EV/pmset.log
pmset -g batt >> $EV/pmset.log
for tol in 10 5; do
  for lc in on off; do
    for r in 1 2 3; do
      tag=velhas_t${tol}_lc${lc}_r${r}
      extra=()
      if [[ $lc == on ]]; then
        extra=(--land-cover mapbiomas-c11 --land-cover-year 2025 --cache $DATA/cache)
      fi
      echo "== $tag $(date -u +%FT%TZ)" >> $EV/pmset.log
      /usr/bin/time -l $RASPUTIN mesh --dem $DEM --domain $DOMAIN --tolerance $tol $extra \
        --binary --out $OUT/velhas_t${tol}_lc${lc}.bin.vtk --stats $EV/$tag.stats.md \
        > $EV/$tag.log 2>&1
      echo "exit $? $tag" >> $EV/pmset.log
    done
  done
done
pmset -g batt >> $EV/pmset.log
echo "end $(date -u +%FT%TZ)" >> $EV/pmset.log
