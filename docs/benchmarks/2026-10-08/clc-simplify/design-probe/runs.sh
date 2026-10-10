#!/bin/zsh
set -e
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/clc-simplify
O=/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify/proto
D=/Users/skavhaug/projects/rasputin_data
K=/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach
# b<band>.geojson in $O come from: python border_apsc.py <band> $O/b<band>.geojson
export RASPUTIN_DATA=$D
pmset -g batt | head -1
for b in 0 20 50 100; do
  for spec in "min 0" "50 25" "50 20" "50 15" "50 0" "20 25" "20 0" "10 25" "10 0"; do
    t=${spec% *}; a=${spec#* }
    if [[ $t == min ]]; then tol=(--tolerance 1e6 --start-min-angle 0); n=b${b}_minimal
    else tol=(--tolerance $t --start-min-angle $a); n=b${b}_tol${t}_a${a}; fi
    s=$(date +%s.%N)
    $W/.venv/bin/rasputin mesh --dem glo30 --cache $D/germany_glo30_cache --out-crs EPSG:25832 \
      --domain $K/fused_catchment_glo30.geojson --features $O/b$b.geojson --features-crs EPSG:25832 \
      --features-map corine --features-tolerance 0 $tol \
      --out $O/$n.vtk --stats $O/${n}_stats.md >/dev/null
    tri=$(grep '^| output triangles' $O/${n}_stats.md | awk -F'|' '{print $3}')
    q=$(grep '^| minimum angle' $O/${n}_stats.md | awk -F'|' '{print $4 "|" $6}')
    echo "$n | $tri | $q"
  done
done
# Extra runs quoted in the increment file, same flags as above unless named:
#   b<0|50>_tol<20|10>_a15                     --start-min-angle 15
#   b50_tol50_a25_g<2|5|10>                    --start-quality-gain 2, 5, 10
#   nolc_tol<50|20|10>_a<0|15|25>              no --features at all
# Table of every run: summarise with
#   for f in $O/*_stats.md; do ...; done   (see meshes.txt for the columns)
