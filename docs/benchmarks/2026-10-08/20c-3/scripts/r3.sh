#!/bin/zsh
# R3 (a) and (b) on aaee898b; meshes in scratch.
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3
D=/Users/skavhaug/projects/rasputin_data; S=/Users/skavhaug/projects/rasputin_scratch/perf-20c3-r3-20261008; mkdir -p $S
LAG=(--dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006 --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --features $D/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine --tolerance 10)
NUM=(--dem $D/DTM10_UTM33_20260925 --domain $D/numedalslagen_outline_nve.geojson --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --tolerance 10)
B=$W/.venv/bin/rasputin
for c in lagan num; do [[ $c == lagan ]] && A=($LAG) || A=($NUM)
  $B mesh $A --out $S/a-$c.vtk --stats $S/a-$c.md > $S/a-$c.log 2>&1; echo "a $c exit $?" >> $S/done.txt; done
for c in lagan num; do [[ $c == lagan ]] && A=($LAG) || A=($NUM)
  $B mesh $A --no-features-merge-same-class --features-tolerance 2 --features-outline-snap 5 --features-repair 0.05 --out $S/b-$c.vtk --stats $S/b-$c.md > $S/b-$c.log 2>&1; echo "b $c exit $?" >> $S/done.txt; done
echo ALLDONE >> $S/done.txt
