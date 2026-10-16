#!/bin/zsh
D=/Users/skavhaug/projects/rasputin_data; S=/Users/skavhaug/projects/rasputin_scratch/perf-20c3-20261008
LAG=(--dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006 --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --features $D/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine --tolerance 10)
NUM=(--dem $D/DTM10_UTM33_20260925 --domain $D/numedalslagen_outline_nve.geojson --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --tolerance 10)
B=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/.venv/bin/rasputin
M=/Users/skavhaug/projects/rasputin/.venv/bin/rasputin
for t in branch master; do
  [[ $t == branch ]] && X=$B || X=$M
  for c in lagan num; do
    [[ $c == lagan ]] && A=($LAG) || A=($NUM)
    /usr/bin/time -p $X mesh $A --out $S/$t-$c.vtk --stats $S/$t-$c.md > $S/$t-$c.log 2>&1
    echo "$t $c exit $?" >> $S/done.txt
  done
done
echo ALLDONE >> $S/done.txt
