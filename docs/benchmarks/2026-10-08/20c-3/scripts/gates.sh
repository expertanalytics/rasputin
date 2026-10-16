#!/bin/zsh
S=/Users/skavhaug/projects/rasputin_scratch/perf-20c3-20261008; D=/Users/skavhaug/projects/rasputin_data
LAG=(--dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006 --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --features $D/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine --tolerance 10)
NUM=(--dem $D/DTM10_UTM33_20260925 --domain $D/numedalslagen_outline_nve.geojson --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --tolerance 10)
B=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3/.venv/bin/rasputin
for v in "tol2:0:0" "tol2snap5:5:0" "tol2rep:0:0.05" "tol2snap5rep:5:0.05"; do
  n=${v%%:*}; r=${v#*:}; snap=${r%%:*}; rep=${r#*:}
  for c in lagan num; do
    [[ $c == lagan ]] && A=($LAG) || A=($NUM)
    $B mesh $A --no-features-merge-same-class --features-tolerance 2 --features-outline-snap $snap --features-repair $rep --out $S/g-$n-$c.vtk --stats $S/g-$n-$c.md > $S/g-$n-$c.log 2>&1
    echo "$n $c exit $?" >> $S/gdone.txt
  done
done
echo ALLDONE >> $S/gdone.txt
