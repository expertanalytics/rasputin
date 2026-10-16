#!/bin/zsh
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3
D=/Users/skavhaug/projects/rasputin_data; S=/Users/skavhaug/projects/rasputin_scratch/perf-20c3-t12-20261008
LAG=(--dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006 --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --features $D/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine --tolerance 10)
NUM=(--dem $D/DTM10_UTM33_20260925 --domain $D/numedalslagen_outline_nve.geojson --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --tolerance 10)
B=$W/.venv/bin/rasputin
cd $W
git checkout 35b7a873 -- src_python/tin_engine/feature_input.py
for c in num lagan; do [[ $c == lagan ]] && A=($LAG) || A=($NUM)
  /usr/bin/time -p $B mesh $A --out $S/t1-$c.vtk --stats $S/t1-$c.md > $S/t1-$c.log 2>&1; echo "t1 $c exit $?" >> $S/done.txt; done
git checkout HEAD -- src_python/tin_engine/feature_input.py
git status --short >> $S/done.txt
for c in num lagan; do [[ $c == lagan ]] && A=($LAG) || A=($NUM)
  /usr/bin/time -p $B mesh $A --out $S/t2-$c.vtk --stats $S/t2-$c.md > $S/t2-$c.log 2>&1; echo "t2 $c exit $?" >> $S/done.txt; done
for c in num lagan; do [[ $c == lagan ]] && A=($LAG) || A=($NUM)
  $B mesh $A --no-features-merge-same-class --features-tolerance 2 --features-outline-snap 5 --features-repair 0 --out $S/g-tol2snap5-$c.vtk --stats $S/g-tol2snap5-$c.md > $S/g-tol2snap5-$c.log 2>&1; echo "gate $c exit $?" >> $S/done.txt; done
echo ALLDONE >> $S/done.txt
