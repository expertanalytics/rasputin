#!/bin/zsh
# R3 (d) and (e) through r3drive.py; meshes in scratch.
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-3; P=$W/docs/benchmarks/2026-10-08/20c-3/scripts/r3drive.py
D=/Users/skavhaug/projects/rasputin_data; S=/Users/skavhaug/projects/rasputin_scratch/perf-20c3-r3-20261008; mkdir -p $S
LAG=(--dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006 --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --features $D/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine --tolerance 10)
NUM=(--dem $D/DTM10_UTM33_20260925 --domain $D/numedalslagen_outline_nve.geojson --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --tolerance 10)
G=(--no-features-merge-same-class --features-tolerance 2 --features-outline-snap 5 --features-repair 0.05)
PY=$W/.venv/bin/python
$PY $P d-lagan EPSG:3006 -- $LAG $G --out $S/d-lagan.vtk > $S/d-lagan.txt 2>&1
$PY $P d-num EPSG:25833 -- $NUM $G --out $S/d-num.vtk > $S/d-num.txt 2>&1
$PY $P e-lagan-merge EPSG:3006 -- $LAG --out $S/e-lagan-merge.vtk > $S/e-merge.txt 2>&1
$PY $P e-lagan-nomerge EPSG:3006 -- $LAG --no-features-merge-same-class --out $S/e-lagan-nomerge.vtk > $S/e-nomerge.txt 2>&1
echo DEDONE >> $S/done.txt
