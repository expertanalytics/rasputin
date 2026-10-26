#!/bin/sh
# 106a's land-cover repair, master against the head, one run each, back to back.
# Usage: sh run.sh CASE   (CASE is kultsjon or sjodal)
# Needs both Release packages built first (bench.py's build(): <tree>/build-bench/pkg).
# Writes OUT/<case>-<side>.{gpkg,stats.json,log}; --stats writes Markdown whatever
# the name. Meshes stay in OUT. sjodal needs OUT/nedre-sjodalsvatn.geojson: the copy
# beside this script (NVE gauge 2.13.0's outline, EPSG:25833, one Polygon).
set -u
REPO=/Users/skavhaug/projects/rasputin/.claude/worktrees/demo-fixes
PY=$REPO/.venv/bin/python
OUT=/Users/skavhaug/projects/rasputin_scratch/perf-106a-out
case "$1" in
  kultsjon) ARGS="--dem glo30 --out-crs EPSG:3006 --domain /Users/skavhaug/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson --land-cover clc2018 --triangles 6000" ;;
  sjodal) ARGS="--dem dtm10-norway --domain $OUT/nedre-sjodalsvatn.geojson --land-cover ar50-innlandet --reservoirs nve-reservoirs --triangles 6000" ;;
esac
for side in master head; do
  if [ $side = master ]; then PKG=/Users/skavhaug/projects/rasputin_scratch/perf-106a-base/build-bench/pkg; else PKG=$REPO/build-bench/pkg; fi
  pmset -g batt | head -1
  /usr/bin/time -p $PY -I $REPO/tools/bench.py _child --pkg $PKG --threads 0 -- mesh \
    --out $OUT/$1-$side.gpkg $ARGS --stats $OUT/$1-$side.stats.json > $OUT/$1-$side.log 2>&1
  echo "$1 $side exit $?"; tail -4 $OUT/$1-$side.log
done
