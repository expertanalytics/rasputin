#!/bin/bash
# Increment 16c acceptance on Bygdin (@perf, 2026-09-29): the reduced
# catchment meshed with the CORINE extract, at 10 m and 1 m. Raw logs go to
# ./logs; meshes go to $SCRATCH and are not committed.
# Usage: SCRATCH=/some/dir bash run.sh   (from the repository root, .venv active)
set -u
HERE=docs/benchmarks/2026-09-29/bygdin-landcover
LOG=$HERE/logs
DEM=../rasputin_data/DTM10_UTM33_20260925
GPKG=../rasputin_data/corine2018_dtm10_utm33.gpkg
DOM=$HERE/bygdin_reduced_t20.geojson
mkdir -p "$LOG" "$SCRATCH"
: > "$LOG/exits.txt"

mesh() {  # $1 tolerance, $2 tag, $3 ascii|binary
  /usr/bin/time -l rasputin mesh --dem $DEM --domain $DOM --tolerance "$1" \
    --features $GPKG --features-layer corine2018 --features-map corine \
    --"$3" --out "$SCRATCH/$2.vtk" --stats "$LOG/$2.stats.md" \
    > "$LOG/$2.out" 2> "$LOG/$2.err"
  echo "$2 exit=$?" | tee -a "$LOG/exits.txt"
}

git rev-parse HEAD > "$LOG/commit.txt"
shasum -a 256 .venv/lib/python3.*/site-packages/tin_engine/_core*.so \
  build-pyext/_core*.so >> "$LOG/commit.txt"
pmset -g batt > "$LOG/pmset_start.txt"
# Three timed binary runs at each tolerance; the median is reported.
for i in 1 2 3; do
  mesh 10 "lc_t10_run$i" binary
  mesh 1 "lc_t1_run$i" binary
done
pmset -g batt > "$LOG/pmset_timed_end.txt"
# One ASCII run each, for tools/bench.py's quality() (it reads ASCII only).
mesh 10 lc_t10_ascii ascii
mesh 1 lc_t1_ascii ascii
pmset -g batt > "$LOG/pmset_end.txt"
rasputin palette corine --out "$HERE/corine_natural.json" > "$LOG/palette.out" 2>&1
echo "palette exit=$?" | tee -a "$LOG/exits.txt"
