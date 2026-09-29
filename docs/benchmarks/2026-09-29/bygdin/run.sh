#!/bin/bash
# Increment 22 acceptance on Bygdin (@perf, 2026-09-29). Runs every timed
# block and writes raw logs to ./logs; meshes and GeoJSONs go to $SCRATCH.
# Usage: SCRATCH=/some/dir bash run.sh   (from the repository root, .venv active)
set -u
HERE=docs/benchmarks/2026-09-29/bygdin
LOG=$HERE/logs
DEM=../rasputin_data/DTM10_UTM33_20260925
GPKG=../rasputin_data/corine2018_dtm10_utm33.gpkg
mkdir -p "$LOG" "$SCRATCH"

catchment() {  # $1 tolerance ("" = default), $2 tag
  local tol=()
  [ -n "$1" ] && tol=(--outline-tolerance "$1")
  /usr/bin/time -l rasputin catchment --dem $DEM --seed 8.5425 61.3512 \
    --lakes $GPKG --lakes-layer corine2018 ${tol[@]+"${tol[@]}"} \
    --out "$SCRATCH/$2.geojson" > "$LOG/$2.out" 2> "$LOG/$2.err"
  echo "$2 exit=$?" | tee -a "$LOG/exits.txt"
}

mesh() {  # $1 domain geojson, $2 tolerance, $3 tag, $4 ascii|binary, rest: extra args
  local dom=$1 tol=$2 tag=$3 fmt=$4; shift 4
  /usr/bin/time -l rasputin mesh --dem $DEM --domain "$dom" --tolerance "$tol" \
    "$@" --"$fmt" --out "$SCRATCH/$tag.vtk" --stats "$LOG/$tag.stats.md" \
    > "$LOG/$tag.out" 2> "$LOG/$tag.err"
  echo "$tag exit=$?" | tee -a "$LOG/exits.txt"
}

# Block 1: the catchment.
pmset -g batt > "$LOG/pmset_catchment_start.txt"
for i in 1 2 3; do catchment "" "catchment_t20_run$i"; done
for t in 0 10 20 50; do catchment "$t" "catchment_t$t"; done
pmset -g batt > "$LOG/pmset_catchment_end.txt"
cp "$SCRATCH/catchment_t20_run1.geojson" "$HERE/bygdin_reduced_t20.geojson"

# Block 2: meshes. Three timed binary runs each, then one ASCII run for the
# quality check (tools/bench.py's quality(), which reads ASCII only).
RED=$SCRATCH/catchment_t20_run1.geojson
FINE=$SCRATCH/catchment_t0.geojson
pmset -g batt > "$LOG/pmset_mesh_start.txt"
for i in 1 2 3; do
  mesh "$RED" 10 "mesh_red_t10_run$i" binary
  mesh "$RED" 1 "mesh_red_t1_run$i" binary
  mesh "$FINE" 10 "mesh_fine_t10_run$i" binary
  mesh "$FINE" 1 "mesh_fine_t1_run$i" binary
  mesh "$RED" 10 "mesh_feat_t10_run$i" binary \
    --features $GPKG --features-layer corine2018 --features-map corine
done
pmset -g batt > "$LOG/pmset_mesh_end.txt"
mesh "$RED" 10 mesh_red_t10_ascii ascii
mesh "$RED" 1 mesh_red_t1_ascii ascii
mesh "$FINE" 10 mesh_fine_t10_ascii ascii
mesh "$FINE" 1 mesh_fine_t1_ascii ascii
mesh "$RED" 10 mesh_feat_t10_ascii ascii \
  --features $GPKG --features-layer corine2018 --features-map corine
pmset -g batt > "$LOG/pmset_quality_end.txt"
