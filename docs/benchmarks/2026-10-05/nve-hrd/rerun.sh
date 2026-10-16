#!/bin/bash
# Step 4's re-runs (@perf): the stations named, at corridor 15 m and 60 m (map
# radius 500 m) and at map radius 250 m and 1000 m (corridor 30 m), through
# rerun.py. Same environment as run.sh. Each run's results.csv is copied here
# as reruns/<corridor>_<radius>.csv.
# Usage, from the repository root:
#   RASPUTIN_DATA=<data root> WORK=<scratch dir> bash rerun.sh STATION [STATION ...]
set -eu
HERE=docs/benchmarks/2026-10-05/nve-hrd
DATA="${RASPUTIN_DATA:?}/nve_hrd"
DEM="$RASPUTIN_DATA/DTM10_UTM33_20260925"
mkdir -p "$HERE/reruns"
for cfg in "15 500" "60 500" "30 250" "30 1000"; do
  read -r corridor radius <<< "$cfg"
  out="${WORK:?}/rerun_${corridor}_${radius}"
  rm -rf "$out"
  .venv/bin/python "$HERE/rerun.py" "$corridor" "$radius" "$DATA" "$DEM" "$out" "$@" \
    2> "$WORK/rerun_${corridor}_${radius}.log" >/dev/null
  cp "$out/results.csv" "$HERE/reruns/${corridor}_${radius}.csv"
  cp "$WORK/rerun_${corridor}_${radius}.log" "$HERE/reruns/${corridor}_${radius}.log"
done
