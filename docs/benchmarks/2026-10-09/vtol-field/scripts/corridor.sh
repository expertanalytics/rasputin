#!/bin/sh
# The Hokksund-Bergen corridor with the ramp, once, for the record (not judged).
# Usage, from the worktree root: corridor.sh RASPUTIN_DATA OUT_DIR EVIDENCE_DIR
set -eu
D=$1; O=$2; E=$3
pmset -g batt > "$E/corridor.power-before.txt"
/usr/bin/time -l .venv/bin/python tools/bench.py _child --pkg "$PWD/build-bench/pkg" --threads 0 -- \
  mesh --dem "$D/DTM10_UTM33_20260925" --domain "$D/banenor_banenettverk/corridor_domain.geojson" \
  --tolerance 20 --tolerance-near "$D/banenor_banenettverk/bergen_line3.geojson" 1 --tolerance-ramp 0 3000 \
  --binary --out "$O/corridor_ramp.vtk" --stats "$E/corridor.stats.md" 2> "$E/corridor.stderr.txt"
pmset -g batt > "$E/corridor.power-after.txt"
