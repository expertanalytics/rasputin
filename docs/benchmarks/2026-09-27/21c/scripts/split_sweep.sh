#!/bin/bash
# split_sweep.sh PKG OUT -- the split phase at 1 and 8 threads (section 6, Q5; THREADS="4 8 16" overrides):
# 3 passes, each a process of 5 refine calls per (domain, threads), interleaved.
# PKG: build-bench/pkg of the unpatched tree (bench.py's build()).
set -e
P=$1; O=$2
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
DRV=$R/docs/benchmarks/2026-09-27/serial-profile/scripts/prof_driver.py
DEM=$R/tests/fixtures/dem_archive/7908_3_10m_z33.tif
Q=$R/docs/benchmarks/2026-09-26/quarter.geojson
: > "$O"
pmset -g batt >> "$O.pmset"
for pass in 1 2 3; do
  for dom in quarter tile; do
    for t in ${THREADS:-1 8}; do
      extra=(); [ $dom = quarter ] && extra=(--domain "$Q")
      "$R/.venv/bin/python" "$DRV" --pkg "$P" --threads $t --repeat 5 -- mesh --dem "$DEM" \
        --tolerance 1 "${extra[@]}" --out "${TMPDIR:-/tmp}/rasputin-21c-split.vtk" --binary 2>&1 \
        | grep '^PHASES' | sed "s/^PHASES /$dom /" >> "$O"
    done
  done
done
pmset -g batt >> "$O.pmset"
