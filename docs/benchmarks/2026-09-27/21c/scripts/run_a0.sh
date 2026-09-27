#!/bin/bash
# run_a0.sh PKG OUTDIR -- A0 section: today, A0, A1 (seed 0), C on both domains; A0 plants on the quarter.
set -e
P=$1; O=$2
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
D=$R/docs/benchmarks/2026-09-27/21c/scripts/sim_driver.py
Q=$R/docs/benchmarks/2026-09-26/quarter.geojson
pmset -g batt > "$O.pmset_before"
for dom in "$Q" tile; do
  for m in today A0 A1 C; do "$R/.venv/bin/python" "$D" "$P" $m "$dom" "$O"; done
done
for pl in norenumber noreserve nocavity; do
  mkdir -p "$O/plant_$pl"
  RASPUTIN_SIM_PLANT=$pl "$R/.venv/bin/python" "$D" "$P" A0 "$Q" "$O/plant_$pl"
done
pmset -g batt > "$O.pmset_after"
