#!/bin/bash
# run_sims.sh PKG OUTDIR -- every simulation of the README, both domains.
# PKG: build-bench/pkg of a worktree of 7f688aa with scripts/sim.patch applied.
set -e
P=$1; O=$2
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
D=$R/docs/benchmarks/2026-09-27/21c/scripts/sim_driver.py
Q=$R/docs/benchmarks/2026-09-26/quarter.geojson
for dom in "$Q" tile; do
  "$R/.venv/bin/python" "$D" "$P" today "$dom" "$O" --foot
  for m in C Cthin Cmis Cthin2 A1; do "$R/.venv/bin/python" "$D" "$P" $m "$dom" "$O"; done
  for s in 1 2 3; do "$R/.venv/bin/python" "$D" "$P" A1 "$dom" "$O" --seed $s; done
done
