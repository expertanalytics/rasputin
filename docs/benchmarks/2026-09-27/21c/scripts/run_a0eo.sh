#!/bin/bash
# run_a0eo.sh PKG OUTDIR -- "A0 with evaluate-once": today and A0 once per domain,
# the three dirty rules A0eot / A0eo / A0eofp three times each (timers), a
# verification pass of each (RASPUTIN_SIM_EOVERIFY: every re-used footprint is
# recomputed and compared), and the evaluate-once plants on the quarter circle
# (on A0eot). A third argument `plants` runs only the plants. PKG: a build with sim.patch, sim_a0.patch and sim_a0_eo.patch applied.
set -e
P=$1; O=$2
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
D=$R/docs/benchmarks/2026-09-27/21c/scripts/sim_driver.py
Q=$R/docs/benchmarks/2026-09-26/quarter.geojson
PY=$R/.venv/bin/python
mkdir -p "$O"
pmset -g batt > "$O/pmset_before${3:+_$3}.txt"
[ "$3" = plants ] || for dom in "$Q" tile; do
  for m in today A0; do $PY "$D" "$P" $m "$dom" "$O"; done
  for rep in 1 2 3; do
    for m in A0eot A0eo A0eofp; do $PY "$D" "$P" $m "$dom" "$O/rep$rep"; done
  done
  for m in A0eot A0eo A0eofp; do RASPUTIN_SIM_EOVERIFY=1 $PY "$D" "$P" $m "$dom" "$O/verify"; done
done
# eo_stalecommit is not run here: it crashes the process (SIGSEGV/SIGBUS, a macOS
# crash dialog each time). Run it only under lldb with the hardened build, as in
# data/a0eo/sanitised/README.txt.
for pl in eo_nodirty eo_ringless eo_nofooted; do
  mkdir -p "$O/plant_$pl"
  RASPUTIN_SIM_PLANT=$pl RASPUTIN_SIM_EOVERIFY=1 perl -e 'alarm shift; exec @ARGV' 600 $PY "$D" "$P" A0eot "$Q" "$O/plant_$pl" \
    > "$O/plant_$pl/stdout.txt" 2>&1 || echo "plant $pl: exit $?" | tee -a "$O/plant_$pl/stdout.txt"
  cat "$O/plant_$pl/stdout.txt" | grep -v Warning
done
pmset -g batt > "$O/pmset_after${3:+_$3}.txt"
