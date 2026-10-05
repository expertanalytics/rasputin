#!/bin/bash
# Step 4's probes (@perf), after run.sh and checks.py; outputs in step4/.
# Usage, from the repository root:
#   RASPUTIN_DATA=<data root> WORK=<scratch dir> bash docs/benchmarks/2026-10-05/nve-hrd/step4.sh
# None of these is timed. The re-runs of the misses are rerun.sh's.
set -eu
HERE=docs/benchmarks/2026-10-05/nve-hrd
DATA="${RASPUTIN_DATA:?}/nve_hrd"
DEM="$RASPUTIN_DATA/DTM10_UTM33_20260925"
PY=".venv/bin/python"
S="$HERE/step4"
mkdir -p "$S" "${WORK:?}/probe"
quiet() { grep -v RuntimeWarning || true; }

bash "$HERE/rerun.sh" 12.215.0 15.49.0 83.2.0 105.1.0 153.1.0 307.7.0
$PY "$HERE/explain.py" "$DATA" "$DEM" "$HERE" 3000 12.215.0 15.49.0 83.2.0 105.1.0 153.1.0 307.7.0 \
  26.29.0 35.9.0 16.66.0 19.79.0 6.10.0 2>&1 | quiet > "$S/explain_3km.txt"
$PY "$HERE/explain.py" "$DATA" "$DEM" "$HERE" 12000 105.1.0 153.1.0 2>&1 | quiet > "$S/explain_12km.txt"
$PY "$HERE/divide.py" "$DATA" "$DEM" "$HERE" 12.215.0 15.49.0 307.7.0 26.29.0 35.9.0 16.66.0 \
  19.79.0 6.10.0 2>&1 | quiet > "$S/divide.txt"
: > "$S/chain_counts_by_window.txt"
for h in 2600 4700 6000 8000 10000 12000; do
  echo "half $h" >> "$S/chain_counts_by_window.txt"
  $PY "$HERE/chain_counts.py" "$DATA" "$DEM" $h 30 105.1.0 153.1.0 2>&1 | quiet | grep -E "==|placed" \
    >> "$S/chain_counts_by_window.txt"
done
$PY "$HERE/chain_counts.py" "$DATA" "$DEM" 6000 30 2.279.0 2>&1 | quiet > "$S/chain_counts_2.279.0.txt"
$PY "$HERE/trace_exit.py" "$DATA" "$DEM" 12000 4700 105.1.0 2>&1 | quiet > "$S/trace_exit_105.1.0.txt"
for sid in 2.142.0 105.1.0 212.48.0; do
  $PY "$HERE/probe_windows.py" "$DATA" "$DEM" "$WORK/probe/$sid" $sid 2>&1 | quiet > "$S/windows_$sid.txt"
done
$PY "$HERE/window_check.py" "$DATA" "$DEM" "$HERE" 12000 "$HERE/window_check_12km.csv" 2>&1 | quiet > "$S/window_check.log"
$PY "$HERE/bypass.py" "$DATA" "$DEM" "$HERE" 6000 "$HERE/bypass.csv" 2>&1 | quiet
$PY "$HERE/findings.py" "$HERE" 2>&1 | quiet > "$S/uncertain_table.md"
$PY "$HERE/lake_above.py" "$DATA" "$HERE" 2>&1 | quiet > "$S/lake_above.txt"
