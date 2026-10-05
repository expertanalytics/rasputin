#!/bin/bash
# Increment 29's acceptance run over the 140 NVE HRD stations (@perf, 2026-10-05).
# Steps 1-3 and 6 of the increment file's "Acceptance: every covered HRD station".
# Usage, from the repository root (a worktree is fine):
#   RASPUTIN_DATA=<data root> WORK=<scratch dir> bash docs/benchmarks/2026-10-05/nve-hrd/run.sh
# RASPUTIN_DATA is ../rasputin_data beside the main checkout, absolute.
# The venv is the tree's own, installed non-editable (Release, bounds checks on):
#   uv venv --python 3.13 .venv
#   uv pip install --python .venv --reinstall-package rasputin --no-cache ".[codecs]"
# Outputs: the batch writes to $WORK/batch; results.csv, summary.json, the log,
# the memory trace and the catchment files (ours, as catchments/<station>.geojson)
# are copied here. Then checks.py, then analyse.py again, which reads checks.csv.
# On 2026-10-05 the catchment copy, checks.py and the second analyse.py were run
# by hand after this script, with the same commands as its last lines below; the
# copy was checked against $WORK/batch with diff -rq (identical, 124 files).
# The lakes along each reach (lake_query.py, network) feed checks.py's
# lake_above_p column only; that query was made later, during step 4's review.
set -eu
HERE=docs/benchmarks/2026-10-05/nve-hrd
DATA="${RASPUTIN_DATA:?set RASPUTIN_DATA, the data root}/nve_hrd"
DEM="$RASPUTIN_DATA/DTM10_UTM33_20260925"
OUT="${WORK:?set WORK, a scratch directory}/batch"
mkdir -p "$WORK"

# Step 1: the station set (network). Files already there are kept unless --refresh.
.venv/bin/rasputin fetch-stations nve-hrd --out-dir "$DATA"

{ git rev-parse HEAD; .venv/bin/rasputin version; .venv/bin/python --version
  shasum -a 256 .venv/lib/python3.*/site-packages/tin_engine/_core*.so "$DATA"/*
  date; pmset -g batt; } > "$HERE/provenance.txt"
cp "$DATA/manifest.json" "$HERE/manifest.json"

# Step 2: the batch, every station, with --lakes, at the defaults.
# /usr/bin/time -l gives the total wall time and the process's peak RSS; each
# stderr line is time-stamped, and a sampler logs the RSS every 0.2 s, so a
# station's peak is the largest sample between its line and the previous one.
rm -rf "$OUT"
( /usr/bin/time -l .venv/bin/rasputin station-catchments --dem "$DEM" \
    --stations "$DATA/stations.geojson" --rivers "$DATA/rivers.geojson" \
    --reference "$DATA/reference.geojson" --lakes "$DATA/lakes.geojson" \
    --out-dir "$OUT" 2>&1 >/dev/null \
  | while IFS= read -r line; do printf '%s\t%s\n' "$(perl -MTime::HiRes=time -e 'printf "%.3f", time')" "$line"; done \
  > "$WORK/batch.log" ) &
RUN=$!
sleep 1
PID=$(pgrep -n -f "rasputin station-catchments" || true)
: > "$WORK/rss.tsv"
while kill -0 "$RUN" 2>/dev/null; do
  rss=$(ps -o rss= -p "$PID" 2>/dev/null || true)
  [ -n "$rss" ] && printf '%s\t%s\n' "$(perl -MTime::HiRes=time -e 'printf "%.3f", time')" "$rss" >> "$WORK/rss.tsv"
  sleep 0.2
done
wait "$RUN" || true
{ date; pmset -g batt; } >> "$HERE/provenance.txt"
cp "$OUT/results.csv" "$OUT/summary.json" "$WORK/batch.log" "$HERE/"
gzip -c "$WORK/rss.tsv" > "$HERE/rss.tsv.gz"

mkdir -p "$HERE/catchments"
cp "$OUT"/*.geojson "$HERE/catchments/"

# Steps 3 and 6, in this order: the tables, the checks, the tables again with
# the checks (analyse.py reads checks.csv when it is there).
.venv/bin/python "$HERE/analyse.py" "$HERE" > "$HERE/analysis.md"
.venv/bin/python "$HERE/lake_query.py" "$DATA" "$HERE" "$WORK/lakes_reach.geojson" > "$HERE/step4/lake_query.txt"
.venv/bin/python "$HERE/checks.py" "$DATA" "$DEM" "$OUT" "$HERE/checks.csv" "$WORK/lakes_reach.geojson"
.venv/bin/python "$HERE/analyse.py" "$HERE" > "$HERE/analysis.md"
