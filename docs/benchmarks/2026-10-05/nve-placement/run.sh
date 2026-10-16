#!/bin/bash
# Placement figures, before and after (@perf, 2026-10-05; increment 29, PR 2).
# Usage, from the repository root:
#   RASPUTIN_DATA=<data root> MPL=<dir with matplotlib> WORK=<scratch dir> bash run.sh
# RASPUTIN_DATA is ../rasputin_data beside the main checkout (absolute: a
# worktree sits deeper, and a relative path would land elsewhere).
# The venv's python is used, with $MPL on PYTHONPATH (matplotlib is in no
# dependency list, Ola's ruling 2026-10-04). Work files (full-run outlines)
# go to $WORK and are not committed; the small tables are copied here.
set -eu
HERE=docs/benchmarks/2026-10-05/nve-placement
DATA="${RASPUTIN_DATA:?set RASPUTIN_DATA, the data root}/nve_hrd"
DEM="$RASPUTIN_DATA/DTM10_UTM33_20260925"
PY=".venv/bin/python"
export PYTHONPATH="$MPL"
mkdir -p "$WORK"

# The data, once (network): stations, NVE's polygons, ELVIS river lines.
.venv/bin/rasputin fetch-stations nve-hrd --out-dir "$DATA"

git rev-parse HEAD > "$HERE/provenance.txt"
shasum -a 256 .venv/lib/python3.*/site-packages/tin_engine/_core*.so build-pyext/_core*.so \
  "$DATA"/*.geojson >> "$HERE/provenance.txt"
pmset -g batt >> "$HERE/provenance.txt"

# The cases, with catchment --rivers's placement (the nearest line).
$PY $HERE/render.py survey "$DATA" "$DEM" "$WORK/nearest"
$PY $HERE/render.py choose "$DATA" "$DEM" "$WORK/nearest" 2> "$HERE/choose.log"
$PY $HERE/render.py draw "$DATA" "$DEM" "$WORK/nearest" "$HERE" 2.11.0 73.27.0 24.9.0
cp "$WORK/nearest/survey.csv" "$WORK/nearest/cases.json" "$HERE/"

# The same with the batch's tiered placement (PR 4's), for comparison.
$PY $HERE/render.py survey "$DATA" "$DEM" "$WORK/tiers" --tiers
$PY $HERE/render.py choose "$DATA" "$DEM" "$WORK/tiers" --tiers 2> "$HERE/choose_tiers.log"
# Etna (12.70.0): the tiered placement picks another line than the nearest.
$PY $HERE/render.py draw "$DATA" "$DEM" "$WORK/tiers" "$HERE" 12.70.0 --tiers --extras-only
cp "$WORK/tiers/survey.csv" "$HERE/survey_tiers.csv"
cp "$WORK/tiers/cases.json" "$HERE/cases_tiers.json"
pmset -g batt >> "$HERE/provenance.txt"

# The drawn stations through the command itself (catchment --rivers), to
# check that render.py's runs are the command's.
: > "$HERE/cli_confirm.log"
for sid in 101.1.0 2.32.0 2.11.0 73.27.0 24.9.0; do
  xy=$($PY -c "import json,sys; f=[f for f in json.load(open('$DATA/stations.geojson'))['features'] if f['properties']['station']=='$sid'][0]; print(*f['geometry']['coordinates'][:2])")
  echo "== $sid" >> "$HERE/cli_confirm.log"
  .venv/bin/rasputin catchment --dem "$DEM" --seed $xy --seed-crs EPSG:25833 \
    --rivers "$DATA/rivers.geojson" --out "$WORK/cli_$sid.geojson" 2>> "$HERE/cli_confirm.log" >/dev/null \
    || echo "exit $?" >> "$HERE/cli_confirm.log"
done

# full_runs*.json: every full run, less the outline and the chain (kept in $WORK).
for mode in nearest tiers; do
  suffix=""; [ $mode = tiers ] && suffix=_tiers
  $PY -c "
import glob, json, sys
out = {}
for f in sorted(glob.glob(sys.argv[1] + '/runs/*.json')):
    r = json.load(open(f)); r.pop('fine', None); r['chain_nodes'] = len(r.pop('chain', None) or [])
    out[r['station']] = r
open(sys.argv[2], 'w').write(json.dumps(out, indent=1))" "$WORK/$mode" "$HERE/full_runs$suffix.json"
done
