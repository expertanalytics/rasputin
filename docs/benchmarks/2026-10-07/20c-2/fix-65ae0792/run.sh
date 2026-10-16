#!/bin/bash
# 20c-2 re-time after the speed fix 65ae0792 (@perf, 2026-10-07): the --stats phases and the exact
# mesh measures on Lagan, Numedalslagen and the 1 m tile. Copied from ../run.sh. Variants:
#   base   master 41bda81a (20c-1 merged), its default flags
#   head   20c-2 65ae0792, default flags (--start-quality-gain 0)
#   headm1 20c-2 65ae0792, --start-quality-gain -1 (must equal base bit for bit)
#   nofeet 20c-2 65ae0792, --no-constraint-feet (the split-phase limit's "feet off")
# usage: run.sh R1 R2 ...        every (catchment, variant), per repeat, in an order that
#                                reverses on even repeats so drift spreads over all of them
#        run.sh tile R1 R2 ...   the tile only, base, head and nofeet in turn (order reversed on even
#                                repeats), into raw-tile/: the longer series for the split-phase limit's point 1
# Every _core is tools/bench.py's build (pairs.sh builds them first): Release, hardening on,
# <tree>/build-bench/pkg. Stops (exit 3) when not on AC before a run; a run whose power changed
# is renamed *.DISCARD and the script stops.
HERE=$(cd "$(dirname "$0")" && pwd)
HEAD_TREE=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b
PY=$HEAD_TREE/.venv/bin/python
D=/Users/skavhaug/projects/rasputin_data
N=/Users/skavhaug/projects/rasputin_scratch/norway/numedalslagen
RAW=$HERE/raw; [ "$1" = tile ] && RAW=$HERE/raw-tile; MESH=$S/meshes; mkdir -p "$RAW" "$MESH"
LAGAN=(--dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006
  --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson
  --features $D/sweden_corine/lagan_clc2018_3035.geojson --features-crs EPSG:3035
  --features-map corine --tolerance 10 --binary)
NUMED=(--dem $D/DTM10_UTM33_20260925 --domain $N/numedalslagen_outline_nve.geojson
  --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018
  --features-map corine --tolerance 10 --binary)
TILE=(--dem $HEAD_TREE/tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --binary)
one() {  # one CATCHMENT VARIANT REPEAT
  local c=$1 v=$2 r=$3 name=$1-$2-r$3 tree=$HEAD_TREE extra=() args
  case $v in base) tree=$S/base;; headm1) extra=(--start-quality-gain -1);;
    nofeet) extra=(--no-constraint-feet);; esac
  case $c in lagan) args=("${LAGAN[@]}");; numed) args=("${NUMED[@]}");; tile) args=("${TILE[@]}");; esac
  local before; before=$(pmset -g batt)
  if ! grep -q "AC Power" <<<"$before"; then echo "$(date -u +%FT%TZ) STOP before $name: off AC"; exit 3; fi
  echo "$(date -u +%FT%TZ) start $name $(tail -1 <<<"$before")"
  local t0; t0=$(python3 -c 'import time; print(time.time())' 2>/dev/null)
  (cd "$MESH" && "$PY" "$HERE/drive.py" --pkg "$tree/build-bench/pkg" --json "$RAW/$name.json" -- \
     mesh "${args[@]}" "${extra[@]}" --out "$MESH/$c-$v.vtk" --stats "$RAW/${name}_stats.md") \
     > "$RAW/$name.log" 2>&1
  local rc=$?
  local after; after=$(pmset -g batt)
  local t1; t1=$(python3 -c 'import time; print(time.time())' 2>/dev/null)
  printf 'before:\n%s\nafter:\n%s\nexit: %s\nproc_s: %s\n' "$before" "$after" "$rc" \
    "$(python3 -c "print(round($t1-$t0,3))" 2>/dev/null)" > "$RAW/$name.power"
  echo "$(date -u +%FT%TZ) end $name exit $rc $(tail -1 <<<"$after")"
  if [ $rc -ne 0 ]; then echo "FAILED $name"; exit 4; fi
  if ! grep -q "AC Power" <<<"$after"; then
    for f in "$RAW/$name".*  "$RAW/${name}_stats.md"; do mv "$f" "$f.DISCARD"; done
    echo "$(date -u +%FT%TZ) DISCARD $name: power changed during the run"; exit 3
  fi
}
if [ "$1" = tile ]; then
  shift
  for r in "$@"; do
    order=(tile:base tile:head tile:nofeet); [ $((r % 2)) -eq 0 ] && order=(tile:nofeet tile:head tile:base)
    for spec in "${order[@]}"; do
      [ -e "$RAW/${spec%%:*}-${spec##*:}-r$r.json" ] || one "${spec%%:*}" "${spec##*:}" "$r"
    done
  done
  echo DONE; exit 0
fi
for r in "$@"; do
  order=(tile:base tile:head tile:headm1 tile:nofeet numed:base numed:head numed:headm1 numed:nofeet
         lagan:base lagan:head lagan:headm1 lagan:nofeet)
  if [ $((r % 2)) -eq 0 ]; then
    order=(lagan:nofeet lagan:headm1 lagan:head lagan:base numed:nofeet numed:headm1 numed:head numed:base
           tile:nofeet tile:headm1 tile:head tile:base)
  fi
  for spec in "${order[@]}"; do
    [ -e "$RAW/${spec%%:*}-${spec##*:}-r$r.json" ] && continue  # done already (a rerun after a stop)
    one "${spec%%:*}" "${spec##*:}" "$r"
  done
done
echo DONE
