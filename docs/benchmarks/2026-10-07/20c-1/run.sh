#!/bin/bash
# 20c-1 acceptance, both catchments (@perf, 2026-10-07).
# usage: run.sh R1 R2 ...   (repeat numbers; each repeat runs all eight configurations: base and head, each with and
# without --no-constraint-feet (variants base, basenofeet, head, nofeet),
# in an order that reverses on even repeats so drift spreads over all of them)
# base = master ed125121 (scratch worktree), head = worktree-soft-quality 34b546c4.
# Both _core builds are tools/bench.py's: Release, hardening on, <tree>/build-bench/pkg.
# Stops (exit 3) when not on AC before a run; a run whose power changed is renamed *.DISCARD.
HERE=$(cd "$(dirname "$0")" && pwd)
HEAD_TREE=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c1
BASE_TREE=$S/base
PY=$HEAD_TREE/.venv/bin/python
D=/Users/skavhaug/projects/rasputin_data
N=/Users/skavhaug/projects/rasputin_scratch/norway/numedalslagen
RAW=$HERE/raw; MESH=$S/meshes; mkdir -p "$RAW" "$MESH"
LAGAN=(--dem glo30 --cache $D/sweden_glo30_cache --out-crs EPSG:3006
  --domain $D/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson
  --features $D/sweden_corine/lagan_clc2018_3035.geojson --features-crs EPSG:3035
  --features-map corine --tolerance 10 --binary)
NUMED=(--dem $D/DTM10_UTM33_20260925 --domain $N/numedalslagen_outline_nve.geojson
  --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018
  --features-map corine --tolerance 10 --binary)
one() {  # one CATCHMENT VARIANT REPEAT
  local c=$1 v=$2 r=$3 name=$1-$2-r$3 tree=$HEAD_TREE extra=() args
  [ $v = base ] || [ $v = basenofeet ] && tree=$BASE_TREE
  [ $v = nofeet ] || [ $v = basenofeet ] && extra=(--no-constraint-feet)
  if [ $c = lagan ]; then args=("${LAGAN[@]}"); else args=("${NUMED[@]}"); fi
  local before; before=$(pmset -g batt)
  if ! grep -q "AC Power" <<<"$before"; then echo "$(date -u +%FT%TZ) STOP before $name: off AC"; exit 3; fi
  echo "$(date -u +%FT%TZ) start $name $(tail -1 <<<"$before")"
  local t0; t0=$(python3 -c 'import time; print(time.time())')
  (cd "$MESH" && "$PY" "$HERE/drive.py" --pkg "$tree/build-bench/pkg" --json "$RAW/$name.json" -- \
     mesh "${args[@]}" "${extra[@]}" --out "$MESH/$c-$v.vtk" --stats "$RAW/${name}_stats.md") \
     > "$RAW/$name.log" 2>&1
  local rc=$?
  local after; after=$(pmset -g batt)
  local t1; t1=$(python3 -c 'import time; print(time.time())')
  printf 'before:\n%s\nafter:\n%s\nexit: %s\nproc_s: %s\n' "$before" "$after" "$rc" \
    "$(python3 -c "print(round($t1-$t0,3))")" > "$RAW/$name.power"
  echo "$(date -u +%FT%TZ) end $name exit $rc $(tail -1 <<<"$after")"
  if [ $rc -ne 0 ]; then echo "FAILED $name"; exit 4; fi
  if ! grep -q "AC Power" <<<"$after"; then
    for f in "$RAW/$name".*  "$RAW/${name}_stats.md"; do mv "$f" "$f.DISCARD"; done
    echo "$(date -u +%FT%TZ) DISCARD $name: power changed during the run"; exit 3
  fi
}
for r in "$@"; do
  order=(lagan:base lagan:head lagan:nofeet lagan:basenofeet numed:base numed:head numed:nofeet numed:basenofeet)
  if [ $((r % 2)) -eq 0 ]; then order=(numed:basenofeet numed:nofeet numed:head numed:base lagan:basenofeet lagan:nofeet lagan:head lagan:base); fi
  for spec in "${order[@]}"; do
    [ -e "$RAW/${spec%%:*}-${spec##*:}-r$r.json" ] && continue  # done already (a rerun after a stop)
    one "${spec%%:*}" "${spec##*:}" "$r"
  done
done
echo DONE
