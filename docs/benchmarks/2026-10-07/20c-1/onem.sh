#!/bin/bash
# The 1 m tile's --stats phases, base vs 20c-1 (feet on and off), 7 repeats interleaved,
# through drive.py (same builds as run.sh). Default threads. Stops off AC.
HERE=$(cd "$(dirname "$0")" && pwd)
HEAD_TREE=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c1
PY=$HEAD_TREE/.venv/bin/python; RAW=$HERE/raw1m; mkdir -p "$RAW" "$S/meshes"
DEM=$HEAD_TREE/tests/fixtures/dem_archive/7908_3_10m_z33.tif
for r in 1 2 3 4 5 6 7; do
  for v in base head nofeet basenofeet; do
    tree=$HEAD_TREE; extra=()
    [ $v = base ] || [ $v = basenofeet ] && tree=$S/base
    [ $v = nofeet ] || [ $v = basenofeet ] && extra=(--no-constraint-feet)
    before=$(pmset -g batt); grep -q "AC Power" <<<"$before" || { echo "STOP off AC"; exit 3; }
    "$PY" "$HERE/drive.py" --pkg "$tree/build-bench/pkg" --json "$RAW/tile-$v-r$r.json" -- mesh --dem "$DEM" \
      --tolerance 1 --binary "${extra[@]}" --out "$S/meshes/tile-$v.vtk" --stats "$RAW/tile-$v-r${r}_stats.md" > "$RAW/tile-$v-r$r.log" 2>&1 || exit 4
    after=$(pmset -g batt)
    printf 'before:\n%s\nafter:\n%s\nexit: 0\nproc_s: 0\n' "$before" "$after" > "$RAW/tile-$v-r$r.power"
    grep -q "AC Power" <<<"$after" || { echo "DISCARD tile-$v-r$r"; exit 3; }
  done
done
echo DONE
