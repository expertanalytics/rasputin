#!/bin/bash
# usage: stats.sh <side: base|branch> <catchment> <rep>
# base   = 6c729e97's non-editable install, scratch venv S/venv-base   (uv pip install ".[codecs]" from `git archive 6c729e97`).
# branch = 0663efc9's non-editable install, scratch venv S/venv-branch (uv pip install ".[codecs]" from `git archive 0663efc9`).
# Each run calls its venv's own `python` (never a bin/ entry script, which can name another venv's
# interpreter: docs/increments/30c-dem-read-speed.md, section 8) and prints the tin_engine it imported.
# Power: `pmset -g batt` before and after, and every "Using Batt" line `pmset -g log` holds from the
# run's start (a change of power source inside the run). Exit 3 when the run was not on AC throughout.
# Meshes go to the scratchpad S (not committed); stats, stderr and power go to raw/stats/.
WT=/Users/skavhaug/projects/rasputin/.claude/worktrees/dem-read-speed
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/30c
R=$WT/docs/benchmarks/2026-10-06/30c-dem-read/raw/stats
D=/Users/skavhaug/projects/rasputin_data
side=$1; c=$2; r=$3
case $side in base|branch) P=$S/venv-$side/bin/python;; *) echo "side?"; exit 2;; esac
mkdir -p $S/out $R
o=${side}_${c}_r${r}
start=$(date '+%Y-%m-%d %H:%M:%S')
echo "before $start: $(pmset -g batt | tr '\n\t' '  ')" > $R/${o}_power.txt
$P -c 'import sys, tin_engine; print("tin_engine from", tin_engine.__file__, file=sys.stderr, flush=True); sys.argv[0] = "rasputin"; from tin_engine.cli import app; app()' \
 mesh --dem $D/DTM10_UTM33_20260925 \
 --domain /Users/skavhaug/projects/rasputin_scratch/norway/$c/${c}_outline_nve.geojson \
 --features $D/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine \
 --tolerance 10 --binary --out $S/out/${o}.vtk --stats $R/${o}_stats.md 2> $R/${o}_stderr.txt
echo "exit $?" >> $R/${o}_stderr.txt
echo "after $(date '+%Y-%m-%d %H:%M:%S'): $(pmset -g batt | tr '\n\t' '  ')" >> $R/${o}_power.txt
batt=$(pmset -g log | awk -v s="$start" 'substr($0, 1, 19) >= s && tolower($0) ~ /using batt/' | wc -l | tr -d ' ')
echo "pmset log lines 'Using Batt' since start: $batt" >> $R/${o}_power.txt
shasum -a 256 $S/out/${o}.vtk | awk -v o=$o '{print o, $1}' >> $R/vtk_sha256.txt
if [ "$(grep -c "AC Power" $R/${o}_power.txt)" != 2 ] || [ "$batt" != 0 ]; then exit 3; fi
