#!/bin/bash
# usage: stats.sh <side: base|branch> <catchment> <rep>
# base = a7154ec's non-editable install (worktree-bottlenecks/.venv); branch = this worktree's .venv.
# Meshes go to the scratchpad S (not committed); stats and power state go to raw/stats/.
WT=/Users/skavhaug/projects/rasputin/.claude/worktrees
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/30a
R=$WT/landcover-speed/docs/benchmarks/2026-10-06/30a-landcover/raw/stats
side=$1; c=$2; r=$3
case $side in base) V=$WT/bottlenecks/.venv;; branch) V=$WT/landcover-speed/.venv;; esac
mkdir -p $S $R
o=${side}_${c}_r${r}
pmset -g batt | head -1 > $R/${o}_power.txt
$V/bin/rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 \
 --domain /Users/skavhaug/projects/rasputin_scratch/norway/$c/${c}_outline_nve.geojson \
 --features /Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine \
 --tolerance 10 --binary --out $S/${o}.vtk --stats $R/${o}_stats.md 2> $R/${o}_stderr.txt
echo "exit $?" >> $R/${o}_stderr.txt
pmset -g batt | head -1 >> $R/${o}_power.txt
shasum -a 256 $S/${o}.vtk | awk -v o=$o '{print o, $1}' >> $R/vtk_sha256.txt
