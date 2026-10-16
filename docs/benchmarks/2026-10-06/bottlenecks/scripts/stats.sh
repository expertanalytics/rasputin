#!/bin/bash
# usage: stats.sh <catchment> <rep>
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/bottlenecks
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/bn
c=$1; r=$2
pmset -g batt | head -1 > $S/${c}_r${r}_power.txt
$W/.venv/bin/rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 \
 --domain /Users/skavhaug/projects/rasputin_scratch/norway/$c/${c}_outline_nve.geojson \
 --features /Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine \
 --tolerance 10 --binary --out $S/${c}_r${r}.vtk --stats $S/${c}_r${r}_stats.md
pmset -g batt | head -1 >> $S/${c}_r${r}_power.txt
