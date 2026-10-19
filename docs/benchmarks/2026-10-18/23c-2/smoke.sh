#!/bin/zsh
# 23c-2 smoke on real data (point 8): the Velhas piece on ANADEM from the cache,
# EPSG:31983, 10 m, --pieces 64, run from this branch's Release package
# (build-bench/pkg, built by pairs.sh's head runs), then the checks.
# Usage: smoke.sh W S QPY QPKG
#   W: the 23c worktree; S: scratch directory for the piece files;
#   QPY: a Python with matplotlib (the checks); QPKG: a tin_engine package for QPY,
#   used by check_pieces.py for tin_engine.stats.quality only.
set -u
W=$1 S=$2 QPY=$3 QPKG=$4
O=$W/docs/benchmarks/2026-10-18/23c-2
DATA=/Users/skavhaug/projects/rasputin_data/sao_francisco_piece
mkdir -p $S/velhas
rm -rf $S/velhas/velhas.vtk.pieces
pmset -g batt | grep -q "AC Power" || { echo "NOT ON AC"; exit 3; }
pmset -g batt > $O/velhas.pmset; sysctl vm.swapusage >> $O/velhas.pmset
MallocLargeCache=0 /usr/bin/time -l $W/.venv/bin/python $O/smoke_driver.py $W/build-bench/pkg \
  $O/velhas_refine_calls.jsonl -- mesh --dem anadem-v1 --cache $HOME/.cache/rasputin \
  --domain $DATA/bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983 \
  --tolerance 10 --pieces 64 --ascii --out $S/velhas/velhas.vtk > $O/velhas.out 2> $O/velhas.err
echo "mesh exit $?"
pmset -g batt >> $O/velhas.pmset; sysctl vm.swapusage >> $O/velhas.pmset
QUALITY_PKG=$QPKG $QPY $O/check_pieces.py $S/velhas/velhas.vtk.pieces 10 $S/velhas/merged.vtk \
  --shift 30 > $O/check_pieces.json 2> $O/check_pieces.err
echo "check exit $?"
$QPY $O/once_on_outline.py $S/velhas/velhas.vtk.pieces $DATA/bho2017_5k_76949_outline_epsg4674.geojson \
  > $O/once_on_outline.txt 2>&1
echo "outline exit $?"
echo SMOKEDONE
