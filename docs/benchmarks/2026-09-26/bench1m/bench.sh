#!/bin/bash
# bench.sh DOMAIN(quarter|tile) label worktree [extra args...]
# 3 timed --binary runs, then one --ascii run kept as runs/<dom>_<label>.vtk
B=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m
R=/Users/skavhaug/projects/rasputin
PY=$R/.venv/bin/python
DEM=$R/tests/fixtures/dem_archive/7908_3_10m_z33.tif
Q=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/paraview/inputs/quarter.geojson
dom=$1; l=$2; wt=$B/wt/$3; shift 3
D=""; [ $dom = quarter ] && D="--domain $Q"
S=""; [ -d $wt/src_python/tin_engine ] && grep -q -- '"--stats"' $wt/src_python/tin_engine/cli.py && S=1
mkdir -p $B/logs $B/runs
cd $R
for i in 1 2 3; do
  t0=$(perl -MTime::HiRes=time -e 'printf "%.4f", time')
  $PY $B/run.py $wt mesh --dem $DEM $D --tolerance 1 --binary --out $B/logs/tmp.vtk "$@" > $B/logs/${dom}_$l.t$i.out 2>&1
  t1=$(perl -MTime::HiRes=time -e 'printf "%.4f", time')
  echo "BENCH proc_s $(echo "$t1 - $t0" | bc)" >> $B/logs/${dom}_$l.t$i.out
done
ST=""; [ -n "$S" ] && ST="--stats $B/runs/${dom}_$l.stats.md"
$PY $B/run.py $wt mesh --dem $DEM $D --tolerance 1 --ascii --out $B/runs/${dom}_$l.vtk $ST "$@" > $B/logs/${dom}_$l.ascii.out 2>&1
$PY $B/quality.py $B/runs/${dom}_$l.vtk > $B/logs/${dom}_$l.q.json
echo "$dom $l done"
