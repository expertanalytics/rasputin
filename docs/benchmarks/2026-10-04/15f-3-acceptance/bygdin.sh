#!/bin/zsh
# 15f-3 acceptance (2): Bygdin with CORINE at 10 m and 1 m, 16c's inputs and
# options (docs/benchmarks/2026-09-29/bygdin-landcover/run.sh), base c193cb1
# and 15f-3 back to back. Each tree runs from its build-bench/pkg through
# bench.py's child (the editable finder dropped, threads at the CLI default).
# Order B N N B; each block makes two timed --binary runs per tolerance; then
# one --ascii run per side and tolerance for quality() and strip_check.py.
# Usage: bygdin.sh W B S   (S: scratch for meshes; logs go to ./bygdin/)
set -u
W=$1 B=$2 S=$3
PY=$W/.venv/bin/python
R=/Users/skavhaug/projects/rasputin_data
H=$W/docs/benchmarks/2026-09-29/bygdin-landcover
L=$W/docs/benchmarks/2026-10-04/15f-3-acceptance/bygdin
mkdir -p $L $S/bygdin
cd $W
mesh() { # side tree tol tag format
  local side=$1 tree=$2 tol=$3 tag=$4 fmt=$5
  if ! pmset -g batt | grep -q "AC Power"; then echo "NOT ON AC before $tag"; pmset -g batt; exit 3; fi
  pmset -g batt | head -1 > $L/$tag.pmset; sysctl -n vm.swapusage >> $L/$tag.pmset
  /usr/bin/time -l $PY tools/bench.py _child --pkg $tree/build-bench/pkg --threads 0 -- \
    mesh --dem $R/DTM10_UTM33_20260925 --domain $H/bygdin_reduced_t20.geojson --tolerance $tol \
    --features $R/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine \
    --$fmt --out $S/bygdin/$tag.vtk --stats $L/$tag.stats.md > $L/$tag.out 2> $L/$tag.err
  echo "$tag exit=$? $(date -u +%T)" | tee -a $L/exits.txt
  pmset -g batt | head -1 >> $L/$tag.pmset; sysctl -n vm.swapusage >> $L/$tag.pmset
}
block() { # side tree k
  for tol in 10 1; do for i in 1 2; do mesh $1 $2 $tol $1_t${tol}_b$3_run$i binary; done; done
}
: > $L/exits.txt
block base $B 1; block 15f3 $W 1; block 15f3 $W 2; block base $B 2
for tol in 10 1; do mesh base $B $tol base_t${tol}_ascii ascii; mesh 15f3 $W $tol 15f3_t${tol}_ascii ascii; done
echo ALLDONE
