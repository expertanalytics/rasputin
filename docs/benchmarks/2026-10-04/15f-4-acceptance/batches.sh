#!/bin/zsh
# 15f-4 acceptance: base c193cb1 (B), 15f-3 107cdb4 (P), 15f-4 5eeb87a (N).
# Batch 1 B P N, batch 2 N P B. Per slot: bench.py run (both domains, sweep,
# 5 repeats); Bygdin 1 m with CORINE, 3 timed --binary runs under time -l;
# Velhas 5 m (EPSG:31983), 2 timed runs + 1 ascii via velhas.py. Batch 1 also
# writes Bygdin's 1 m ASCII mesh per tree for the hash. Stops off AC.
# Usage: batches.sh W B P S   (W is the 15f-4 worktree, N)
set -u
W=$1 B=$2 P=$3 S=$4
PY=/Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/.venv/bin/python
R=/Users/skavhaug/projects/rasputin_data
H=$W/docs/benchmarks/2026-09-29/bygdin-landcover
E=$W/docs/benchmarks/2026-10-04/15f-4-acceptance
L=$E/bygdin; mkdir -p $L $S/bygdin
cd $W
ac() { pmset -g batt | grep -q "AC Power" || { echo "NOT ON AC before $1"; pmset -g batt; exit 3; } }
bygdin() { # tag tree fmt
  ac $1
  pmset -g batt | head -1 > $L/$1.pmset; sysctl -n vm.swapusage >> $L/$1.pmset
  /usr/bin/time -l $PY tools/bench.py _child --pkg $2/build-bench/pkg --threads 0 -- \
    mesh --dem $R/DTM10_UTM33_20260925 --domain $H/bygdin_reduced_t20.geojson --tolerance 1 \
    --features $R/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine \
    --$3 --out $S/bygdin/$1.vtk --stats $L/$1.stats.md > $L/$1.out 2> $L/$1.err
  echo "$1 exit=$? $(date -u +%T)"
  pmset -g batt | head -1 >> $L/$1.pmset; sysctl -n vm.swapusage >> $L/$1.pmset
}
slot() { # side tree batch baseline
  local side=$1 tree=$2 b=$3 base=$4
  ac $side$b; echo "=== $side b$b start $(date -u +%FT%TZ) $(sysctl -n vm.swapusage)"; pmset -g batt
  $PY tools/bench.py run --label 15f-4-$side-b$b --tree $tree --hardening on --repeats 5 \
     --mesh-dir $S/meshes/$side-b$b ${=base}
  echo "=== bench $side b$b exit $? $(date -u +%T)"
  for i in 1 2 3; do bygdin ${side}_b${b}_run$i $tree binary; done
  [ $b = 1 ] && bygdin ${side}_ascii $tree ascii
  $PY $E/velhas.py run $side-b$b $tree $E/velhas/$side-b$b $S/velhas/$side-b$b --tolerances 5 --repeats 2
  echo "=== $side b$b end $(date -u +%FT%TZ)"; pmset -g batt
}
D=$W/docs/benchmarks/2026-10-04
slot base $B 1 ""
slot 15f3 $P 1 "--baseline $D/15f-4-base-b1"
slot 15f4 $W 1 "--baseline $D/15f-4-15f3-b1"
slot 15f4 $W 2 "--baseline $D/15f-4-15f3-b1"
slot 15f3 $P 2 "--baseline $D/15f-4-15f4-b2"
slot base $B 2 "--baseline $D/15f-4-base-b1"
echo ALLDONE; pmset -g batt
