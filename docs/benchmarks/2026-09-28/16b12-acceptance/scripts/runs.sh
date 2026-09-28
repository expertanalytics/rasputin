#!/bin/zsh
# Parts 1, 2 and 4 of the 16b-1/2 acceptance, through the real CLI (`rasputin mesh`).
# Usage: runs.sh <part> <repeats>. Meshes stay in the scratchpad; the logs and
# --stats reports go to the evidence directory.
set -u
cd /Users/skavhaug/projects/rasputin
SP=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12
D=docs/benchmarks/2026-09-28/16b12-acceptance
DATA=/Users/skavhaug/projects/rasputin_data
EU=$DATA/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg
NO=$DATA/corine2018_dtm10_utm33.gpkg
TILE=tests/fixtures/dem_archive/7908_3_10m_z33.tif
QC=docs/benchmarks/2026-09-26/quarter.geojson
EXT=tests/fixtures/corine/clc2018_7908_3.gpkg
CATCH=$D/catchment/catchment_4326.geojson
SQ=$D/catchment/square48_25833.geojson
T6603=$DATA/DTM10_UTM33_20260925/6603_4_10m_z33.tif
mkdir -p $SP/meshes $D/stats
one() {  # one <name> <rep> <args...>
  local name=$1 rep=$2; shift 2
  local out=$SP/meshes/$name.vtk st=$D/stats/$name-r$rep.md
  echo "### $name r$rep $(date +%T) | $(pmset -g batt | head -1 | cut -d"'" -f2)"
  echo "\$ rasputin mesh $*"
  /usr/bin/time -l .venv/bin/rasputin mesh "$@" --out $out --stats $st 2> $SP/err.txt >/dev/null
  echo "exit $?"
  grep -v '^ \+[0-9]\+  [a-z]' $SP/err.txt | grep -v 'real\|user\|sys' 
  grep 'maximum resident set size' $SP/err.txt
  grep -E '^\| (triangles|vertices)|achieved|^\| [0-9]' $st | head -12
  [ $rep = 1 ] && .venv/bin/python $D/scripts/check_vtk.py $out
}
part=$1; reps=$2
for rep in $(seq ${FIRST:-1} $reps); do
case $part in
qc) for tol in 1 10; do
      one qc-t$tol-none $rep --dem $TILE --domain $QC --tolerance $tol
      one qc-t$tol-corine $rep --dem $TILE --domain $QC --tolerance $tol --features $EXT --features-map corine
    done ;;
ola) for tol in 1 10; do
      one ola-t$tol-none $rep --dem $DATA/DTM10_UTM33_20260925 --domain $CATCH --tolerance $tol
      one ola-t$tol-eu3035 $rep --dem $DATA/DTM10_UTM33_20260925 --domain $CATCH --tolerance $tol --features $EU --features-layer U2018_CLC2018_V2020_20u1 --features-map corine
      one ola-t$tol-no25833 $rep --dem $DATA/DTM10_UTM33_20260925 --domain $CATCH --tolerance $tol --features $NO --features-layer corine2018 --features-map corine
    done ;;
sq) for tol in 1 10; do
      one sq-t$tol-none-q25 $rep --dem $T6603 --domain $SQ --tolerance $tol
      one sq-t$tol-none-q0 $rep --dem $T6603 --domain $SQ --tolerance $tol --start-min-angle 0
      one sq-t$tol-eu-q25 $rep --dem $T6603 --domain $SQ --tolerance $tol --features $EU --features-layer U2018_CLC2018_V2020_20u1 --features-map corine
      one sq-t$tol-eu-q0 $rep --dem $T6603 --domain $SQ --tolerance $tol --start-min-angle 0 --features $EU --features-layer U2018_CLC2018_V2020_20u1 --features-map corine
    done ;;
fix)  # 2026-09-29: the Europe-file cases again, after d58d693's query fix
      for tol in 1 10; do
        one ola-t$tol-eu3035-fix $rep --dem $DATA/DTM10_UTM33_20260925 --domain $CATCH --tolerance $tol --features $EU --features-layer U2018_CLC2018_V2020_20u1 --features-map corine
      done
      one sq-t1-eu-q25-fix $rep --dem $T6603 --domain $SQ --tolerance 1 --features $EU --features-layer U2018_CLC2018_V2020_20u1 --features-map corine ;;
esac
done
echo "### done $(date +%T) | $(pmset -g batt | head -1 | cut -d"'" -f2)"
