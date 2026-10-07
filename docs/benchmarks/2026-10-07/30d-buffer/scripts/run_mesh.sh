#!/bin/zsh
# One `rasputin mesh --stats` run through launch.py.
# usage: run_mesh.sh <tag> <catchment> <features> [launch.py flags...]
#   catchment: lagan | ljungan_flasjo | numedalslagen
#   features:  geojson (the CORINE window cut for the catchment, EPSG:3035)
#              gpkg    (the European CORINE GeoPackage, EPSG:3035, as Ola's mesh_sweden.sh)
#              gpkg33  (Numedalslagen only: corine2018_dtm10_utm33.gpkg, as the 2026-10-06 profiles)
# Stats, stderr, power before and after and the .vtk's SHA-256 go to raw/runs/;
# the .vtk goes to $SCRATCH/out (not committed). PKG and SCRATCH come from the environment.
set -u
D=${0:A:h}; R=${D:h}/raw/runs
DATA=/Users/skavhaug/projects/rasputin_data
PY=/Users/skavhaug/projects/rasputin/.venv/bin/python
export RASPUTIN_DATA=$DATA
tag=$1; c=$2; f=$3; shift 3
mkdir -p $R $SCRATCH/out
o=${tag}_${c}_${f}
EU=$DATA/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg
case $f in
  geojson) F=(--features $DATA/sweden_corine/${c}_clc2018_3035.geojson --features-crs EPSG:3035 --features-map corine);;
  gpkg) F=(--features $EU --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine);;
  gpkg33) F=(--features $DATA/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine);;
esac
case $c in
  lagan) M=(--dem glo30 --cache $DATA/sweden_glo30_cache --out-crs EPSG:3006 --domain $DATA/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson);;
  ljungan_flasjo) M=(--dem glo30 --cache $DATA/sweden_glo30_cache --out-crs EPSG:3006 --domain $DATA/sweden_smhi_svar/ljungan_flasjo_svar2022_3006.geojson);;
  numedalslagen) M=(--dem $DATA/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin_scratch/norway/numedalslagen/numedalslagen_outline_nve.geojson);;
esac
pmset -g batt | head -1 > $R/${o}_power.txt
PKG=$PKG $PY $D/launch.py "$@" -- mesh $M $F --tolerance 10 --binary \
  --out $SCRATCH/out/${o}.vtk --stats $R/${o}_stats.md 2> $R/${o}_stderr.txt
echo "exit $?" >> $R/${o}_stderr.txt
pmset -g batt | head -1 >> $R/${o}_power.txt
shasum -a 256 $SCRATCH/out/${o}.vtk | awk -v o=$o '{print o, $1}' >> $R/vtk_sha256.txt
