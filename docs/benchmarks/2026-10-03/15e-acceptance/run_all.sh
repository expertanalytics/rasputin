#!/bin/bash
# The 15e memory runs: sub-basin 761 at 10 m, base (7810cf8) and branch, as
# shipped and with MallocLargeCache=0, two rounds, interleaved. A measurement
# script. Usage: run_all.sh VENV_BASE VENV_NEW DATA OUT_MESH_DIR
set -u
VB=$1; VN=$2; D=$3; OUT=$4
cd "$(dirname "$0")"; mkdir -p runs "$OUT"
WKT="$(cat basin_out_crs.wkt)"
for r in 1 2; do
  for cfg in base new; do
    for nlc in "" _nolargecache; do
      tag=761_t10_${cfg}${nlc}_r$r
      PY=$VB/bin/python; [ $cfg = new ] && PY=$VN/bin/python
      env_=(); [ -n "$nlc" ] && env_=(MallocLargeCache=0)
      env ${env_[@]+"${env_[@]}"} $PY run_phases.py runs $tag -- mesh --dem anadem-v1 --cache $D/cache \
        --domain $D/sao_francisco_piece/bho2017_50k_level3/761_epsg4674.geojson \
        --out-crs "$WKT" --tolerance 10 --binary --out $OUT/$tag.vtk --stats runs/$tag.stats.md
      $PY counts.py runs $tag $OUT/$tag.vtk
    done
  done
done
