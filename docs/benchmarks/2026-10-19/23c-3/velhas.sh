#!/bin/zsh
# 23c-3's Velhas runs: ANADEM from the cache, EPSG:31983, default threads, --binary,
# each from a tree's Release package (build-bench/pkg, built by pairs.sh) under
# /usr/bin/time -l with MallocLargeCache=0, through velhas_driver.py.
# Usage: velhas.sh TREE LABEL S TOL PIECES [THREADS]   (PIECES 0: uncut; THREADS: --threads, named _th<N>)
#   Python: the 23c-3 worktree's venv for every tree (bench.py builds with it).
#   Writes logs/<LABEL>_t<TOL>_p<PIECES>.{err,out,jsonl,pmset} here and the
#   mesh to S/<LABEL>_t<TOL>_p<PIECES>/; the index.json is copied beside the logs.
set -u
T=$1 L=$2 S=$3 TOL=$4 P=$5 TH=${6:-}
O=${0:A:h}
DATA=/Users/skavhaug/projects/rasputin_data/sao_francisco_piece
N=${L}_t${TOL}_p${P}${TH:+_th$TH}
rm -rf $S/$N; mkdir -p $O/logs $S/$N
pmset -g batt | grep -q "AC Power" || { echo "NOT ON AC before $N"; exit 3; }
CUT=(); [ $P != 0 ] && CUT=(--pieces $P)
[ -n "$TH" ] && CUT+=(--threads $TH)
{ date -u +%FT%TZ; pmset -g batt; sysctl vm.swapusage vm.loadavg; } > $O/logs/$N.pmset
MallocLargeCache=0 caffeinate -i /usr/bin/time -l /Users/skavhaug/projects/rasputin/.claude/worktrees/23c-3/.venv/bin/python $O/velhas_driver.py $T/build-bench/pkg \
  $O/logs/$N.jsonl -- mesh --dem anadem-v1 --cache $HOME/.cache/rasputin \
  --domain $DATA/bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983 \
  --tolerance $TOL $CUT --binary --out $S/$N/velhas.vtk \
  > $O/logs/$N.out 2> $O/logs/$N.err
E=$?
{ date -u +%FT%TZ; pmset -g batt; sysctl vm.swapusage vm.loadavg; } >> $O/logs/$N.pmset
[ -f $S/$N/velhas.vtk.pieces/index.json ] && cp $S/$N/velhas.vtk.pieces/index.json $O/logs/$N.index.json
echo "$N exit $E"
