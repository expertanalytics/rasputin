#!/bin/zsh
# Increment 46 PR 6 acceptance: `rasputin snow` with cast shadows on Bygdin
# (DTM10, tolerance-10 mesh, 122 days) and on the synthetic 2 m mesh
# (91 days), at 1 and 8 threads; then the shadow parts timed one by one.
# Usage: run_pr6.sh <venv bin dir> <work dir with bygdin.vtk, mesh2.vtk, synth.tif> <DTM10 dir> <raw dir> <repeats>
set -u
B=$1; W=$2; DTM=$3; RAW=$4; REPEATS=$5
HERE=${0:A:h}
cd "$W" || exit 1
TSV=$RAW/runs.tsv
print -- "case\tthreads\trepeat\twall_s\tmax_rss_bytes\tpeak_footprint_bytes\tout_bytes\texit\tpower\tload" > $TSV
power() { pmset -g batt | head -1 | sed "s/.*'\(.*\)'.*/\1/" }
load() { sysctl -n vm.loadavg | awk '{print $2}' }
one() {  # case threads repeat args...
  local c=$1 t=$2 i=$3; shift 3
  local o="out-$c-t$t"
  rm -rf "$o"
  local p=$(power) l=$(load)
  /usr/bin/time -l "$B/rasputin" snow "$@" --out "$o" --threads "$t" \
    > "$RAW/log-$c-t$t-r$i.txt" 2> "$RAW/time-$c-t$t-r$i.txt"
  local rc=$?
  local wall=$(awk '/ real /{print $1}' "$RAW/time-$c-t$t-r$i.txt")
  local rss=$(awk '/maximum resident set size/{print $1}' "$RAW/time-$c-t$t-r$i.txt")
  local peak=$(awk '/peak memory footprint/{print $1}' "$RAW/time-$c-t$t-r$i.txt")
  local size=$(find "$o" -type f -exec stat -f %z {} + | awk '{s+=$1} END{print s}')
  print -- "$c\t$t\t$i\t$wall\t$rss\t$peak\t$size\t$rc\t$p\t$l" >> $TSV
  print -r -- "$c t$t r$i: $wall s rc $rc ($p)"
}
for i in $(seq 1 "$REPEATS"); do
  for t in 1 8; do
    one bygdin-on $t $i --mesh bygdin.vtk --dem "$DTM" --shadows on --days 122
    one synth2-on $t $i --mesh mesh2.vtk --dem synth.tif --shadows on
  done
done
for t in 1 8; do
  one bygdin-off $t 1 --mesh bygdin.vtk --dem "$DTM" --shadows off --days 122
  one synth2-off $t 1 --mesh mesh2.vtk --shadows off
done
for i in $(seq 1 "$REPEATS"); do
  for t in 1 8; do
    "$B/python" "$HERE/components.py" bygdin.vtk "$DTM" 122 7 $t "$RAW/components-bygdin-t$t-r$i.json" > /dev/null
    print -r -- "components bygdin t$t r$i ($(power))"
    "$B/python" "$HERE/components.py" mesh2.vtk synth.tif 91 7 $t "$RAW/components-synth2-t$t-r$i.json" > /dev/null
    print -r -- "components synth2 t$t r$i ($(power))"
  done
done
"$B/python" "$HERE/search_reach.py" bygdin.vtk "$DTM" 2000 "$RAW/search_reach.json" > "$RAW/search_reach.out"
print -r -- "search_reach done ($(power))"
