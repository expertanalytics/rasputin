#!/bin/zsh
# Increment 46 PR 7: `rasputin snow` on Bygdin (DTM10, tolerance-10 mesh) at
# the new defaults (153 days, 1 April to 1 September, warming held from
# 1 August, one bin), shadows on and off, 8 threads, interleaved repeats.
# Usage: run_pr7.sh <venv bin dir> <bygdin.vtk> <DTM10 dir> <out dir for "on"> <out dir for "off"> <raw dir> <repeats>
set -u
B=$1; MESH=$2; DTM=$3; ON=$4; OFF=$5; RAW=$6; REPEATS=$7
TSV=$RAW/runs.tsv
print -- "case\tthreads\trepeat\twall_s\tmax_rss_bytes\tpeak_footprint_bytes\tout_bytes\tframes\texit\tpower\tload" > $TSV
power() { pmset -g batt | head -1 | sed "s/.*'\(.*\)'.*/\1/" }
load() { sysctl -n vm.loadavg | awk '{print $2}' }
one() {  # case out repeat args...
  local c=$1 o=$2 i=$3; shift 3
  rm -rf "$o"
  local p=$(power) l=$(load)
  /usr/bin/time -l "$B/rasputin" snow --mesh "$MESH" --dem "$DTM" "$@" --out "$o" --threads 8 \
    > "$RAW/log-$c-r$i.txt" 2> "$RAW/time-$c-r$i.txt"
  local rc=$?
  local wall=$(awk '/ real /{print $1}' "$RAW/time-$c-r$i.txt")
  local rss=$(awk '/maximum resident set size/{print $1}' "$RAW/time-$c-r$i.txt")
  local peak=$(awk '/peak memory footprint/{print $1}' "$RAW/time-$c-r$i.txt")
  local size=$(find "$o" -type f -exec stat -f %z {} + | awk '{s+=$1} END{print s}')
  local frames=$(ls "$o" | grep -c '^swe[0-9]*\.var$')
  print -- "$c\t8\t$i\t$wall\t$rss\t$peak\t$size\t$frames\t$rc\t$p\t$l" >> $TSV
  print -r -- "$c r$i: $wall s rc $rc frames $frames ($p)"
}
for i in $(seq 1 "$REPEATS"); do
  one bygdin-on "$ON" $i --shadows on
  one bygdin-off "$OFF" $i --shadows off
done
