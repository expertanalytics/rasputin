#!/bin/zsh
# Increment 46 PR 3 acceptance: `rasputin snow` on the 1 m synthetic mesh,
# default spring, at 1 and 8 threads, interleaved, REPEATS times each.
# Usage: run_snow.sh <rasputin executable> <work dir holding mesh1.vtk> <repeats>
# Writes one line per run to <work dir>/runs.tsv; keeps the last run's output
# of each thread count as out-t<N>/.
set -u
R=$1; W=$2; REPEATS=$3
cd "$W" || exit 1
print -r -- "threads\trepeat\twall_s\tmax_rss_bytes\tpeak_footprint_bytes\tout_bytes\texit" > runs.tsv
for i in $(seq 1 "$REPEATS"); do
  for t in 1 8; do
    rm -rf "out-t$t"
    /usr/bin/time -l "$R" snow --mesh mesh1.vtk --out "out-t$t" --threads "$t" \
      > "log-t$t-r$i.txt" 2> "time-t$t-r$i.txt"
    rc=$?
    wall=$(awk '/ real /{print $1}' "time-t$t-r$i.txt")
    rss=$(awk '/maximum resident set size/{print $1}' "time-t$t-r$i.txt")
    peak=$(awk '/peak memory footprint/{print $1}' "time-t$t-r$i.txt")
    size=$(find "out-t$t" -type f -exec stat -f %z {} + | awk '{s+=$1} END{print s}')
    print -r -- "$t\t$i\t$wall\t$rss\t$peak\t$size\t$rc" >> runs.tsv
  done
done
cat runs.tsv
