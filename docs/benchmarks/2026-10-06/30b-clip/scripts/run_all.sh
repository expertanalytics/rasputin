#!/bin/bash
# Three repeats per catchment, base and branch alternated (odd repeat: base first; even: branch first).
# Stops at the first run not on AC power.
D=$(cd "$(dirname "$0")" && pwd)
R=$D/../raw/stats
for r in 1 2 3; do
  for c in numedalslagen skiensvassdraget; do
    if [ $((r % 2)) = 1 ]; then a=base; b=branch; else a=branch; b=base; fi
    for s in $a $b; do
      $D/stats.sh $s $c $r
      if grep -qv "AC Power" $R/${s}_${c}_r${r}_power.txt; then echo "NOT ON AC: $s $c $r"; exit 3; fi
      echo "$s $c r$r: $(tail -1 $R/${s}_${c}_r${r}_stderr.txt)"
    done
  done
done
echo done
