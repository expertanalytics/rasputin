#!/bin/bash
# Three repeats per catchment, base and branch alternated (odd repeat: base first; even: branch first).
D=$(dirname "$0")
for r in 1 2 3; do
  for c in numedalslagen skiensvassdraget; do
    if [ $((r % 2)) = 1 ]; then a=base; b=branch; else a=branch; b=base; fi
    $D/stats.sh $a $c $r; $D/stats.sh $b $c $r
  done
done
echo done
