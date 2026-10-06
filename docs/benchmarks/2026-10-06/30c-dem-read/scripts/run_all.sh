#!/bin/bash
# Three repeats per catchment, base and branch alternated (odd repeat: base first; even: branch first).
# A run not on AC throughout is moved to raw/discarded/ and the script stops (exit 3);
# rerunning it skips the runs already kept, so it resumes where it stopped.
D=$(cd "$(dirname "$0")" && pwd)
R=$D/../raw/stats
X=$D/../raw/discarded
for r in 1 2 3; do
  for c in numedalslagen skiensvassdraget; do
    if [ $((r % 2)) = 1 ]; then a=base; b=branch; else a=branch; b=base; fi
    for s in $a $b; do
      o=${s}_${c}_r${r}
      [ -f $R/${o}_stats.md ] && continue
      $D/stats.sh $s $c $r
      if [ $? = 3 ]; then
        mkdir -p $X; t=$(date +%H%M%S)
        for f in $R/${o}_*; do mv $f $X/$(basename $f .${f##*.})_$t.${f##*.}; done
        grep -v "^$o " $R/vtk_sha256.txt > $R/vtk.tmp; mv $R/vtk.tmp $R/vtk_sha256.txt
        echo "NOT ON AC THROUGHOUT: $o (moved to raw/discarded/, suffix $t)"; exit 3
      fi
      echo "$o: $(tail -1 $R/${o}_stderr.txt); $(head -1 $R/${o}_power.txt | cut -c1-60)"
    done
  done
done
echo done
