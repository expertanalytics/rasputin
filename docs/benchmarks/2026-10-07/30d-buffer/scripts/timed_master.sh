#!/bin/zsh
# Step 1 and 2: master f81b20b7 unpatched, untraced. Per catchment, one warm-up
# (tag w0) and three timed repeats (r1-r3), the two feature inputs alternated
# (odd repeat: first input first; even: second first).
D=${0:A:h}
for c in lagan:geojson:gpkg ljungan_flasjo:geojson:gpkg numedalslagen:gpkg33:gpkg; do
  n=${c%%:*}; rest=${c#*:}; a=${rest%%:*}; b=${rest##*:}
  for r in w0 r1 r2 r3; do
    if [[ $r == r2 ]]; then x=$b; y=$a; else x=$a; y=$b; fi
    $D/run_mesh.sh master_$r $n $x; $D/run_mesh.sh master_$r $n $y
  done
done
echo done
