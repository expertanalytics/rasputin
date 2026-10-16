#!/bin/zsh
# Step 3: master f81b20b7 with the candidate fixes patched in by launch.py,
# untraced. Per catchment and feature input, one warm-up (w0) and three timed
# repeats (r1-r3); the variants rotate their order each repeat.
#   p1      buffer once (shared)
#   p12     + the features region grown from the convex hull
#   p123k*  + the outline grown in pieces of 1000, 500 or 250 edges
# usage: timed_candidates.sh [catchment:input ...]  (default: Lagan and Ljungan, both inputs)
# A run whose stats file exists is skipped, so an interrupted sweep resumes in order.
D=${0:A:h}; R=${D:h}/raw/runs
variants=(p1:1:1000 p12:1,2:1000 p123k1000:1,2,3:1000 p123k500:1,2,3:500 p123k250:1,2,3:250)
pairs=(${@:-lagan:geojson lagan:gpkg ljungan_flasjo:geojson ljungan_flasjo:gpkg})
for cf in $pairs; do
  c=${cf%%:*}; f=${cf##*:}
  i=0
  for r in w0 r1 r2 r3; do
    for j in 1 2 3 4 5; do
      v=${variants[$(( (j + i - 1) % 5 + 1 ))]}
      tag=${v%%:*}; rest=${v#*:}; p=${rest%%:*}; k=${rest##*:}
      [[ -e $R/${tag}_${r}_${c}_${f}_stats.md ]] && continue
      $D/run_mesh.sh ${tag}_$r $c $f --patch $p --pieces $k
    done
    i=$((i + 1))
  done
done
echo done
