#!/bin/bash
# usage: runs.sh <tag> <side: master|branch> <args...>; writes $OUT/<tag>-<side>.{vtk,json,err,code}
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/lakes-overrule
OUT=$S/out; mkdir -p $OUT
tag=$1; side=$2; shift 2
if [ "$side" = master ]; then cmd=("$W/.venv/bin/python" "$S/bin/run_master.py"); else cmd=("$W/.venv/bin/rasputin"); fi
cd $S
/usr/bin/time -p "${cmd[@]}" "$@" --out $OUT/$tag-$side.vtk --record $OUT/$tag-$side.json --stats $OUT/$tag-$side.md > $OUT/$tag-$side.out 2> $OUT/$tag-$side.err
echo $? > $OUT/$tag-$side.code
echo "$tag $side exit $(cat $OUT/$tag-$side.code)"
