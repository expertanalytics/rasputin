# usage: run_b.sh ROUNDS ; interleaves master and branch, networking denied by sandbox-exec
D=/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925
BOX="799750 7849750 900250 7950250"
pmset -g batt
for r in $(seq 1 $1); do
  for side in base new; do
    if [ $side = base ]; then T=/Users/skavhaug/.claude/jobs/85c14e7c/tmp/master; else T=/Users/skavhaug/projects/rasputin/.claude/worktrees/23a-1; fi
    echo "== $side round $r $(date +%T)"
    /usr/bin/time -l sandbox-exec -p '(version 1)(allow default)(deny network*)' \
      /Users/skavhaug/projects/rasputin/.venv/bin/python /Users/skavhaug/projects/rasputin/.claude/worktrees/23a-1/tools/bench.py _child --pkg $T/build-bench/pkg --threads 0 -- \
      mesh --dem $D --bbox $BOX --tolerance 1 --out /Users/skavhaug/.claude/jobs/85c14e7c/tmp/b/${side}_$r.vtk --binary \
      > /Users/skavhaug/.claude/jobs/85c14e7c/tmp/b/${side}_$r.out 2> /Users/skavhaug/.claude/jobs/85c14e7c/tmp/b/${side}_$r.err
    echo "exit=$?"
    tail -25 /Users/skavhaug/.claude/jobs/85c14e7c/tmp/b/${side}_$r.err | grep -E 'real|maximum resident|refine_s'
  done
done
pmset -g batt
