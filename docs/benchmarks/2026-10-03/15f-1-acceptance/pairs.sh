#!/bin/zsh
set -u
PY=/Users/skavhaug/projects/rasputin/.venv/bin/python
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/agent-a0cb49bea16623b07; B=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/base-390b516; S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad
D=$W/docs/benchmarks/2026-10-03
cd $W
TH=1,2,4,6,8,10,16,20
run() { # label tree baseline
  local lab=$1 tree=$2 base=$3
  echo "=== $lab start $(date -u +%FT%TZ)"; pmset -g batt
  PYTHONPATH=$tree/build-bench/pkg $PY tools/bench.py run --label $lab --tree $tree --threads $TH --repeats 5 \
     --mesh-dir $S/meshes/$lab ${=base}
  echo "=== $lab exit $? $(date -u +%FT%TZ)"
}
run 15f-1-base-390b516 $B ""
run 15f-1 $W "--baseline $D/15f-1-base-390b516"
run 15f-1-r2 $W "--baseline $D/15f-1-base-390b516"
run 15f-1-base-390b516-r2 $B "--baseline $D/15f-1-r2"
echo ALLDONE
