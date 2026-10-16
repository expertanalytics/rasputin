#!/bin/bash
# 20c-2 re-time after the speed fix 65ae0792 (@perf, 2026-10-07): the 1 m benchmark and thread sweep
# (tools/bench.py run, defaults: 7908_3 tile and quarter, tolerance 1, threads default and 1 to 20,
# 5 repeats, hardening on), order ABBAAB: base = master 41bda81a (20c-1 merged; scratch worktree),
# head = 20c-2 65ae0792 (this worktree). Copied from ../pairs.sh, the master-2764ef71 build dropped.
# bench.py run rebuilds each tree (Release) into <tree>/build-bench before it measures.
# After each run it prints the tin_engine and _core files the build tree's pkg loads, and the
# _core sha256 (run.json's build.so_sha256 is the same file's). Stops off AC.
S=/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c2b
W=/Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality-2
HERE=$W/docs/benchmarks/2026-10-07/20c-2/fix-65ae0792
cd $W
probe() {  # probe TREE
  .venv/bin/python - "$1/build-bench/pkg" <<'EOF'
import hashlib, sys
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, sys.argv[1])
import tin_engine, tin_engine._core as c
print("tin_engine.__file__ =", tin_engine.__file__)
print("_core.__file__ =", c.__file__, "sha256", hashlib.sha256(open(c.__file__, "rb").read()).hexdigest())
EOF
}
for spec in base:r1 head:r1 head:r2 base:r2 base:r3 head:r3; do
  who=${spec%%:*}; r=${spec##*:}
  if [ $who = base ]; then tree=$S/base; else tree=$W; fi
  [ -e "$HERE/bench/b1m-$who-$r/run.json" ] && continue  # done already (a rerun after a stop)
  if ! pmset -g batt | grep -q "AC Power"; then echo "$(date -u +%FT%TZ) STOP: off AC"; exit 3; fi
  echo "$(date -u +%FT%TZ) start b1m-$who-$r $(pmset -g batt | tail -1)"
  .venv/bin/python tools/bench.py run --label b1m-$who-$r --tree $tree --out-root $S/out --mesh-dir $S/meshes1m 2>&1 | tail -2
  rc=${PIPESTATUS[0]}
  echo "$(date -u +%FT%TZ) end b1m-$who-$r exit $rc $(pmset -g batt | tail -1)"
  probe $tree
  if [ $rc -ge 3 ]; then echo "FAILED b1m-$who-$r"; exit 4; fi  # 0, 1, 2: ACCEPTED, REGRESSION, NO BASELINE
  if ! pmset -g batt | grep -q "AC Power"; then echo "$(date -u +%FT%TZ) STOP after b1m-$who-$r: off AC; discard it"; exit 3; fi
  mkdir -p "$HERE/bench/b1m-$who-$r"
  cp $S/out/*/b1m-$who-$r/* "$HERE/bench/b1m-$who-$r/"
done
echo DONE
