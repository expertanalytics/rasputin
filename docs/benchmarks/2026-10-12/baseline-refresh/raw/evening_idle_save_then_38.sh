#!/bin/zsh
# Evening, 2026-10-12. Wait until the 1-minute load average is under 2.0 on
# two samples 30 s apart (at most 30 minutes), then:
#   run 4: the quick check on this tree (master f09b9863 merged in) with
#          --save-baseline, judged against the stored 357808b1 baseline;
#   run 5: increment 38's tree (worktrees/threads) judged against the
#          baseline run 4 just saved. Nothing else of ours runs beside them.
set -u
TREE=/Users/skavhaug/projects/rasputin/.claude/worktrees/baseline-1012
TREE38=/Users/skavhaug/projects/rasputin/.claude/worktrees/threads
DIR=$TREE/docs/benchmarks/2026-10-12/baseline-refresh
LOG4=$DIR/quick-run-4-idle-save.log
LOG5=$DIR/quick-run-5-increment-38.log
cd "$TREE" || exit 3
deadline=$(( $(date +%s) + 30 * 60 ))
low=0
echo "# run 4 (wait for load < 2.0 twice, then --save-baseline), started $(date '+%Y-%m-%d %H:%M:%S %z')" > "$LOG4"
while [ "$(date +%s)" -lt "$deadline" ]; do
  load=$(sysctl -n vm.loadavg | awk '{print $2}')
  echo "$(date +%H:%M:%S) load1 $load" >> "$LOG4"
  if [ "$(echo "$load < 2.0" | bc)" = 1 ]; then low=$((low + 1)); else low=0; fi
  if [ "$low" -ge 2 ]; then
    { echo "## pmset before"; pmset -g batt; echo "## uptime"; uptime;
      echo "## top CPU users"; ps -eo pcpu,comm | sort -rn | head -6; } >> "$LOG4"
    RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data .venv/bin/python tools/bench_quick.py run --save-baseline >> "$LOG4" 2>&1
    echo "exit $?" >> "$LOG4"
    { echo "## pmset after"; pmset -g batt; echo "## uptime after"; uptime; } >> "$LOG4"
    echo "# run 5 (increment 38, --tree $TREE38, judged against the baseline run 4 saved), started $(date '+%Y-%m-%d %H:%M:%S %z')" > "$LOG5"
    { echo "## tree 38 head"; git -C "$TREE38" rev-parse HEAD; echo "## pmset before"; pmset -g batt; echo "## uptime"; uptime;
      echo "## top CPU users"; ps -eo pcpu,comm | sort -rn | head -6; } >> "$LOG5"
    RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data .venv/bin/python tools/bench_quick.py run --tree "$TREE38" >> "$LOG5" 2>&1
    echo "exit $?" >> "$LOG5"
    { echo "## pmset after"; pmset -g batt; echo "## uptime after"; uptime; } >> "$LOG5"
    exit 0
  fi
  sleep 30
done
echo "GAVE UP: the load never fell under 2.0 on two samples before $(date '+%H:%M:%S')" >> "$LOG4"
exit 4
