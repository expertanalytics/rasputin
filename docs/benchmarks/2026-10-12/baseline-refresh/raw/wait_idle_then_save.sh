#!/bin/zsh
# Wait until the 1-minute load average is under 1.5 on two samples 30 s apart
# (at most 45 minutes), then run the quick check with --save-baseline once.
set -u
TREE=/Users/skavhaug/projects/rasputin/.claude/worktrees/baseline-1012
LOG=$TREE/docs/benchmarks/2026-10-12/baseline-refresh/quick-run-3-idle-save.log
cd "$TREE" || exit 3
deadline=$(( $(date +%s) + 45 * 60 ))
low=0
echo "# run 3 (wait for load < 1.5, then --save-baseline), started $(date '+%Y-%m-%d %H:%M:%S %z')" > "$LOG"
while [ "$(date +%s)" -lt "$deadline" ]; do
  load=$(sysctl -n vm.loadavg | awk '{print $2}')
  echo "$(date +%H:%M:%S) load1 $load" >> "$LOG"
  if [ "$(echo "$load < 1.5" | bc)" = 1 ]; then low=$((low + 1)); else low=0; fi
  if [ "$low" -ge 2 ]; then
    { echo "## pmset before"; pmset -g batt; echo "## uptime"; uptime; } >> "$LOG"
    RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data .venv/bin/python tools/bench_quick.py run --save-baseline >> "$LOG" 2>&1
    echo "exit $?" >> "$LOG"
    { echo "## pmset after"; pmset -g batt; echo "## uptime after"; uptime; } >> "$LOG"
    exit 0
  fi
  sleep 30
done
echo "GAVE UP: the load never fell under 1.5 on two samples before $(date '+%H:%M:%S')" >> "$LOG"
exit 4
