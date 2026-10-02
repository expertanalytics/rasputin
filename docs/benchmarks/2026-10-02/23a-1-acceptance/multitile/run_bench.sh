set -x
cd /Users/skavhaug/projects/rasputin/.claude/worktrees/23a-1
pmset -g batt
/Users/skavhaug/projects/rasputin/.venv/bin/python /Users/skavhaug/projects/rasputin/.claude/worktrees/23a-1/tools/bench.py run --label 23a1-base-5a57793 --tree /Users/skavhaug/.claude/jobs/85c14e7c/tmp/master --out-root /Users/skavhaug/.claude/jobs/85c14e7c/tmp/bench
echo base_exit=$?
pmset -g batt
/Users/skavhaug/projects/rasputin/.venv/bin/python /Users/skavhaug/projects/rasputin/.claude/worktrees/23a-1/tools/bench.py run --label 23a1-3b729fb --tree /Users/skavhaug/projects/rasputin/.claude/worktrees/23a-1 --out-root /Users/skavhaug/.claude/jobs/85c14e7c/tmp/bench --baseline /Users/skavhaug/.claude/jobs/85c14e7c/tmp/bench/2026-10-02/23a1-base-5a57793
echo new_exit=$?
pmset -g batt
