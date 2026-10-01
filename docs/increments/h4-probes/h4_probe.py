#!/usr/bin/env python3
"""h4 probe: feed the live guards the PreToolUse events they receive, by day and at night.

Evidence for docs/increments/h4-guard-fixes.md §2, not production code. It copies
the hooks into a temporary repository with tests/python/harness_fixtures.py
(the U1 fixture), so the real harness state is never touched, and prints one
line per (case, mode): the decision each guard returns, or `pass`.

    python3 docs/increments/h4-probes/h4_probe.py [<git-common-dir>/harness/queue-*.jsonl ...]

Each queue file given adds its `act` lines as Bash cases (labelled `queue:<n>`).
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests" / "python"))
import harness_fixtures as hf  # noqa: E402

AWAY = "aw" + "ay.py"  # spelled in two parts only so this file's own text stays plain

#: (label, tool, payload, agent_type or None). Bash payload is the command.
CASES: list[tuple[str, str, str, str | None]] = [
    # 1: a subagent's own progress file, Write tool and Bash forms.
    ("fp1-write", "Write", ".claude/current-task/reviewer-231000.md", "reviewer"),
    ("fp1-bash", "Bash",
     "printf 'ask: x\\n' > .claude/current-task/reviewer-231000.md; "
     "cat .claude/REQUIRED-READING.md", "reviewer"),
    # 2: the words "git push" inside text written to a file, or a PR body file.
    ("fp2-heredoc", "Bash",
     "cat > .claude/current-task/session.md <<'EOF'\n"
     "ASK OLA: git push -u origin x: push + open PR\nEOF", None),
    ("fp2-body-file", "Bash",
     "cat > $CLAUDE_JOB_DIR/tmp/body.md <<'EOF'\nAfter merge, git push the tag.\nEOF", None),
    ("fp2-commit-msg", "Bash", "git commit -q -F $CLAUDE_JOB_DIR/tmp/msg.txt "
     "&& echo 'next: git push once Ola says yes'", None),
    # 3: rule-file names in the text, an ordinary file written.
    ("fp3-sed", "Bash",
     "sed -i '' 's/a/b/' docs/increments/23-basin-scale.md; "
     "python3 tools/check_citations.py | tail -3", "architect"),
    ("fp3-heredoc", "Bash",
     "python3 - <<'EOF'\np='docs/increments/23-basin-scale.md'; s=open(p).read()\n"
     "s=s.replace('see CLAUDE.md', 'see docs/increments/README.md')\n"
     "open(p,'w').write(s)\nEOF", "architect"),
    ("fp3-commit-msg", "Bash",
     "printf 'docs: cite CLAUDE.md section 2\\n' > $CLAUDE_JOB_DIR/tmp/m.txt "
     "&& git commit -q -F $CLAUDE_JOB_DIR/tmp/m.txt", None),
    # 4: reads treated as writes.
    ("fp4-git-show", "Bash",
     "git show agents-web-tools:.claude/settings.json | python3 -c "
     "\"import json,sys; print(json.load(sys.stdin).keys())\"", "reviewer"),
    ("fp4-diff-to-tmp", "Bash",
     "python3 tools/check_citations.py > /dev/null; git diff -U0 HEAD > "
     "$CLAUDE_JOB_DIR/tmp/diff.txt; python3 - <<'EOF'\nprint(1)\nEOF", "developer"),
    # 5: the away script's name in text that does not run it.
    # 5: the away script's name in text that does not run it (real refused commands).
    ("fp5-grep", "Bash",
     f"grep -n -B3 -A12 'def _open_tty\\|/dev/tty' tools/{AWAY} | head -60; which -a python3",
     None),
    ("fp5-echo-memory", "Bash",
     f"echo '- each night is a trial of {AWAY} + guards; report idle time' >> MEMORY.md", None),
    ("fp5-printf-task", "Bash",
     f"printf 'Ask: spec amendment for {AWAY} /dev/tty defect; commit, no push.\\n' "
     "> .claude/current-task/architect-224000.md", "architect"),
    # True positives that must stay as they are.
    ("tp-push", "Bash", "git push -u origin x 2>&1 | tail -2", None),
    ("tp-rebase", "Bash", "cd w && git rebase master 2>&1 | tail -3", None),
    ("tp-claude-md-sed", "Bash", "sed -i '' 's/a/b/' CLAUDE.md", None),
    ("tp-readme-heredoc", "Bash",
     "cat > docs/increments/README.md <<'EOF'\nx\nEOF", None),
    ("tp-settings-py", "Bash",
     "python3 -c \"open('.claude/settings.json','w').write('{}')\"", None),
    ("tp-check-cp", "Bash", "cp /tmp/x.py tools/check_citations.py", None),
    ("tp-edit-claude", "Edit", "CLAUDE.md", None),
    ("tp-run-away", "Bash", f"python3 tools/{AWAY} --back", None),
    ("tp-harness-write", "Bash", "echo {} > .git/" + "harness/unattended.json", None),
]


def decide(repo: Path, tool: str, payload: str, agent: str | None) -> str:
    extra = {"agent_id": "p1", "agent_type": agent} if agent else {}
    if tool == "Bash":
        event = hf.bash_event(repo, payload, **extra)
        hooks = (hf.GUARD_GOVERNANCE, hf.GUARD_PUSH)
    else:
        event = hf.file_event(repo, tool, payload, **extra)
        hooks = (hf.GUARD_GOVERNANCE,)
    out = []
    for hook in hooks:
        found = hf.pretool_decision(hf.run_script(repo, hook, event))
        if found is not None:
            name = hook.rsplit("/", 1)[-1].removesuffix(".py").removeprefix("guard_")
            out.append(f"{name}:{found[0]}")
    return " ".join(out) or "pass"


def main(argv: list[str]) -> int:
    cases = list(CASES)
    for queue in argv:
        for n, line in enumerate(Path(queue).read_text().splitlines()):
            entry = json.loads(line)
            cases.append((f"queue:{n}", "Bash", entry["act"], entry.get("agent_type")))
    with tempfile.TemporaryDirectory() as tmp:
        repo = hf.make_repo(Path(tmp) / "repo")
        for mode in ("off", "on"):
            if mode == "on":
                hf.set_mode(repo, "on")
            for label, tool, payload, agent in cases:
                print(f"{mode:3} {label:18} {decide(repo, tool, payload, agent)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
