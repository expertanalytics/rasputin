#!/usr/bin/env python3
"""PostToolUse: run the gates after a commit, and report their OUTPUT.

Two failures this answers, both from the same root. The agent has more than once
read a command's SHAPE instead of its CONTENT: it committed with `ruff` failing
because the summary looked like a pass, and -- on 2026-09-23, minutes after
writing a principle about exactly this -- it probed the dependency gate with
`2>/dev/null` and kept only `exit=1`, discarding the findings that said WHICH
dependency was caught. Both times the information was on the screen or one flag
away.

A hook cannot make an agent read. It can make the reading unnecessary by putting
the gates' own words into the transcript after every commit, unsuppressed and
unsummarised.

This fires AFTER the commit, which is deliberate. A pre-commit block would stall work on a WIP
commit, and there is no pre-commit slot for a PostToolUse check anyway.
A commit is cheap to amend; the cost of finding out in CI is not.

NOT AUTHORITATIVE. CLAUDE.md is explicit that CI decides, and a whole branch once
merged with CI red on every commit. `mypy` and the C++ suite are omitted here
because they need a build; green from this hook means "the GATES below passed on
these files", which is strictly less than `gh pr checks`.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

HOOK_ROOT = Path(__file__).resolve().parent.parent.parent


def git(cwd: Path, *args: str) -> str:
    try:
        return subprocess.run(
            ["git", *args], cwd=cwd, capture_output=True, text=True, timeout=10
        ).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        return ""


def committed_tree(event: dict) -> Path:
    """The tree the commit was made in, which is not always this file's tree.

    The settings run `$CLAUDE_PROJECT_DIR/.claude/hooks/...`, and in a session
    isolated in a worktree (.claude/worktrees/<name>) that variable is the MAIN
    checkout. On activation (2026-09-29) this hook therefore gated the main
    checkout's branch after a commit made in a worktree. The event's `cwd` is
    where the command ran; its git toplevel is the tree that was committed.
    """
    cwd = Path(event.get("cwd") or HOOK_ROOT)
    top = git(cwd, "rev-parse", "--show-toplevel") if cwd.is_dir() else ""
    return Path(top) if top else HOOK_ROOT


def find_ruff(root: Path) -> str:
    """The tree's .venv, else the main checkout's, else PATH.

    A worktree has no .venv of its own; looking only in root/.venv made every
    worktree commit report both ruff gates as DID NOT RUN (found on activation,
    2026-09-29). If none is found the literal path is returned, so the gate
    still reports DID NOT RUN rather than passing silently.
    """
    own = root / ".venv" / "bin" / "ruff"
    common = git(root, "rev-parse", "--path-format=absolute", "--git-common-dir")
    main = Path(common).parent / ".venv" / "bin" / "ruff" if common else own
    for candidate in (own, main):
        if candidate.is_file():
            return str(candidate)
    return shutil.which("ruff") or str(own)


def gates(root: Path) -> tuple[list[str], ...]:
    """Fast, no build required, no network.

    Anything needing the extension belongs in CI -- a hook slow enough to
    resent is a hook that gets removed.

    ruff is invoked from a .venv directly and NOT as `uv run ruff`. `uv run`
    regenerates uv.lock, and uv.lock is gitignored on this branch precisely
    because 1016 lines of it rode into a commit that way. A hook scheduling
    that command on every commit would recreate the artifact forever. CI runs
    plain `ruff check .` and `ruff format --check .` (.github/workflows/main.yaml).
    """
    ruff = find_ruff(root)
    return (
        ["python3", "tools/check_prohibited_deps.py"],
        ["python3", "tools/check_detria_boundary.py"],
        ["python3", "tools/check_citations.py"],
        [ruff, "check", "."],
        [ruff, "format", "--check", "."],
    )


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    if event.get("tool_name") != "Bash":
        return 0
    command = (event.get("tool_input", {}) or {}).get("command", "")
    if "git commit" not in command and "git merge" not in command:
        return 0

    root = committed_tree(event)
    failures = []
    for gate in gates(root):
        try:
            done = subprocess.run(gate, cwd=root, capture_output=True, text=True, timeout=180)
        except (OSError, subprocess.TimeoutExpired) as exc:
            # A gate that cannot run is a gate that did not pass. Reporting it as
            # silence is how `ruff: command not found` reads as green.
            failures.append(f"$ {' '.join(gate)}\n  DID NOT RUN: {exc}")
            continue
        if done.returncode != 0:
            body = (done.stdout + done.stderr).strip()
            failures.append(f"$ {' '.join(gate)}  (exit {done.returncode})\n{body}")

    if not failures:
        return 0

    print(
        "The commit is in, and these gates are red. Fix and amend -- do not "
        "push or open a PR on this.\n\n" + "\n\n".join(failures)
        + "\n\nThis is not the full set: mypy and ctest need a build, and CI "
          "decides. Run `gh pr checks <pr>` before calling anything merge-ready.",
        file=sys.stderr,
    )
    return 2  # 2 feeds stderr back to the agent.


if __name__ == "__main__":
    sys.exit(main())
