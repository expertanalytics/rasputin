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
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HOOK_ROOT = Path(__file__).resolve().parent.parent.parent
# Appended, never inserted first: a tools/ file named after a stdlib module
# would otherwise replace that module for every later import (h19 §3).
sys.path.append(str(HOOK_ROOT / "tools"))
try:
    import shell_scan
except ImportError:
    shell_scan = None

#: The git subcommands that make commits on the checked-out branch (h19 §3.1).
COMMITTING = frozenset({"commit", "merge", "pull", "cherry-pick", "revert", "am", "rebase"})
#: git's own options whose next word is their argument; copied from guard_push.py, not
#: imported, so that one hook failing to load cannot take the other with it.
GIT_TAKES_ARG = frozenset(
    {"-C", "-c", "--git-dir", "--work-tree", "--namespace", "--attr-source", "--config-env"}
)
#: The text test for a line shell_scan cannot read; `(?![\w-])` so `git merge-base` is no merge.
TEXT_TEST = re.compile(
    r"\bgit(\s+-\S+(\s+[^-\s]\S*)?)*\s+(commit|merge|pull|cherry-pick|revert|am|rebase)(?![\w-])"
)
NOT_READ = "this command line could not be read"
NO_SCAN = "tools/shell_scan.py could not be imported, so this command line was not read"


def git(cwd: Path, *args: str) -> str:
    try:
        return subprocess.run(
            ["git", *args], cwd=cwd, capture_output=True, text=True, timeout=10
        ).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        return ""


def step(here: Path | str, word: str, known: bool = True) -> Path | str:
    """The directory after `cd word` from `here`; a str is unknown, the first unknown word."""
    if not known or word == "-" or "$" in word or "`" in word:
        return here if isinstance(here, str) else word
    path = Path(os.path.expanduser(word))
    if path.is_absolute():
        return path
    return here / path if isinstance(here, Path) else here


def git_dirs(command: str, cwd: Path) -> list[Path | str] | None:
    """The directory each committing git command on the line runs in, in order:
    a Path where h19 §3.2 resolves it, the word as written where it does not.
    [] when no command commits; None when shell_scan cannot read the line
    or failed to import (the hook then uses the text test).

    Pure: joins paths and expands `~`, nothing else. The shell of a subagent
    returns to the session directory between calls, so a worktree is reached
    only through `cd` or `git -C` on the line itself; following them is how the
    committed tree is found (D3 in the 2026-10-06 retrospective).
    """
    simples = shell_scan.parse(command) if shell_scan else None
    if simples is None:
        return None
    here: Path | str = cwd
    found: list[Path | str] = []
    for argv in (simple.argv for simple in simples if simple.argv):
        name, rest = shell_scan.base(argv[0]), argv[1:]
        if name in ("cd", "pushd", "popd"):
            words = [w for w in rest if w == "-" or not w.startswith("-")]
            here = step(here, words[0], name != "popd") if words else step(here, name, False)
        elif name == "git":
            where, pinned = here, False  # pinned: --git-dir or --work-tree seen, -C cannot undo it
            while rest and rest[0].startswith("-"):
                option, takes = rest[0], rest[0] in GIT_TAKES_ARG and len(rest) > 1
                if option == "-C" and takes and not pinned:
                    where = step(where, rest[1])
                elif option.split("=")[0] in ("--git-dir", "--work-tree"):
                    where, pinned = (f"{option} {rest[1]}" if takes else option), True
                rest = rest[2:] if takes else rest[1:]
            if rest and rest[0] in COMMITTING:
                found.append(where)
    return found


def committed_trees(event: dict) -> list[Path | str]:
    """The distinct work-tree tops to gate, in first-seen order; a str is a
    NOT CHECKED entry (h19 §3.2). [] means the hook does nothing.

    Not the hook's own tree: the settings run `$CLAUDE_PROJECT_DIR/.claude/hooks/...`,
    which is the MAIN checkout even for a commit made in a worktree.
    """
    cwd = Path(event.get("cwd") or HOOK_ROOT)
    command = (event.get("tool_input", {}) or {}).get("command", "")
    dirs = git_dirs(command, cwd)
    notes: list[Path | str] = []
    if dirs is None:  # Unreadable: gate cwd's tree as before h19, and say it may be the wrong one.
        dirs = [cwd] if TEXT_TEST.search(command) else []
        notes = [NOT_READ if shell_scan else NO_SCAN] if dirs else []
    trees: list[Path | str] = []
    for where in dirs:
        top = ""
        if isinstance(where, Path) and where.is_dir():
            top = git(where, "rev-parse", "--show-toplevel")
        entry = Path(top) if top else str(where)
        if entry not in trees:
            trees.append(entry)
    # The note says gates ran in cwd; true only when cwd resolved to a work tree.
    return trees + notes if all(isinstance(t, Path) for t in trees) else trees


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
    blocks = []
    for root in committed_trees(event):
        if root in (NOT_READ, NO_SCAN):
            blocks.append(f"NOT CHECKED: {root}, so the gates ran in the session's directory, "
                          "which may not be the tree committed; run them there yourself.")
        elif isinstance(root, str):
            blocks.append(f"NOT CHECKED: git ran in `{root}`, which this hook cannot resolve "
                          "to a work tree; run the gates there yourself.")
        elif failures := run_gates(root):
            blocks.append(f"Gates in {root}:\n" + "\n\n".join(failures))

    if not blocks:
        return 0

    print(
        "The commit is in, and below are the gates that are red or the trees this hook did "
        "not check. Fix and amend -- do not push or open a PR on this.\n\n" + "\n\n".join(blocks)
        + "\n\nThis is not the full set: mypy and ctest need a build, and CI "
          "decides. Run `gh pr checks <pr>` before calling anything merge-ready.",
        file=sys.stderr,
    )
    return 2  # 2 feeds stderr back to the agent.


def run_gates(root: Path) -> list[str]:
    """Each gate that failed or could not run in `root`, with its own output."""
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
    return failures


if __name__ == "__main__":
    sys.exit(main())
