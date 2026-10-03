"""Shared fixture for the unattended-mode harness suites (h3 U1).

`docs/increments/h3-unattended-u1.md` §4 fixes the fixture: a temporary git
repository with one commit, holding copies of the hooks and tools at their real
relative paths. A hook copy computes its checkout root from its own `__file__`,
so its state dir is `<tmp>/.git/harness` and the real repository's is never
read or written. Nothing here runs the real `tools/away.py`.

`tools/` is not a package, so pure functions are loaded from their path, as
`test_session_state.py` does.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REAL = Path(__file__).resolve().parents[2]

#: The scripts §4 copies into the temporary repository, at the same paths.
COPIED = (
    ".claude/hooks/guard_push.py",
    ".claude/hooks/guard_governance.py",
    ".claude/hooks/guard_unattended.py",
    "tools/harness_mode.py",
    "tools/away.py",
    "tools/session_state.py",
    # h4: the guards judge a command's targets through this parser
    # (docs/increments/h4-guard-fixes.md §3). Copied when present, like the rest.
    "tools/shell_scan.py",
    # h8: the recap prints the size table through this module
    # (docs/increments/h8-window-and-recap.md §3.6).
    "tools/rule_sizes.py",
)

GUARD_PUSH = ".claude/hooks/guard_push.py"
GUARD_GOVERNANCE = ".claude/hooks/guard_governance.py"
GUARD_UNATTENDED = ".claude/hooks/guard_unattended.py"

#: MAX_HOURS x BUFFER, in the constants of §3.1.
MAX_WINDOW = timedelta(hours=72 * 1.2)


def git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
        env=clean_env(),
    ).stdout


def clean_env() -> dict[str, str]:
    """The test's environment without GIT_* variables, which would redirect git."""
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


def make_repo(root: Path) -> Path:
    """A git repository with the copied scripts and two plain files, in one commit.

    `CLAUDE.md` is a governed file and `notes.txt` is not, so a suite can touch
    one of each; everything starts committed, so the working tree starts clean.
    """
    root.mkdir(parents=True, exist_ok=True)
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    for relative in COPIED:
        source = REAL / relative
        if source.exists():
            (root / relative).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, root / relative)
    (root / "CLAUDE.md").write_text("# rules\n")
    (root / "notes.txt").write_text("notes\n")
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "root")
    return root


def add_worktree(main: Path, path: Path, branch: str) -> Path:
    """A linked worktree of `main` at `path`, on a new branch (h8 §3.1)."""
    git(main, "worktree", "add", "-q", "-b", branch, str(path))
    return path


def commit_at(repo: Path, moment: datetime, message: str) -> str:
    """An empty commit with committer and author date `moment`; its short hash."""
    stamp = moment.isoformat(timespec="seconds")
    env = {**clean_env(), "GIT_COMMITTER_DATE": stamp, "GIT_AUTHOR_DATE": stamp}
    subprocess.run(
        ["git", "-C", str(repo), "commit", "-q", "--allow-empty", "-m", message],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    return git(repo, "log", "-1", "--format=%h").strip()


def state_dir(repo: Path) -> Path:
    return repo / ".git" / "harness"


def flag_path(repo: Path) -> Path:
    return state_dir(repo) / "unattended.json"


def queue_path(repo: Path) -> Path:
    return state_dir(repo) / "queue.jsonl"


def iso(moment: datetime) -> str:
    return moment.isoformat(timespec="seconds")


def write_flag(repo: Path, content: str | dict[str, Any]) -> Path:
    """Write the flag file as given; a dict is serialised, a str written verbatim."""
    path = flag_path(repo)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content if isinstance(content, str) else json.dumps(content))
    return path


def flag_window(
    repo: Path, since: datetime, until: datetime, keep_awake_pid: int | None = None
) -> Path:
    return write_flag(
        repo,
        {
            "since": iso(since),
            "until": iso(until),
            "set_by": "away",
            "keep_awake_pid": keep_awake_pid,
        },
    )


def flag_on(repo: Path, now: datetime | None = None, keep_awake_pid: int | None = None) -> datetime:
    """The flag §4 calls `on`: until = now + 1 h. Returns `until`."""
    now = now or datetime.now(UTC).replace(microsecond=0)
    until = now + timedelta(hours=1)
    flag_window(repo, now, until, keep_awake_pid)
    return until


def flag_expired(repo: Path, now: datetime | None = None) -> datetime:
    """The flag §4 calls `expired`: since = now - 2 h, until = now - 1 s."""
    now = now or datetime.now(UTC).replace(microsecond=0)
    until = now - timedelta(seconds=1)
    flag_window(repo, now - timedelta(hours=2), until)
    return until


def flag_broken(repo: Path, kind: str, now: datetime | None = None) -> None:
    """Each `broken` flag of §4: not JSON, a directory, naive time, reversed, too long."""
    now = now or datetime.now(UTC).replace(microsecond=0)
    if kind == "not-json":
        write_flag(repo, "{not json")
    elif kind == "directory":
        flag_path(repo).mkdir(parents=True)
    elif kind == "naive":
        naive = now.replace(tzinfo=None)
        write_flag(repo, {"since": iso(naive), "until": iso(naive + timedelta(hours=1))})
    elif kind == "reversed":
        flag_window(repo, now, now)
    elif kind == "too-long":
        flag_window(repo, now, now + MAX_WINDOW + timedelta(seconds=61))
    else:  # pragma: no cover - a typo in a test, not a case
        raise ValueError(kind)


BROKEN_KINDS = ("not-json", "directory", "naive", "reversed", "too-long")


def set_mode(repo: Path, mode: str) -> None:
    """Put the repository's flag in `off`, `on`, `expired` or `broken` (not JSON)."""
    if mode == "on":
        flag_on(repo)
    elif mode == "expired":
        flag_expired(repo)
    elif mode == "broken":
        flag_broken(repo, "not-json")
    elif mode != "off":  # pragma: no cover
        raise ValueError(mode)


def queue_lines(repo: Path) -> list[dict[str, Any]]:
    path = queue_path(repo)
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def run_script(
    repo: Path, relative: str, stdin: str | dict[str, Any], *args: str
) -> subprocess.CompletedProcess[str]:
    """Run a copied script as the harness does: its path, the event on stdin, cwd = repo."""
    script = repo / relative
    assert script.exists(), f"{relative} is missing from the copy"
    text = stdin if isinstance(stdin, str) else json.dumps(stdin)
    return subprocess.run(
        [sys.executable, str(script), *args],
        input=text,
        capture_output=True,
        text=True,
        cwd=repo,
        env=clean_env(),
        timeout=60,
        check=False,
    )


def bash_event(repo: Path, command: str, **extra: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": "Bash",
        "tool_input": {"command": command},
        "cwd": str(repo),
        **extra,
    }


def file_event(repo: Path, tool: str, path: str, **extra: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": tool,
        "tool_input": {"file_path": path, "content": "x"},
        "cwd": str(repo),
        **extra,
    }


SUBAGENT = {"agent_id": "a1b2c3", "agent_type": "tester"}


def pretool_decision(result: subprocess.CompletedProcess[str]) -> tuple[str, str] | None:
    """(permissionDecision, reason) from a PreToolUse hook's stdout, or None if silent.

    Claude Code reads a hook's JSON only on exit 0, so any output that is to take
    effect must come with exit status 0.
    """
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    if not result.stdout.strip():
        return None
    output = json.loads(result.stdout)["hookSpecificOutput"]
    assert output["hookEventName"] == "PreToolUse"
    return output["permissionDecision"], output["permissionDecisionReason"]


def load_tool(name: str) -> ModuleType:
    """Import `tools/<name>.py` from the real checkout, or fail naming what is missing."""
    path = REAL / "tools" / f"{name}.py"
    if not path.exists():
        pytest.fail(f"tools/{name}.py is missing from the checkout")
    spec = importlib.util.spec_from_file_location(f"u1_{name}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # dataclasses resolve annotations through it
    spec.loader.exec_module(module)
    return module


class Tool:
    """`tools/<name>.py`, loaded on first attribute access.

    Loading lazily, inside the test body, makes a missing module a test failure
    naming the file rather than a fixture error, and keeps collection working.
    """

    def __init__(self, name: str) -> None:
        self._name = name
        self._module: ModuleType | None = None

    def __getattr__(self, attr: str) -> Any:
        # pytest probes module-level objects during collection (`__test__`,
        # `pytestmark`, ...); those probes must not load the module.
        if attr.startswith(("_", "pytest")):
            raise AttributeError(attr)
        if self._module is None:
            self._module = load_tool(self._name)
        return getattr(self._module, attr)
