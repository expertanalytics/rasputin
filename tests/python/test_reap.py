"""`tools/reap.py`: which processes it lists, over a fake process table.

Pins the selection rule only; nothing here signals a real process. The two
cases from 2026-10-30: a process left in a worktree whose branch merged (pid
37142, 23 h on #400), and processes left in a removed worktree (89071 to 89074).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

pytestmark = pytest.mark.harness

TOOL = Path(__file__).resolve().parents[2] / "tools" / "reap.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("reap", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reap = _load()

REPO = "/r/rasputin"
WT = f"{REPO}/.claude/worktrees"
HOUR = 3600
TREES = [
    reap.Tree(REPO, merged=False),
    reap.Tree(f"{WT}/live", merged=False),
    reap.Tree(f"{WT}/done", merged=True),
    reap.Tree("/r/rasputin-side", merged=False),
]


def P(pid: int, cwd: str, age: int = 60, ppid: int = 1, command: str = "python3 x.py"):  # noqa: N802
    return reap.Proc(pid=pid, ppid=ppid, age=age, cwd=cwd, command=command)


def picked(procs, me: int = 1000, only: str | None = None) -> dict[int, str]:
    return {p.pid: why for p, why in reap.select(procs, TREES, REPO, me=me, only=only)}


def test_merged_branch_is_listed_whatever_its_age() -> None:
    assert picked([P(37142, f"{WT}/done/src")]) == {37142: "branch merged"}


def test_removed_worktree_is_listed() -> None:
    procs = [P(89071, f"{WT}/gone"), P(89072, f"{WT}/gone/tests (deleted)"),
             P(89073, f"{WT}/live (deleted)")]  # fmt: skip
    assert picked(procs) == dict.fromkeys((89071, 89072, 89073), "worktree removed")


def test_old_project_process_is_listed_and_young_or_foreign_one_is_not() -> None:
    procs = [
        P(1, REPO, age=7 * HOUR),
        P(2, f"{WT}/live/build", age=6 * HOUR + 1),
        P(3, "/r/rasputin-side", age=9 * HOUR),
        P(4, f"{WT}/live", age=6 * HOUR),
        P(5, "/r/rasputin-other", age=99 * HOUR),
        P(6, "/elsewhere", age=99 * HOUR),
    ]
    assert picked(procs) == {1: "over 6 h", 2: "over 6 h", 3: "over 6 h"}


def test_never_lists_itself_its_ancestors_a_claude_session_or_a_terminal() -> None:
    procs = [
        P(1000, f"{WT}/done", ppid=999),
        P(999, f"{WT}/done", ppid=998),
        P(998, REPO, ppid=1, age=30 * HOUR, command="/usr/local/bin/claude --resume"),
        P(997, REPO, age=30 * HOUR, command="claude"),
        P(996, f"{WT}/done", ppid=1000),
        reap.Proc(995, 1, 30 * HOUR, f"{WT}/done", "-zsh", tty="ttys003"),
        reap.Proc(994, 1, 30 * HOUR, f"{WT}/gone", "sleep 9", tty="??"),
    ]
    assert picked(procs) == {996: "branch merged", 994: "worktree removed"}


def test_worktree_option_limits_to_that_tree() -> None:
    procs = [P(1, f"{WT}/done"), P(2, f"{WT}/gone"), P(3, f"{WT}/done-too"),
             P(4, f"{WT}/done/sub", age=9 * HOUR)]  # fmt: skip
    assert picked(procs, only=f"{WT}/done") == {1: "branch merged", 4: "branch merged"}
    assert picked(procs, only=f"{WT}/gone/") == {2: "worktree removed"}


def test_ps_age_and_listing_lines() -> None:
    assert [reap.seconds(t) for t in ("00:07", "12:34", "01:02:03", "4-01:00:00")] == [
        7,
        754,
        3723,
        4 * 86400 + 3600,
    ]
    lines = reap.listing([(P(42, f"{WT}/gone", age=90061, command="sleep 9"), "worktree removed")])
    assert lines == [f"1. pid 42, 1d 01:01, worktree removed, in {WT}/gone: sleep 9"]
    assert reap.listing([]) == ["none"]
