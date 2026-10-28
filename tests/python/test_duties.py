"""`tools/duties.py`: the `@orchestrator` runs owed, computed from git and the away log.

Spec: docs/increments/h23-duties.md. Section 7 lists these tests by number;
each block below carries its number. docs/increments/h34-lean-process.md
replaces the research round by the weekly retrospective and adds the lessons
log (`logged`, `lessons_path`), and owes a post-merge check only for a merge
with a `docs/increments/` path (a full-lane merge); its section 4 tests are the
`h34` blocks at the end, and the older tests it amends say so. The pure
functions (`window_ends`, `owed`, `render`) get synthetic inputs; `gather`,
`main` and the recap run
against git repositories and away logs built under `tmp_path`, never the
real away state. `tools/` is not a package, so the module is loaded from its
path, inside each test, so a missing module fails the test naming the file.
"""

from __future__ import annotations

import io
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import UTC, datetime, timedelta, timezone
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import harness_fixtures
from harness_fixtures import Tool

pytestmark = pytest.mark.harness

duties = Tool("duties")

#: The section 3 fault line: the exception's name, then its message.
UNKNOWN = re.compile(r"^Duties: unknown \(\w+: .*\)$")


def utc(day: str, hm: str = "00:00", s: int = 0) -> datetime:
    """`utc("2026-10-17", "09:40")` is 2026-10-17 09:40 UTC."""
    return datetime.fromisoformat(f"{day}T{hm}:{s:02d}+00:00")


def inputs(**fields: Any) -> Any:
    """`duties.Inputs` with nothing owed, overridden by `fields`."""
    base: dict[str, Any] = {
        "retro": utc("2026-10-17", "06:40"),
        "logged": None,
        "merges": [],
        "windows": [],
        "unmerged": [],
    }
    return duties.Inputs(**{**base, **fields})


def merge(when: datetime, label: str, *paths: str) -> Any:
    return duties.Merge(when=when, label=label, paths=tuple(paths))


NOW = utc("2026-10-17", "12:00")


# ---------------------------------------------------------------- test 1, window_ends
#
# Shaped on the real log's lines 39-40 (two enters, the first ends at its own
# until, before the next since), 46-47 (two enters, the first cut off by the
# next since) and 81-83 (an enter closed by its back, then one still running).

EVENTS: list[dict[str, Any]] = [
    # Malformed, first, so skipping them cannot change any neighbour's end.
    {"event": "enter", "since": "not a time", "until": "2026-10-06T10:00:00+00:00"},
    {"event": "enter"},
    {"event": "back"},
    {"event": "sideways", "since": "2026-10-06T08:00:00+00:00",
     "until": "2026-10-06T09:00:00+00:00"},
    # Lines 39-40.
    {"event": "enter", "since": "2026-10-06T12:15:00+00:00", "until": "2026-10-06T13:27:00+00:00"},
    {"event": "enter", "since": "2026-10-06T22:27:04+00:00", "until": "2026-10-07T08:03:04+00:00"},
    # Lines 46-47.
    {"event": "enter", "since": "2026-10-08T16:06:15+00:00", "until": "2026-10-09T01:42:15+00:00"},
    {"event": "enter", "since": "2026-10-08T22:26:14+00:00", "until": "2026-10-09T08:02:14+00:00"},
    # Lines 81-83.
    {"event": "enter", "since": "2026-10-16T23:24:03+00:00", "until": "2026-10-17T11:24:03+00:00",
     "grants": ["worktree-legend"]},
    {"event": "back", "since": "2026-10-16T23:24:03+00:00", "until": "2026-10-17T06:45:45+00:00"},
    {"event": "enter", "since": "2026-10-17T12:16:28+00:00", "until": "2026-10-18T05:04:28+00:00",
     "grants": ["worktree-retro-1017"]},
]  # fmt: skip


def _window(since: str, end: str) -> Any:
    return duties.Window(start=datetime.fromisoformat(since), end=datetime.fromisoformat(end))


def test_window_ends_reads_back_pairs_cut_offs_and_skips_running_and_malformed() -> None:
    now = utc("2026-10-17", "20:00")
    assert duties.window_ends(EVENTS, now) == [
        # An enter with no back ends at its own until when the next since is later.
        _window("2026-10-06T12:15:00+00:00", "2026-10-06T13:27:00+00:00"),
        _window("2026-10-06T22:27:04+00:00", "2026-10-07T08:03:04+00:00"),
        # ... and at the next enter's since when that is earlier than its until.
        _window("2026-10-08T16:06:15+00:00", "2026-10-08T22:26:14+00:00"),
        _window("2026-10-08T22:26:14+00:00", "2026-10-09T08:02:14+00:00"),
        # A back ends the window it names, at its until.
        _window("2026-10-16T23:24:03+00:00", "2026-10-17T06:45:45+00:00"),
        # The last enter ends after now: still running, so dropped.
    ]


def test_window_ends_of_an_empty_log_is_empty() -> None:
    assert duties.window_ends([], NOW) == []


# ---------------------------------------------------------------- test 2, D2 threshold

R = utc("2026-10-17", "06:40")


def _ended(hours: float, end: datetime) -> Any:
    return duties.Window(start=end - timedelta(hours=hours), end=end)


def test_a_window_of_3_h_59_min_ending_after_r_owes_nothing() -> None:
    end = R + timedelta(hours=1)
    window = duties.Window(start=end - timedelta(hours=3, minutes=59), end=end)
    assert duties.owed(inputs(retro=R, windows=[window]), NOW) == []


def test_a_window_of_4_h_0_min_ending_after_r_owes_a_morning_check() -> None:
    end = utc("2026-10-17", "07:40", 30)
    assert duties.owed(inputs(retro=R, windows=[_ended(4, end)]), NOW) == [
        "morning check (1 away window ended since, latest 2026-10-17 07:40 UTC)"
    ]


def test_a_9_h_window_ending_one_minute_before_r_owes_nothing() -> None:
    window = _ended(9, R - timedelta(minutes=1))
    assert duties.owed(inputs(retro=R, windows=[window]), NOW) == []


def test_a_window_ending_after_now_owes_nothing_yet() -> None:
    window = _ended(9, NOW + timedelta(minutes=1))
    assert duties.owed(inputs(retro=R, windows=[window]), NOW) == []


def test_two_windows_count_and_name_the_latest_end() -> None:
    first, second = _ended(8, R + timedelta(hours=1)), _ended(5, R + timedelta(hours=2))
    assert duties.owed(inputs(retro=R, windows=[first, second]), NOW) == [
        "morning check (2 away windows ended since, latest 2026-10-17 08:40 UTC)"
    ]


# ---------------------------------------------------------------- test 3, D1

CODE = "src_python/tin_engine/x.py"
#: A full-lane merge carries its increment file; only such a merge owes a
#: post-merge check (h34 section 3; a code-only merge is fast-lane, h34 test 9).
INCREMENT = "docs/increments/46-x.md"
FULL = (CODE, INCREMENT)


def test_a_merge_one_second_after_r_owes_a_post_merge_check() -> None:
    one = merge(R + timedelta(seconds=1), "#280", *FULL)
    assert duties.owed(inputs(retro=R, merges=[one]), NOW) == [
        "post-merge check (1 merge since the last log entry, 2026-10-17 06:40 UTC: #280)"
    ]


def test_a_merge_one_second_before_r_owes_nothing() -> None:
    one = merge(R - timedelta(seconds=1), "#279", *FULL)
    assert duties.owed(inputs(retro=R, merges=[one]), NOW) == []


@pytest.mark.parametrize(
    "paths",
    [
        pytest.param(("docs/retrospectives/x.md",), id="retrospective-alone"),
        pytest.param(
            ("ROADMAP.md", "docs/increments/46-x.md", "docs/retrospectives/x.md"), id="shape-of-268"
        ),
    ],
)
def test_a_check_recording_merge_owes_nothing(paths: tuple[str, ...]) -> None:
    one = merge(R + timedelta(hours=1), "#268", *paths)
    assert duties.owed(inputs(retro=R, merges=[one]), NOW) == []


@pytest.mark.parametrize(
    "paths",
    [
        pytest.param(("ROADMAP.md", "docs/increments/46-x.md"), id="shape-of-237"),
        pytest.param(("docs/retrospectives/x.md", *FULL), id="retrospective-and-code"),
    ],
)
def test_a_merge_that_is_not_check_recording_owes(paths: tuple[str, ...]) -> None:
    one = merge(R + timedelta(hours=1), "#237", *paths)
    assert duties.owed(inputs(retro=R, merges=[one]), NOW) == [
        "post-merge check (1 merge since the last log entry, 2026-10-17 06:40 UTC: #237)"
    ]


def test_five_merges_list_the_newest_three_and_the_rest_as_a_count() -> None:
    # Newest first, the order `git log` gives.
    five = [merge(R + timedelta(minutes=10 * n), f"#{n}", *FULL) for n in (5, 4, 3, 2, 1)]
    assert duties.owed(inputs(retro=R, merges=five), NOW) == [
        "post-merge check (5 merges since the last log entry, 2026-10-17 06:40 UTC: #5, #4, #3, +2)"
    ]


def test_the_last_retrospective_time_prints_in_utc_whatever_its_offset() -> None:
    oslo = datetime(2026, 10, 17, 11, 40, tzinfo=timezone(timedelta(hours=2)))
    one = merge(oslo + timedelta(minutes=5), "#280", *FULL)
    assert duties.owed(inputs(retro=oslo, merges=[one]), NOW) == [
        "post-merge check (1 merge since the last log entry, 2026-10-17 09:40 UTC: #280)"
    ]


def test_with_no_baseline_on_record_a_merge_owes_and_says_so() -> None:
    """h34 amends: no retrospective and no log entry; the weekly item is owed too."""
    one = merge(R, "#280", *FULL)
    assert duties.owed(inputs(retro=None, logged=None, merges=[one]), NOW) == [
        "post-merge check (no log entry on record)",
        "weekly retrospective (none on record)",
    ]


# ---------------------------------------------------------------- test 4, D3
#
# Removed by h34: the research round is no longer its own duty; the weekly
# retrospective replaces it (the h34 test 6 block).


# ---------------------------------------------------------------- test 5, render
#
# Section 3's golden line, with one change: its R (09:40) is after its away
# window's end (06:45), which by section 2 would not owe a morning check, so R
# here is 06:40. h34 amends it: R is the last log entry, and with no
# retrospective on record the weekly item takes the research item's place.

GOLDEN = (
    "Duties owed: post-merge check (1 merge since the last log entry, 2026-10-17 06:40 UTC:"
    " #280); morning check (1 away window ended since, latest 2026-10-17 06:45 UTC); weekly"
    " retrospective (none on record); push 6 retrospective branches (worktree-retro-1004c,"
    " worktree-retro-1005, +4)"
)
UNMERGED = [  # unsorted, so the line's order is the tool's
    "worktree-retro-1017",
    "worktree-retro-1005",
    "worktree-retro-1004c",
    "worktree-retro-1015",
    "worktree-retro-1005b",
    "worktree-retro-1008",
]


def _golden_inputs() -> Any:
    return inputs(
        retro=None,
        logged=R,
        merges=[merge(utc("2026-10-17", "10:15"), "#280", *FULL)],
        windows=[_window("2026-10-16T23:24:03+00:00", "2026-10-17T06:45:45+00:00")],
        unmerged=UNMERGED,
    )


@pytest.mark.parametrize("place", ["recap", "status"])
def test_render_prints_the_golden_line_in_the_recap_and_the_status_line(place: str) -> None:
    assert duties.render(duties.owed(_golden_inputs(), NOW), place) == GOLDEN


def test_render_in_the_hook_leaves_out_the_push_item() -> None:
    hook = GOLDEN.removesuffix(
        "; push 6 retrospective branches (worktree-retro-1004c, worktree-retro-1005, +4)"
    )
    assert duties.render(duties.owed(_golden_inputs(), NOW), "hook") == hook


def test_render_of_the_push_item_alone_is_silent_in_the_hook_only() -> None:
    items = duties.owed(inputs(unmerged=["worktree-retro-1005"]), NOW)
    assert items == ["push 1 retrospective branch (worktree-retro-1005)"]
    assert duties.render(items, "hook") is None
    line = "Duties owed: push 1 retrospective branch (worktree-retro-1005)"
    assert duties.render(items, "recap") == line
    assert duties.render(items, "status") == line


def test_two_unmerged_branches_list_both_and_no_count() -> None:
    items = duties.owed(inputs(unmerged=["b-retro", "a-retro"]), NOW)
    assert items == ["push 2 retrospective branches (a-retro, b-retro)"]


@pytest.mark.parametrize(
    ("place", "line"), [("hook", None), ("recap", None), ("status", "Duties: none owed")]
)
def test_render_with_nothing_owed(place: str, line: str | None) -> None:
    items = duties.owed(inputs(), NOW)
    assert items == []
    assert duties.render(items, place) == line


# ---------------------------------------------------------------- test 6, gather
#
# One repository, grown in four stages; after each, `gather` is read again.
# Commit and merge dates are fixed; `refs/remotes/origin/master` is set by
# `update-ref`, as a fetch would leave it.

T0 = utc("2026-10-01", "08:00")
T1 = utc("2026-10-02", "09:00")  # retrospective first added, branch retro-a
T2 = utc("2026-10-03", "09:00")  # the same path re-added, branch retro-b
TM1 = utc("2026-10-03", "11:00")  # code merge #11
TM1B = utc("2026-10-03", "12:00")  # code merge with no PR number
TM2 = utc("2026-10-04", "09:00")  # retrospective-only merge #12 (retro-b)
T3 = utc("2026-10-05", "09:00")  # docs/retrospectives/next.md, branch agenda
T4 = utc("2026-10-06", "09:00")  # a check only on origin/master
RETRO = "docs/retrospectives/2026-10-02-a.md"
#: The code merges are full-lane: each carries its increment file (h34 section 3).
INC11 = {"docs/increments/11-code.md": "# 11\n"}
INC12 = {"docs/increments/12-other.md": "# 12\n"}


def git(repo: Path, *args: str) -> str:
    return harness_fixtures.git(repo, *args)


def _dated(moment: datetime) -> dict[str, str]:
    stamp = moment.isoformat(timespec="seconds")
    return {**harness_fixtures.clean_env(), "GIT_COMMITTER_DATE": stamp, "GIT_AUTHOR_DATE": stamp}


def commit(repo: Path, moment: datetime, subject: str, files: dict[str, str], *extra: str) -> str:
    """Write `files`, commit them at `moment` with `subject`; the full sha."""
    for relative, text in files.items():
        (repo / relative).parent.mkdir(parents=True, exist_ok=True)
        (repo / relative).write_text(text)
        git(repo, "add", relative)
    subprocess.run(
        ["git", "-C", str(repo), "commit", "-q", *extra, "-m", subject],
        check=True, capture_output=True, text=True, env=_dated(moment),
    )  # fmt: skip
    return git(repo, "rev-parse", "HEAD").strip()


def merge_into_master(repo: Path, moment: datetime, branch: str, subject: str) -> str:
    git(repo, "checkout", "-q", "master")
    subprocess.run(
        ["git", "-C", str(repo), "merge", "-q", "--no-ff", "-m", subject, branch],
        check=True, capture_output=True, text=True, env=_dated(moment),
    )  # fmt: skip
    return git(repo, "rev-parse", "--short", "HEAD").strip()


def fetched(repo: Path, rev: str = "master") -> None:
    git(repo, "update-ref", "refs/remotes/origin/master", rev)


def new_repo(root: Path) -> Path:
    root.mkdir()
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    git(root, "config", "commit.gpgsign", "false")
    commit(root, T0, "root", {"ROADMAP.md": "# roadmap\n", "src/code.py": "x = 0\n"})
    return root


@pytest.fixture(autouse=True)
def _no_outer_git(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """`gather` runs git in this process: no GIT_* redirection, no repository above tmp_path."""
    for name in [n for n in os.environ if n.startswith("GIT_")]:
        monkeypatch.delenv(name)
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))


def _post_merge(items: list[str]) -> list[str]:
    return [item for item in items if item.startswith("post-merge check")]


def test_gather_reads_retrospective_times_merges_and_unmerged_branches(tmp_path: Path) -> None:
    repo = new_repo(tmp_path / "g")
    root = git(repo, "rev-parse", "HEAD").strip()
    # Stage 1: a check on retro-a at T1, carried onto retro-b at T2; two code merges.
    git(repo, "checkout", "-q", "-b", "retro-a", root)
    commit(repo, T1, "Retrospective a (@orchestrator)", {RETRO: "# a\n"})
    git(repo, "checkout", "-q", "-b", "retro-b", root)
    commit(repo, T2, "Retrospective a, carried (@orchestrator)", {RETRO: "# a\n"})
    git(repo, "checkout", "-q", "-b", "code", root)
    commit(repo, TM1 - timedelta(hours=1), "Code (@developer)", {"src/code.py": "x = 1\n", **INC11})
    merge_into_master(repo, TM1, "code", "Merge pull request #11 from x/code")
    git(repo, "checkout", "-q", "-b", "feat2", "master")
    commit(
        repo,
        TM1B - timedelta(minutes=30),
        "Other (@developer)",
        {"src/other.py": "y = 1\n", **INC12},
    )
    bare = merge_into_master(repo, TM1B, "feat2", "Merge branch 'feat2'")
    fetched(repo)

    stage1 = duties.gather(repo, True)
    assert stage1.retro == T1  # the earliest add of the path, not the carried one
    assert stage1.logged is None  # no lessons log beside the repository
    assert sorted(stage1.unmerged) == ["retro-a", "retro-b"]
    d1 = "post-merge check (2 merges since the last log entry, 2026-10-02 09:00 UTC: "
    assert _post_merge(duties.owed(stage1, TM2)) == [f"{d1}{bare}, #11)"]
    assert list(duties.gather(repo, False).unmerged) == []

    # Stage 2: retro-b merged by a retrospective-only merge, then fetched.
    merge_into_master(repo, TM2, "retro-b", "Merge pull request #12 from x/retro-b")
    fetched(repo)
    stage2 = duties.gather(repo, True)
    assert stage2.retro == T1
    assert list(stage2.unmerged) == []  # retro-a's path is on origin/master now
    assert _post_merge(duties.owed(stage2, T3)) == [f"{d1}{bare}, #11)"]

    # Stage 3: the agenda file, by a commit that is not a check, leaves R alone ...
    git(repo, "checkout", "-q", "-b", "agenda", "master")
    commit(repo, T3, "Retrospective agenda (Ola)", {"docs/retrospectives/next.md": "# next\n"})
    assert duties.gather(repo, True).retro == T1
    # ... and the same commit as a check moves it.
    commit(repo, T3, "Retrospective agenda (@orchestrator)", {}, "--amend")
    stage3 = duties.gather(repo, True)
    assert stage3.retro == T3
    assert _post_merge(duties.owed(stage3, T4)) == []

    # Stage 4: a check that only origin/master holds, as after a merge whose
    # local branch was deleted.
    git(repo, "checkout", "-q", "--detach", "master")
    weekly = "docs/retrospectives/2026-10-06-weekly.md"
    only_remote = commit(repo, T4, "Weekly retrospective (@orchestrator)", {weekly: "# r\n"})
    git(repo, "checkout", "-q", "master")
    fetched(repo, only_remote)
    assert git(repo, "branch", "--contains", only_remote).strip() == ""
    stage4 = duties.gather(repo, True)
    assert stage4.retro == T4


def test_gather_without_origin_master_reads_master(tmp_path: Path) -> None:
    repo = new_repo(tmp_path / "g")
    commit(repo, T1, "Retrospective a (@orchestrator)", {RETRO: "# a\n"})
    git(repo, "checkout", "-q", "-b", "code")
    commit(repo, TM1 - timedelta(hours=1), "Code (@developer)", {"src/code.py": "x = 1\n", **INC11})
    merge_into_master(repo, TM1, "code", "Merge pull request #11 from x/code")

    found = duties.gather(repo, False)
    assert found.retro == T1
    assert _post_merge(duties.owed(found, TM2)) == [
        "post-merge check (1 merge since the last log entry, 2026-10-02 09:00 UTC: #11)"
    ]


def test_gather_reads_the_away_log_and_skips_a_line_that_is_not_json(tmp_path: Path) -> None:
    repo = new_repo(tmp_path / "g")
    fetched(repo)
    log = harness_fixtures.state_dir(repo) / "windows.jsonl"
    log.parent.mkdir(parents=True)
    lines = [json.dumps(event) for event in EVENTS[-3:-1]]  # an enter and its back
    log.write_text("\n".join(["{not json", *lines]) + "\n")
    assert list(duties.gather(repo, False).windows) == [
        _window("2026-10-16T23:24:03+00:00", "2026-10-17T06:45:45+00:00")
    ]


def test_gather_with_no_away_log_has_no_windows(tmp_path: Path) -> None:
    repo = new_repo(tmp_path / "g")
    fetched(repo)
    assert not (harness_fixtures.state_dir(repo) / "windows.jsonl").exists()
    assert list(duties.gather(repo, False).windows) == []


# ---------------------------------------------------------------- test 7, main


def _module() -> ModuleType:
    return harness_fixtures.load_tool("duties")


def _status_stdin(directory: Path) -> io.StringIO:
    return io.StringIO(json.dumps({"workspace": {"project_dir": str(directory)}}))


def test_main_outside_a_repository_prints_the_unknown_line_and_returns_0(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    plain = tmp_path / "plain"
    plain.mkdir()
    assert _module().main(["--status"], _status_stdin(plain)) == 0
    out = capsys.readouterr().out.splitlines()
    assert len(out) == 1 and UNKNOWN.match(out[0]), out


def test_main_with_status_stdin_that_is_not_json_prints_the_unknown_line(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert _module().main(["--status"], io.StringIO("{not json")) == 0
    out = capsys.readouterr().out.splitlines()
    assert len(out) == 1 and UNKNOWN.match(out[0]), out


def test_main_cuts_the_fault_message_at_120_characters(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    module = _module()

    def broken(repo: Path, unmerged: bool) -> Any:
        raise RuntimeError("x" * 500)

    monkeypatch.setattr(module, "gather", broken)
    assert module.main(["--status"], _status_stdin(tmp_path)) == 0
    assert capsys.readouterr().out == f"Duties: unknown (RuntimeError: {'x' * 120})\n"


def _run_copy(root: Path, stdin: str, *args: str) -> subprocess.CompletedProcess[str]:
    """`tools/duties.py` as the hook and the status line run it, from inside `root`."""
    env = {**harness_fixtures.clean_env(), "CLAUDE_PROJECT_DIR": str(root)}
    env["GIT_CEILING_DIRECTORIES"] = str(root.parent)
    return subprocess.run(
        [sys.executable, str(root / "tools" / "duties.py"), *args],
        input=stdin, capture_output=True, text=True, cwd=root, env=env, timeout=60, check=False,
    )  # fmt: skip


def test_the_hook_outside_a_repository_prints_the_unknown_line_and_exits_0(tmp_path: Path) -> None:
    plain = tmp_path / "plain"
    (plain / "tools").mkdir(parents=True)
    for name in ("duties.py", "harness_mode.py"):
        source = harness_fixtures.REAL / "tools" / name
        assert source.exists(), f"tools/{name} is missing from the checkout"
        shutil.copy2(source, plain / "tools" / name)
    result = _run_copy(plain, "")
    assert result.returncode == 0, result.stderr
    out = result.stdout.splitlines()
    assert len(out) == 1 and UNKNOWN.match(out[0]), out


# ---------------------------------------------------------------- test 8, the recap


@pytest.fixture
def harness_repo(tmp_path: Path) -> Path:
    repo = harness_fixtures.make_repo(tmp_path / "repo")
    assert (repo / "tools" / "duties.py").exists(), "tools/duties.py is missing from the copy"
    return repo


def _check_now(repo: Path, name: str, branch: str | None = None) -> None:
    """A check commit, dated now, adding `docs/retrospectives/<name>`."""
    if branch:
        git(repo, "checkout", "-q", "-b", branch)
    commit(repo, datetime.now(UTC), f"Retrospective {name} (@orchestrator)",
           {f"docs/retrospectives/{name}": "# r\n"})  # fmt: skip
    if branch:
        git(repo, "checkout", "-q", "master")


def _recap(repo: Path) -> list[str]:
    result = harness_fixtures.run_script(repo, "tools/session_state.py", "")
    assert result.returncode == 0, result.stderr
    return result.stdout.splitlines()


def _after_last_landed(lines: list[str]) -> str:
    landed = next(i for i, line in enumerate(lines) if line.startswith("Last landed:"))
    return lines[landed + 1]


def test_the_recap_prints_the_duties_line_after_last_landed(harness_repo: Path) -> None:
    """h34 amends: no retrospective on record owes the weekly one."""
    fetched(harness_repo)
    line = _after_last_landed(_recap(harness_repo))
    assert line == "Duties owed: weekly retrospective (none on record)"


def test_the_recap_prints_no_duties_line_when_none_is_owed(harness_repo: Path) -> None:
    """h34 amends: any retrospective, not only a research one, ends the weekly duty."""
    _check_now(harness_repo, "2026-10-20-plain.md")
    fetched(harness_repo)
    lines = _recap(harness_repo)
    assert not any(line.startswith("Duties") for line in lines), lines
    assert _after_last_landed(lines) == "In flight:"


def test_the_recap_with_a_broken_duties_module_prints_one_line_and_the_rest(
    harness_repo: Path,
) -> None:
    (harness_repo / "tools" / "duties.py").write_text('raise RuntimeError("planted")\n')
    lines = _recap(harness_repo)
    assert _after_last_landed(lines) == "Duties: unknown (RuntimeError: planted)"
    assert len([line for line in lines if line.startswith("Duties")]) == 1
    assert "Next on ROADMAP.md:" in lines
    assert "== .claude/current-task/ ==" in lines


# ---------------------------------------------------------------- the hook leaves out D4


def test_the_hook_is_silent_on_a_push_item_the_status_line_shows(harness_repo: Path) -> None:
    _check_now(harness_repo, "2026-10-20-plain.md")
    fetched(harness_repo)
    _check_now(harness_repo, "2026-10-20-unpushed.md", branch="retro-unpushed")
    hook = _run_copy(harness_repo, "")
    assert (hook.returncode, hook.stdout) == (0, ""), hook.stderr
    status = _run_copy(harness_repo, _status_stdin(harness_repo).getvalue(), "--status")
    assert status.returncode == 0, status.stderr
    assert status.stdout == "Duties owed: push 1 retrospective branch (retro-unpushed)\n"


# ================================================================ h34, the lessons log
#
# docs/increments/h34-lean-process.md section 4, tests 1 to 7 (test 8 is in
# test_brief.py). Fixture times are UTC.

NOW34 = utc("2026-10-27", "18:00")
LOG = (
    "# lessons\n"
    "## 2026-10-26T09:00Z merge #380\n"
    "- a finding, with its evidence\n"
    "## 2026-10-27T12:00Z night\n"
    "## 2026-10-27T13:00 merge #381\n"  # no `Z`: not a heading the log counts
    "## not a heading\n"
)


# ---------------------------------------------------------------- h34 test 1, logged


def test_h34_1_logged_is_the_latest_well_formed_stamp() -> None:
    found = duties.logged(LOG, NOW34)
    assert found == utc("2026-10-27", "12:00")
    assert found is not None and found.utcoffset() == timedelta(0)


def test_h34_1_logged_of_empty_text_is_none() -> None:
    assert duties.logged("", NOW34) is None


@pytest.mark.parametrize(
    "heading",
    [
        pytest.param("## 2026-10-27T14:00Z nightly", id="night-then-a-letter"),
        pytest.param("## 2026-10-27T14:00Z merge 382", id="merge-without-hash"),
        pytest.param("## 2026-10-27T14:00Z review", id="other-kind"),
        pytest.param(" ## 2026-10-27T14:00Z night", id="indented"),
        pytest.param("### 2026-10-27T14:00Z night", id="deeper-heading"),
    ],
)
def test_h34_1_a_heading_off_the_pattern_is_ignored(heading: str) -> None:
    assert duties.logged(f"{LOG}{heading}\n", NOW34) == utc("2026-10-27", "12:00")


# ---------------------------------------------------------------- h34 test 2, future stamps


def test_h34_2_a_stamp_later_than_now_is_ignored() -> None:
    later = "## 2026-10-28T09:00Z night\n"  # 15 hours after NOW34
    assert duties.logged(LOG + later, NOW34) == utc("2026-10-27", "12:00")
    assert duties.logged(later, NOW34) is None


def test_h34_2_a_stamp_equal_to_now_counts() -> None:
    """ "Not later than `now`": the boundary is kept."""
    assert duties.logged("## 2026-10-27T18:00Z merge #383\n", NOW34) == NOW34


# ---------------------------------------------------------------- h34 tests 3 and 4, post-merge

MERGES34 = [  # newest first, as `git log` gives them
    merge(utc("2026-10-27", "13:00"), "#382", *FULL),
    merge(utc("2026-10-27", "11:00"), "#381", *FULL),
]
AFTER_NOON = "post-merge check (1 merge since the last log entry, 2026-10-27 12:00 UTC: #382)"
RECENT34 = utc("2026-10-21", "10:00")  # a retrospective under a week old at NOW34


def test_h34_3_merges_after_the_log_entry_are_owed() -> None:
    """A week-old retrospective: the log entry is the baseline, and the weekly item is owed."""
    found = inputs(retro=utc("2026-10-20", "10:00"), logged=utc("2026-10-27", "12:00"),
                   merges=MERGES34)  # fmt: skip
    assert duties.owed(found, NOW34) == [AFTER_NOON, "weekly retrospective (last 2026-10-20)"]


def test_h34_4_with_no_log_entry_the_retrospective_is_the_baseline() -> None:
    found = inputs(retro=utc("2026-10-27", "12:00"), logged=None, merges=MERGES34)
    assert duties.owed(found, NOW34) == [AFTER_NOON]


def test_h34_4_the_later_of_the_two_is_the_baseline() -> None:
    found = inputs(retro=utc("2026-10-27", "12:00"), logged=utc("2026-10-27", "10:00"),
                   merges=MERGES34)  # fmt: skip
    assert duties.owed(found, NOW34) == [AFTER_NOON]


def test_h34_4_with_neither_no_log_entry_is_on_record() -> None:
    found = inputs(retro=None, logged=None, merges=MERGES34)
    assert duties.owed(found, NOW34) == [
        "post-merge check (no log entry on record)",
        "weekly retrospective (none on record)",
    ]


# ---------------------------------------------------------------- h34 test 9, fast lane


@pytest.mark.parametrize(
    "paths",
    [
        pytest.param((CODE,), id="python-code"),
        pytest.param(("include/terrain/cdt/x.hpp", "tests/cpp/x.cpp"), id="core-code-and-tests"),
        pytest.param(("tools/duties.py", "ROADMAP.md"), id="tool-and-roadmap"),
        pytest.param(("docs/increments.md", "docs/x.md"), id="docs-but-not-increments"),
    ],
)
def test_h34_9_a_merge_with_no_increment_path_owes_no_post_merge_check(
    paths: tuple[str, ...],
) -> None:
    """A fast-lane merge, after the baseline, owes nothing; the full-lane one
    beside it still does, alone."""
    fast = merge(utc("2026-10-27", "14:00"), "#383", *paths)
    found = inputs(retro=RECENT34, logged=utc("2026-10-27", "12:00"),
                   merges=[fast, *MERGES34])  # fmt: skip
    assert duties.owed(found, NOW34) == [AFTER_NOON]
    alone = inputs(retro=RECENT34, logged=utc("2026-10-27", "12:00"), merges=[fast])
    assert duties.owed(alone, NOW34) == []


def test_h34_9_with_no_baseline_a_fast_lane_merge_alone_owes_no_post_merge_check() -> None:
    found = inputs(retro=None, logged=None, merges=[merge(NOW34, "#383", CODE)])
    assert duties.owed(found, NOW34) == ["weekly retrospective (none on record)"]


# ---------------------------------------------------------------- h34 test 5, morning

RECENT = utc("2026-10-21", "10:00")  # a retrospective under a week old at NOW34


def test_h34_5_a_window_ending_before_the_log_entry_owes_nothing() -> None:
    window = _ended(5, utc("2026-10-27", "11:00"))
    found = inputs(retro=RECENT, logged=utc("2026-10-27", "12:00"), windows=[window])
    assert duties.owed(found, NOW34) == []


def test_h34_5_a_window_ending_after_the_log_entry_owes_the_morning_check() -> None:
    window = _ended(5, utc("2026-10-27", "13:00"))
    found = inputs(retro=RECENT, logged=utc("2026-10-27", "12:00"), windows=[window])
    assert duties.owed(found, NOW34) == [
        "morning check (1 away window ended since, latest 2026-10-27 13:00 UTC)"
    ]


# ---------------------------------------------------------------- h34 test 6, weekly


@pytest.mark.parametrize(
    ("retro", "items"),
    [
        pytest.param(utc("2026-10-19", "10:00"), ["weekly retrospective (last 2026-10-19)"],
                     id="8-days"),
        pytest.param(NOW34 - timedelta(days=7), ["weekly retrospective (last 2026-10-20)"],
                     id="7-days-exactly"),
        pytest.param(NOW34 - timedelta(days=7) + timedelta(seconds=1), [], id="7-days-less-1-s"),
        pytest.param(RECENT, [], id="6-days"),
        pytest.param(None, ["weekly retrospective (none on record)"], id="none"),
    ],
)  # fmt: skip
def test_h34_6_the_weekly_retrospective_is_owed_a_week_after_the_last(
    retro: datetime | None, items: list[str]
) -> None:
    assert duties.owed(inputs(retro=retro), NOW34) == items


def test_h34_6_no_item_names_a_research_round() -> None:
    window = _ended(5, utc("2026-10-27", "13:00"))
    found = inputs(retro=None, merges=MERGES34, windows=[window], unmerged=["retro-x"])
    items = duties.owed(found, NOW34)
    assert len(items) == 4, items
    assert not [item for item in items if "research" in item], items


def test_h34_6_inputs_has_no_research_field() -> None:
    assert "research" not in duties.Inputs.__dataclass_fields__
    assert "logged" in duties.Inputs.__dataclass_fields__


# ---------------------------------------------------------------- h34 test 7, lessons_path


def _log_beside(main: Path) -> Path:
    """`<parent of the main checkout>/rasputin_scratch/lessons.md`, resolved."""
    return (main.parent / "rasputin_scratch" / "lessons.md").resolve()


def test_h34_7_lessons_path_is_the_same_from_the_repository_and_a_worktree(
    tmp_path: Path,
) -> None:
    (tmp_path / "checkouts").mkdir()
    main = new_repo(tmp_path / "checkouts" / "rasputin")
    tree = harness_fixtures.add_worktree(main, tmp_path / "elsewhere" / "wt", "wt")
    assert duties.lessons_path(main).resolve() == _log_beside(main)
    assert duties.lessons_path(tree).resolve() == _log_beside(main)
    assert duties.lessons_path(tree / "src").resolve() == _log_beside(main)


def test_h34_7_gather_reads_the_log_from_a_worktree_and_starts_the_merge_query_there(
    tmp_path: Path,
) -> None:
    """The log beside the main checkout is read from a worktree; merges at or
    before its entry are not queried (`git log --since` uses the baseline)."""
    (tmp_path / "checkouts").mkdir()
    main = new_repo(tmp_path / "checkouts" / "rasputin")
    git(main, "checkout", "-q", "-b", "code")
    commit(main, TM1 - timedelta(hours=1), "Code (@developer)", {"src/code.py": "x = 1\n", **INC11})
    merge_into_master(main, TM1, "code", "Merge pull request #11 from x/code")
    git(main, "checkout", "-q", "-b", "code2")
    commit(main, TM2 - timedelta(hours=1), "More (@developer)", {"src/more.py": "y = 1\n", **INC12})
    merge_into_master(main, TM2, "code2", "Merge pull request #12 from x/code2")
    fetched(main)
    tree = harness_fixtures.add_worktree(main, tmp_path / "elsewhere" / "wt", "wt")
    assert [m.label for m in duties.gather(tree, False).merges] == ["#12", "#11"]

    log = _log_beside(main)
    log.parent.mkdir(parents=True)
    log.write_text(
        "## 2026-10-03T12:00Z merge #11\n"  # between the two merges
        "## 2999-01-01T00:00Z night\n"  # later than now: ignored
    )
    found = duties.gather(tree, False)
    assert found.logged == utc("2026-10-03", "12:00")
    assert [m.label for m in found.merges] == ["#12"]
    assert _post_merge(duties.owed(found, utc("2026-10-05"))) == [
        "post-merge check (1 merge since the last log entry, 2026-10-03 12:00 UTC: #12)"
    ]
