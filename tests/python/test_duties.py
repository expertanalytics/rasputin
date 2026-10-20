"""`tools/duties.py`: the `@orchestrator` runs owed, computed from git and the away log.

Spec: docs/increments/h23-duties.md. Section 7 lists these tests by number;
each block below carries its number. The pure functions (`window_ends`,
`owed`, `render`) get synthetic inputs; `gather`, `main` and the recap run
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
        "research": utc("2026-10-16"),
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


def test_a_merge_one_second_after_r_owes_a_post_merge_check() -> None:
    one = merge(R + timedelta(seconds=1), "#280", CODE)
    assert duties.owed(inputs(retro=R, merges=[one]), NOW) == [
        "post-merge check (1 merge since the last retrospective, 2026-10-17 06:40 UTC: #280)"
    ]


def test_a_merge_one_second_before_r_owes_nothing() -> None:
    one = merge(R - timedelta(seconds=1), "#279", CODE)
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
        pytest.param(("docs/retrospectives/x.md", CODE), id="retrospective-and-code"),
    ],
)
def test_a_merge_that_is_not_check_recording_owes(paths: tuple[str, ...]) -> None:
    one = merge(R + timedelta(hours=1), "#237", *paths)
    assert duties.owed(inputs(retro=R, merges=[one]), NOW) == [
        "post-merge check (1 merge since the last retrospective, 2026-10-17 06:40 UTC: #237)"
    ]


def test_five_merges_list_the_newest_three_and_the_rest_as_a_count() -> None:
    # Newest first, the order `git log` gives.
    five = [merge(R + timedelta(minutes=10 * n), f"#{n}", CODE) for n in (5, 4, 3, 2, 1)]
    assert duties.owed(inputs(retro=R, merges=five), NOW) == [
        "post-merge check (5 merges since the last retrospective, 2026-10-17 06:40 UTC:"
        " #5, #4, #3, +2)"
    ]


def test_the_last_retrospective_time_prints_in_utc_whatever_its_offset() -> None:
    oslo = datetime(2026, 10, 17, 11, 40, tzinfo=timezone(timedelta(hours=2)))
    one = merge(oslo + timedelta(minutes=5), "#280", CODE)
    assert duties.owed(inputs(retro=oslo, merges=[one]), NOW) == [
        "post-merge check (1 merge since the last retrospective, 2026-10-17 09:40 UTC: #280)"
    ]


def test_with_no_retrospective_on_record_a_merge_owes_and_says_so() -> None:
    one = merge(R, "#280", CODE)
    assert duties.owed(inputs(retro=None, merges=[one]), NOW) == [
        "post-merge check (no retrospective on record)"
    ]


# ---------------------------------------------------------------- test 4, D3


def test_no_research_file_owes_a_research_round_none_on_record() -> None:
    assert duties.owed(inputs(research=None), NOW) == ["research round (none on record)"]


def test_a_research_file_6_days_23_h_old_owes_nothing() -> None:
    research = NOW - timedelta(days=6, hours=23)
    assert duties.owed(inputs(research=research), NOW) == []


def test_a_research_file_7_days_1_h_old_owes_and_names_its_date() -> None:
    research = NOW - timedelta(days=7, hours=1)
    assert duties.owed(inputs(research=research), NOW) == ["research round (last 2026-10-10)"]


# ---------------------------------------------------------------- test 5, render
#
# Section 3's golden line, with one change: its R (09:40) is after its away
# window's end (06:45), which by section 2 would not owe a morning check, so R
# here is 06:40. Every other character is section 3's.

GOLDEN = (
    "Duties owed: post-merge check (1 merge since the last retrospective, 2026-10-17 06:40 UTC:"
    " #280); morning check (1 away window ended since, latest 2026-10-17 06:45 UTC); research"
    " round (none on record); push 6 retrospective branches (worktree-retro-1004c,"
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
        retro=R,
        research=None,
        merges=[merge(utc("2026-10-17", "10:15"), "#280", CODE)],
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
    commit(repo, TM1 - timedelta(hours=1), "Code (@developer)", {"src/code.py": "x = 1\n"})
    merge_into_master(repo, TM1, "code", "Merge pull request #11 from x/code")
    git(repo, "checkout", "-q", "-b", "feat2", "master")
    commit(repo, TM1B - timedelta(minutes=30), "Other (@developer)", {"src/other.py": "y = 1\n"})
    bare = merge_into_master(repo, TM1B, "feat2", "Merge branch 'feat2'")
    fetched(repo)

    stage1 = duties.gather(repo, True)
    assert stage1.retro == T1  # the earliest add of the path, not the carried one
    assert stage1.research is None
    assert sorted(stage1.unmerged) == ["retro-a", "retro-b"]
    d1 = "post-merge check (2 merges since the last retrospective, 2026-10-02 09:00 UTC: "
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
    research = "docs/retrospectives/2026-10-06-research-weekly.md"
    only_remote = commit(repo, T4, "Research round (@orchestrator)", {research: "# r\n"})
    git(repo, "checkout", "-q", "master")
    fetched(repo, only_remote)
    assert git(repo, "branch", "--contains", only_remote).strip() == ""
    stage4 = duties.gather(repo, True)
    assert stage4.retro == T4
    assert stage4.research == T4


def test_gather_without_origin_master_reads_master(tmp_path: Path) -> None:
    repo = new_repo(tmp_path / "g")
    commit(repo, T1, "Retrospective a (@orchestrator)", {RETRO: "# a\n"})
    git(repo, "checkout", "-q", "-b", "code")
    commit(repo, TM1 - timedelta(hours=1), "Code (@developer)", {"src/code.py": "x = 1\n"})
    merge_into_master(repo, TM1, "code", "Merge pull request #11 from x/code")

    found = duties.gather(repo, False)
    assert found.retro == T1
    assert _post_merge(duties.owed(found, TM2)) == [
        "post-merge check (1 merge since the last retrospective, 2026-10-02 09:00 UTC: #11)"
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
    _check_now(harness_repo, "2026-10-20-plain.md")  # a check, but not a research round
    fetched(harness_repo)
    line = _after_last_landed(_recap(harness_repo))
    assert line == "Duties owed: research round (none on record)"


def test_the_recap_prints_no_duties_line_when_none_is_owed(harness_repo: Path) -> None:
    _check_now(harness_repo, "2026-10-20-research-x.md")
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
    assert "== .claude/current-task/session.md ==" in lines


# ---------------------------------------------------------------- the hook leaves out D4


def test_the_hook_is_silent_on_a_push_item_the_status_line_shows(harness_repo: Path) -> None:
    _check_now(harness_repo, "2026-10-20-research-x.md")
    fetched(harness_repo)
    _check_now(harness_repo, "2026-10-20-unpushed.md", branch="retro-unpushed")
    hook = _run_copy(harness_repo, "")
    assert (hook.returncode, hook.stdout) == (0, ""), hook.stderr
    status = _run_copy(harness_repo, _status_stdin(harness_repo).getvalue(), "--status")
    assert status.returncode == 0, status.stderr
    assert status.stdout == "Duties owed: push 1 retrospective branch (retro-unpushed)\n"
