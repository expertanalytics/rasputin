"""The recap `tools/session_state.py` opens every round with (retrospective rule 5).

Rule 5 (docs/retrospectives/2026-09-27-increments-14-to-20b.md): every new
round or increment opens with a structured recap and a roadmap reminder,
generated from files on disk. These tests pin the four parts of that recap:
the last thing landed, what is in flight, the decisions waiting on Ola, and the
next three ROADMAP items. `tools/` is not a package, so the module is loaded
from its path.
"""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path
from types import ModuleType

import pytest

TOOL = Path(__file__).resolve().parents[2] / "tools" / "session_state.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("session_state", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


session_state = _load()

ROADMAP = """\
# Roadmap

| # | What it is | Status | Record |
|---|---|---|---|
| 1 | Predicates | shipped (`d0946bc`) | `docs/increments/01-predicates.md` |
| 5d | The corner graze | designed, **unscheduled** | `docs/increments/05d.md` |
| 9 | Gallery output | designed | `docs/increments/09.md` |
| 20 | Quality start | landed with 20b as interim | `docs/increments/20.md` |
| 20c | Soft quality criterion | to design | `docs/increments/20.md` |
| 21 | Mosaic | planned | none |
| 22 | Sizing field | planned | none |

## Not a table row

| 99 | a pipe-bearing line outside the table | designed | x |
"""


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True, text=True
    ).stdout


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A repository with master, one merge into it, and a branch ahead of it."""
    _git(tmp_path, "init", "-q", "-b", "master")
    _git(tmp_path, "config", "user.email", "t@example.invalid")
    _git(tmp_path, "config", "user.name", "T")
    _git(tmp_path, "commit", "-q", "--allow-empty", "-m", "root")
    _git(tmp_path, "checkout", "-q", "-b", "feature")
    _git(tmp_path, "commit", "-q", "--allow-empty", "-m", "feature work")
    _git(tmp_path, "checkout", "-q", "master")
    _git(tmp_path, "merge", "-q", "--no-ff", "feature", "-m", "Merge pull request #1 feature")
    _git(tmp_path, "checkout", "-q", "-b", "next")
    _git(tmp_path, "commit", "-q", "--allow-empty", "-m", "red: next suite")
    _git(tmp_path, "commit", "-q", "--allow-empty", "-m", "green: next impl")
    return tmp_path


def test_roadmap_next_skips_shipped_landed_and_unscheduled_rows() -> None:
    rows = session_state.roadmap_next(ROADMAP)
    assert [row.split(":")[0] for row in rows] == ["9", "20c", "21"]


def test_roadmap_next_keeps_the_status_so_the_reader_can_judge_it() -> None:
    assert session_state.roadmap_next(ROADMAP)[1] == ("20c: Soft quality criterion [to design]")


def test_roadmap_next_honours_n_and_ignores_rows_outside_the_first_table() -> None:
    rows = session_state.roadmap_next(ROADMAP, n=10)
    assert [row.split(":")[0] for row in rows] == ["9", "20c", "21", "22"]


def test_roadmap_next_on_a_file_without_a_table_is_empty() -> None:
    assert session_state.roadmap_next("# Roadmap\n\nNo table yet.\n") == []


def test_pending_decisions_collects_ask_ola_lines_from_every_task_file(
    tmp_path: Path,
) -> None:
    (tmp_path / "session.md").write_text("SESSION: did a thing.\nASK OLA: fix before landing?\n")
    (tmp_path / "tester-101010.md").write_text("ask ola: keep the slow case?\nother\n")
    (tmp_path / "developer-111111.md").write_text("nothing pending\n")
    assert session_state.pending_decisions(tmp_path) == [
        "session.md: ASK OLA: fix before landing?",
        "tester-101010.md: ask ola: keep the slow case?",
    ]


def test_pending_decisions_survives_a_missing_directory(tmp_path: Path) -> None:
    assert session_state.pending_decisions(tmp_path / "absent") == []


def test_last_landed_is_the_newest_merge_reachable_from_head(repo: Path) -> None:
    landed = session_state.last_landed(repo)
    assert landed.endswith("Merge pull request #1 feature")
    assert "green: next impl" not in landed


def test_last_landed_without_a_merge_says_so(tmp_path: Path) -> None:
    _git(tmp_path, "init", "-q", "-b", "master")
    _git(
        tmp_path,
        "-c",
        "user.email=t@e.invalid",
        "-c",
        "user.name=T",
        "commit",
        "-q",
        "--allow-empty",
        "-m",
        "root",
    )
    assert session_state.last_landed(tmp_path) == "(no merge commit reachable from HEAD)"


def test_in_flight_names_the_branch_and_its_commits_ahead_of_base(repo: Path) -> None:
    lines = session_state.in_flight(repo, base="master")
    assert lines[0] == "branch next, 2 commit(s) ahead of master"
    assert [line.split(" ", 1)[1] for line in lines[1:]] == [
        "green: next impl",
        "red: next suite",
    ]


def test_in_flight_with_an_unknown_base_reports_rather_than_raises(repo: Path) -> None:
    lines = session_state.in_flight(repo, base="no-such-branch")
    assert lines == ["branch next (base no-such-branch not found)"]


def test_default_base_prefers_the_remote_branch_when_it_exists(repo: Path) -> None:
    # A local master lags the remote one after a PR merged on GitHub; counting
    # "ahead of master" against it lists already-merged commits as in flight.
    _git(repo, "update-ref", "refs/remotes/origin/master", "master")
    assert session_state.default_base(repo) == "origin/master"


def test_default_base_falls_back_to_local_master(repo: Path) -> None:
    assert session_state.default_base(repo) == "master"
