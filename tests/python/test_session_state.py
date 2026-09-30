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
import json
import subprocess
from pathlib import Path
from types import ModuleType

import pytest

import harness_fixtures

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


# ------------------------------------------------ unattended mode in the recap
#
# h3 U1, T18 (docs/increments/h3-unattended-u1.md §3.8, §4). The copy of the
# script runs in a temporary repository (harness_fixtures), so its state dir is
# that repository's `.git/harness` and the real queue is never read.


@pytest.fixture
def harness_repo(tmp_path: Path) -> Path:
    return harness_fixtures.make_repo(tmp_path / "repo")


def _recap(repo: Path) -> list[str]:
    result = harness_fixtures.run_script(repo, "tools/session_state.py", "")
    assert result.returncode == 0, result.stderr
    return result.stdout.splitlines()


def _header(lines: list[str]) -> list[str]:
    """The lines printed before `== recap ==`."""
    return lines[: lines.index("== recap ==")]


def _queue_entry(n: int) -> dict[str, object]:
    return {
        "at": f"2026-09-30T{n:02d}:00:00+00:00",
        "branch": "feat",
        "cwd": "/r",
        "agent_type": None,
        "agent_id": None,
        "hook": "guard_push",
        "act": f"cmd-{n:02d}",
        "why": "w",
    }


def _write_queue(repo: Path, entries: list[str]) -> None:
    path = harness_fixtures.queue_path(repo)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(line + "\n" for line in entries))


def test_recap_while_unattended_opens_with_the_until(harness_repo: Path) -> None:
    until = harness_fixtures.flag_on(harness_repo)
    assert _header(_recap(harness_repo)) == [
        f"UNATTENDED until {until:%Y-%m-%d %H:%M} UTC. Guarded acts are refused and queued: "
        "record each as an ASK OLA line and continue."
    ]


def test_recap_with_a_broken_flag_opens_with_the_failure(harness_repo: Path) -> None:
    harness_fixtures.flag_broken(harness_repo, "not-json")
    [line] = _header(_recap(harness_repo))
    assert line.startswith("UNATTENDED FLAG UNREADABLE (")
    assert line.endswith("): every guarded act is refused. Ola: python3 tools/away.py --back.")


def test_recap_with_an_expired_flag_says_the_mode_ended(harness_repo: Path) -> None:
    harness_fixtures.flag_expired(harness_repo)
    [line] = _header(_recap(harness_repo))
    assert line.startswith("Unattended mode ended at ")
    assert line.endswith(". Ola: python3 tools/away.py --back prints and archives the queue.")


def test_recap_while_attended_has_no_header_and_no_harness_sections(harness_repo: Path) -> None:
    lines = _recap(harness_repo)
    assert _header(lines) == []
    assert not any(line.startswith("Queued while unattended") for line in lines)
    assert "Uncommitted rule-file changes:" not in lines


def test_recap_shows_the_newest_ten_queue_lines_after_waiting_on_ola(
    harness_repo: Path,
) -> None:
    _write_queue(harness_repo, [json.dumps(_queue_entry(n)) for n in range(12)])
    lines = _recap(harness_repo)
    waiting = next(i for i, line in enumerate(lines) if line.startswith("Waiting on Ola:"))
    at = lines.index("Queued while unattended (12):")
    assert at > waiting
    shown = lines[at + 1 : at + 11]
    assert [line.rsplit(" ", 1)[1] for line in shown] == [f"cmd-{n:02d}" for n in range(2, 12)]
    assert shown[-1] == "  2026-09-30T11:00:00+00:00 feat guard_push main: cmd-11"
    assert lines[at + 11] == "  ... 2 older"
    assert not any("cmd-00" in line or "cmd-01" in line for line in lines)


def test_recap_counts_unreadable_queue_lines(harness_repo: Path) -> None:
    _write_queue(harness_repo, [json.dumps(_queue_entry(1)), "{not json"])
    lines = _recap(harness_repo)
    assert "Queued while unattended (1):" in lines or "Queued while unattended (2):" in lines
    assert "  (1 unreadable lines)" in lines


def test_recap_lists_uncommitted_rule_files_in_every_worktree(
    harness_repo: Path, tmp_path: Path
) -> None:
    second = tmp_path / "second-tree"
    harness_fixtures.git(harness_repo, "worktree", "add", "-q", "-b", "side", str(second))
    (harness_repo / "CLAUDE.md").write_text("# changed here\n")
    (second / "CLAUDE.md").write_text("# changed there\n")
    (harness_repo / "notes.txt").write_text("not a rule\n")
    lines = _recap(harness_repo)
    at = lines.index("Uncommitted rule-file changes:")
    section: list[str] = []
    for line in lines[at + 1 :]:
        if not line.startswith("  "):
            break
        section.append(line)
    assert len(section) == 2
    assert any("/repo" in line and line.endswith("CLAUDE.md") for line in section)
    assert any("/second-tree" in line and line.endswith("CLAUDE.md") for line in section)
    assert not any("notes.txt" in line for line in section)
