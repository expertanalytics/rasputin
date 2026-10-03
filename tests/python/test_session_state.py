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
import os
import shlex
import shutil
import subprocess
import sys
import time
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
    # h8 §3.2 (Ola's ruling on item 3): matching is case-sensitive, so the
    # lowercase line no longer counts.
    assert session_state.pending_decisions(tmp_path) == [
        "session.md: ASK OLA: fix before landing?",
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


@pytest.mark.parametrize("branch", ["null", "missing"])
def test_recap_prints_a_null_or_missing_branch_as_no_branch_and_never_none(
    harness_repo: Path, branch: str
) -> None:
    # §3.8 (review 1): the same text `away.py --back` groups such entries under.
    entry = {key: value for key, value in _queue_entry(3).items() if key != "branch"}
    if branch == "null":
        entry["branch"] = None
    _write_queue(harness_repo, [json.dumps(entry)])
    lines = _recap(harness_repo)
    at = lines.index("Queued while unattended (1):")
    assert lines[at + 1] == "  2026-09-30T03:00:00+00:00 (no branch) guard_push main: cmd-03"
    assert not any("None" in line for line in lines)


def test_recap_has_no_uncommitted_rule_file_section(harness_repo: Path) -> None:
    # §3.8 (review 1): the scan moved to U3, so a modified rule file is not listed.
    (harness_repo / "CLAUDE.md").write_text("# changed here\n")
    lines = _recap(harness_repo)
    assert "Uncommitted rule-file changes:" not in lines
    assert not any(line.endswith(" CLAUDE.md") for line in lines)


# ------------------------------------------------ T21: a broken harness
#
# §3.8 (review 1): a harness fault costs the recap one header line, never the
# recap. The queue section is omitted; everything else prints; exit 0.

UNKNOWN_TAIL = (
    "): the recap could not read the flag or the queue. Guarded acts may be refused; "
    "record each refusal as an ASK OLA line and continue."
)


def _break_harness(repo: Path, how: str) -> None:
    module = repo / "tools" / "harness_mode.py"
    if how == "missing":
        module.unlink()
    else:
        with module.open("a") as out:
            out.write('\n\ndef read_mode(*a, **k):\n    raise RuntimeError("planted")\n')


@pytest.mark.parametrize("how", ["missing", "raising"])
def test_recap_with_a_broken_harness_says_so_and_prints_the_rest(
    harness_repo: Path, how: str
) -> None:
    harness_fixtures.flag_on(harness_repo)
    _write_queue(harness_repo, [json.dumps(_queue_entry(1))])
    tasks = harness_repo / ".claude" / "current-task"
    tasks.mkdir(parents=True)
    (tasks / "session.md").write_text("ASK OLA: may I push?\n")
    _break_harness(harness_repo, how)

    lines = _recap(harness_repo)

    if how == "missing":
        assert lines[0].startswith("UNATTENDED STATE UNKNOWN (ModuleNotFoundError: ")
        assert lines[0].endswith(UNKNOWN_TAIL)
    else:
        assert lines[0] == f"UNATTENDED STATE UNKNOWN (RuntimeError: planted{UNKNOWN_TAIL}"
    assert lines[1] == "== recap =="
    waiting = next(i for i, line in enumerate(lines) if line.startswith("Waiting on Ola:"))
    assert lines[waiting + 1] == "  session.md: ASK OLA: may I push?"
    assert "Next on ROADMAP.md:" in lines
    assert not any(line.startswith("Queued while unattended") for line in lines)


# ================================================================ h8
#
# docs/increments/h8-window-and-recap.md §3.1-§3.7 and test plan §6. Fixture
# repositories are made under `tmp_path.resolve()`, so the paths git prints
# (realpaths, /private/var on macOS) and the ones the tests build agree.

ASK_EMPTY = "WARNING: empty ASK OLA line"


def _tasks(checkout: Path) -> Path:
    tasks = checkout / ".claude" / "current-task"
    tasks.mkdir(parents=True, exist_ok=True)
    return tasks


def _porcelain_worktrees(repo: Path) -> list[Path]:
    """The `worktree <path>` lines of `git worktree list --porcelain`, as printed."""
    text = harness_fixtures.git(repo, "worktree", "list", "--porcelain")
    return [
        Path(line[len("worktree ") :]) for line in text.splitlines() if line.startswith("worktree ")
    ]


# ---------------------------------------------------------------- test 1, §3.2

COUNTED = ("ASK OLA: a", "- ASK OLA: b", "* ASK OLA: c", "  + ASK OLA: d", "\t-\tASK OLA: e")
NOT_COUNTED = ("ask ola: x", "Note: ASK OLA: x", "**ASK OLA:** x", "1. ASK OLA: x", "ASK OLA x")


def test_ask_ola_counts_only_at_the_start_after_an_optional_bullet(tmp_path: Path) -> None:
    (tmp_path / "session.md").write_text("\n".join((*COUNTED, *NOT_COUNTED)) + "\n")
    assert session_state.pending_decisions(tmp_path) == [
        "session.md: ASK OLA: a",
        "session.md: - ASK OLA: b",
        "session.md: * ASK OLA: c",
        "session.md: + ASK OLA: d",
        "session.md: -\tASK OLA: e",
    ]


@pytest.mark.parametrize("line", NOT_COUNTED)
def test_a_line_that_only_mentions_ask_ola_is_not_counted(tmp_path: Path, line: str) -> None:
    (tmp_path / "tester-101010.md").write_text(line + "\n")
    assert session_state.pending_decisions(tmp_path) == []


def test_an_empty_ask_ola_line_is_listed_as_a_warning_with_its_line_number(
    tmp_path: Path,
) -> None:
    (tmp_path / "session.md").write_text("NOW: x\nASK OLA:\nASK OLA: real\n- ASK OLA:   \n")
    assert session_state.pending_decisions(tmp_path) == [
        f"session.md:2: {ASK_EMPTY}",
        "session.md: ASK OLA: real",
        f"session.md:4: {ASK_EMPTY}",
    ]


def test_a_counted_line_is_listed_in_full_however_long(tmp_path: Path) -> None:
    line = "ASK OLA: " + "q" * 400
    (tmp_path / "session.md").write_text(line + "\n")
    assert session_state.pending_decisions(tmp_path) == [f"session.md: {line}"]


# ---------------------------------------------------------------- test 2, §3.1


@pytest.fixture
def forest(tmp_path: Path) -> tuple[Path, Path, Path]:
    """A main checkout, a worktree inside it and one outside it."""
    root = tmp_path.resolve()
    main = harness_fixtures.make_repo(root / "main")
    inside = harness_fixtures.add_worktree(main, main / ".claude" / "worktrees" / "w1", "w1")
    outside = harness_fixtures.add_worktree(main, root / "elsewhere" / "w2", "w2")
    return main, inside, outside


def test_main_checkout_is_the_main_path_from_anywhere(forest: tuple[Path, Path, Path]) -> None:
    main, inside, outside = forest
    for start in (main, inside, outside):
        assert session_state.main_checkout(start).resolve() == main.resolve()


def test_checkouts_lists_the_main_checkout_first_then_every_worktree(
    forest: tuple[Path, Path, Path],
) -> None:
    main, inside, outside = forest
    for start in (main, inside, outside):
        found = [p.resolve() for p in session_state.checkouts(start)]
        assert found[0] == main.resolve()
        assert sorted(found[1:]) == sorted([inside.resolve(), outside.resolve()])


def test_checkouts_skips_a_worktree_whose_directory_was_deleted(
    forest: tuple[Path, Path, Path],
) -> None:
    main, inside, outside = forest
    shutil.rmtree(outside)
    assert [p.resolve() for p in session_state.checkouts(main)] == [
        main.resolve(),
        inside.resolve(),
    ]


def test_outside_a_repository_both_fall_back_to_the_given_path(tmp_path: Path) -> None:
    plain = tmp_path.resolve() / "plain"
    plain.mkdir()
    assert session_state.main_checkout(plain) == plain
    assert session_state.checkouts(plain) == [plain]


def test_a_bare_repository_is_its_own_main_checkout(tmp_path: Path) -> None:
    bare = tmp_path.resolve() / "bare.git"
    subprocess.run(
        ["git", "init", "-q", "--bare", str(bare)], check=True, env=harness_fixtures.clean_env()
    )
    assert session_state.main_checkout(bare) == bare


# ---------------------------------------------------------------- test 3, §3.2


def _decide(forest: tuple[Path, Path, Path]) -> list[str]:
    """Decisions in all three checkouts; the expected list in checkout order."""
    main, inside, outside = forest
    (_tasks(main) / "session.md").write_text("NOW: n\nQUEUE: q\nASK OLA: main?\n")
    (_tasks(inside) / "tester-101010.md").write_text("ASK OLA: inside?\n")
    (_tasks(outside) / "developer-121212.md").write_text("- ASK OLA: outside?\n")
    by_path = {
        main.resolve(): ["session.md: ASK OLA: main?"],
        inside.resolve(): [".claude/worktrees/w1: tester-101010.md: ASK OLA: inside?"],
    }
    for listed in _porcelain_worktrees(main):
        if listed.resolve() == outside.resolve():
            by_path[outside.resolve()] = [f"{listed}: developer-121212.md: - ASK OLA: outside?"]
    order = [p.resolve() for p in _porcelain_worktrees(main)]
    return [line for p in order for line in by_path[p]]


def test_all_decisions_lists_every_checkout_with_its_prefix(
    forest: tuple[Path, Path, Path],
) -> None:
    expected = _decide(forest)
    assert len(expected) == 3
    for start in forest:
        assert session_state.all_decisions(start) == expected


def test_the_recap_run_from_a_worktree_reads_the_main_checkouts_task_files(
    forest: tuple[Path, Path, Path],
) -> None:
    main, inside, _ = forest
    expected = _decide(forest)
    (_tasks(inside) / "session.md").write_text("NOW: the worktree's own ask\n")
    lines = _recap(inside)
    waiting = lines.index("Waiting on Ola:")
    assert [line.strip() for line in lines[waiting + 1 : waiting + 4]] == expected
    at = lines.index("== .claude/current-task/session.md ==")
    assert lines[at + 1 : at + 4] == ["NOW: n", "QUEUE: q", "ASK OLA: main?"]
    assert "NOW: the worktree's own ask" not in lines
    assert (main / ".claude" / "current-task" / "session.md").exists()


# ---------------------------------------------------------------- test 4, §3.3

NOT_A = "not a NOW, QUEUE or ASK OLA line: "
GOOD = "NOW: build h8\nQUEUE: red, green\nASK OLA: push?\n\nASK OLA: merge?\n"
GOOD_BULLETED = "- NOW: build h8\n* QUEUE: red\n  + ASK OLA: push?\n\t- ASK OLA: merge?\n"


@pytest.mark.parametrize("text", [GOOD, GOOD_BULLETED])
def test_a_well_formed_session_md_gives_no_warning(text: str) -> None:
    assert session_state.session_format(text) == []


@pytest.mark.parametrize(
    ("text", "warning"),
    [
        ("QUEUE: q\n", "session.md: 0 NOW lines; exactly one expected"),
        ("NOW: n\n", "session.md: 0 QUEUE lines; exactly one expected"),
        ("NOW: n\nQUEUE: a\nQUEUE: b\n", "session.md: 2 QUEUE lines; exactly one expected"),
        ("NOW: a\nNOW: b\nQUEUE: q\n", "session.md: 2 NOW lines; exactly one expected"),
        ("NOW:   \nQUEUE: q\n", "session.md:1: empty NOW line"),
        ("NOW: n\n- QUEUE:\n", "session.md:2: empty QUEUE line"),
        ("NOW: n\nQUEUE: q\nASK OLA: x\nRulings: Ola said y\n",
         f"session.md:4: {NOT_A}Rulings: Ola said y"),
    ],
)  # fmt: skip
def test_each_fault_gives_its_warning_once(text: str, warning: str) -> None:
    assert session_state.session_format(text) == [warning]


def test_a_lowercase_kind_is_not_a_line_of_that_kind() -> None:
    assert session_state.session_format("now: x\nQUEUE: q\n") == [
        "session.md: 0 NOW lines; exactly one expected",
        f"session.md:1: {NOT_A}now: x",
    ]


def test_the_not_a_line_warning_quotes_the_first_60_characters() -> None:
    line = "Rulings: " + "r" * 100
    assert session_state.session_format(f"NOW: n\nQUEUE: q\n{line}\n") == [
        f"session.md:3: {NOT_A}{line[:60]}"
    ]


def test_an_empty_ask_ola_line_is_not_warned_twice() -> None:
    assert session_state.session_format("NOW: n\nQUEUE: q\nASK OLA:\n") == []


def test_the_warnings_come_in_the_order_of_the_design() -> None:
    long_ask = "ASK OLA: " + "y" * 300  # 309 characters
    text = f"Rulings: z\nNOW:\nNOW: a\n{long_ask}\nQUEUE: q\nQUEUE: r\n"
    assert session_state.session_format(text) == [
        "session.md: 2 NOW lines; exactly one expected",
        "session.md: 2 QUEUE lines; exactly one expected",
        "session.md:2: empty NOW line",
        f"session.md:1: {NOT_A}Rulings: z",
        "session.md:4: 309 characters; at most 300",
    ]


def test_the_name_is_the_one_given() -> None:
    assert session_state.session_format("QUEUE: q\n", name="x.md") == [
        "x.md: 0 NOW lines; exactly one expected"
    ]


def test_the_line_limit_is_300() -> None:
    assert session_state.MAX_LINE == 300


def test_a_line_of_exactly_300_characters_passes() -> None:
    line = "NOW: " + "a" * 295
    assert len(line) == 300
    assert session_state.session_format(f"{line}\nQUEUE: q\n") == []


def test_a_line_of_301_characters_warns_with_its_number_and_length() -> None:
    line = "NOW: " + "a" * 296
    assert session_state.session_format(f"QUEUE: q\n{line}\n") == [
        "session.md:2: 301 characters; at most 300"
    ]


def test_trailing_whitespace_does_not_count_toward_the_limit() -> None:
    line = "NOW: " + "a" * 295 + "   \t "
    assert session_state.session_format(f"{line}\nQUEUE: q\n") == []


def test_a_long_ask_ola_line_warns_too() -> None:
    line = "- ASK OLA: " + "z" * 300  # 311 characters
    assert session_state.session_format(f"NOW: n\nQUEUE: q\n{line}\n") == [
        "session.md:3: 311 characters; at most 300"
    ]


def _block(lines: list[str], heading: str) -> list[str]:
    """The indented lines under `heading`, stripped."""
    at = lines.index(heading)
    block = []
    for line in lines[at + 1 :]:
        if not line.startswith("  "):
            break
        block.append(line.strip())
    return block


def test_the_recap_prints_format_warnings_for_a_faulty_main_session_md(
    harness_repo: Path,
) -> None:
    (_tasks(harness_repo) / "session.md").write_text("ASK OLA: push?\nRulings: x\n")
    lines = _recap(harness_repo)
    waiting = lines.index("Waiting on Ola:")
    assert lines[waiting + 1] == "  session.md: ASK OLA: push?"
    assert lines[waiting + 2] == "session.md format:"
    assert _block(lines, "session.md format:") == [
        "session.md: 0 NOW lines; exactly one expected",
        "session.md: 0 QUEUE lines; exactly one expected",
        f"session.md:2: {NOT_A}Rulings: x",
    ]


def test_the_recap_prints_no_format_block_for_a_good_or_absent_session_md(
    harness_repo: Path,
) -> None:
    assert "session.md format:" not in _recap(harness_repo)
    (_tasks(harness_repo) / "session.md").write_text(GOOD)
    assert "session.md format:" not in _recap(harness_repo)


def test_the_recap_does_not_format_check_a_worktrees_session_md(
    forest: tuple[Path, Path, Path],
) -> None:
    main, inside, _ = forest
    (_tasks(main) / "session.md").write_text(GOOD)
    (_tasks(inside) / "session.md").write_text("Rulings: x\nASK OLA: from w1?\n")
    for start in (main, inside):
        lines = _recap(start)
        assert "session.md format:" not in lines
        assert "  .claude/worktrees/w1: session.md: ASK OLA: from w1?" in lines


def test_the_recap_caps_the_format_block_at_four_lines(harness_repo: Path) -> None:
    faults = "".join(f"Ruling {n}\n" for n in range(6))
    (_tasks(harness_repo) / "session.md").write_text(f"NOW: n\nQUEUE: q\n{faults}")
    assert _block(_recap(harness_repo), "session.md format:") == [
        *(f"session.md:{n + 3}: {NOT_A}Ruling {n}" for n in range(4)),
        "... 2 more",
    ]


# ---------------------------------------------------------------- test 5, §3.5
#
# The wrapper line as observed (§3.5): a shell that sources a snapshot under
# ~/.claude/shell-snapshots/ and evals the command.

JOBS = "Running Bash-tool jobs, any session (Claude Code shell wrappers):"
NOT_RECOGNISED = "(wrapper format not recognised)"


def wrapper(command: str) -> str:
    return (
        "/bin/zsh -c source /Users/o/.claude/shell-snapshots/snapshot-zsh-1-ab.sh 2>/dev/null"
        " || true && setopt NO_EXTENDED_GLOB 2>/dev/null || true && "
        f"eval '{command}' < /dev/null && pwd -P >| /tmp/claude-1a2b-cwd"
    )


def ps_line(pid: int, ppid: int, etime: str, command: str) -> str:
    """One line of `ps -e -ww -o pid=,ppid=,etime=,command=`."""
    return f"{pid:5} {ppid:5} {etime:>11} {command}"


PS_TEXT = "\n".join(
    [
        ps_line(1, 0, "10-01:02:03", "/sbin/launchd"),
        ps_line(48413, 58828, "09:41", wrapper("pytest   tests/python")),
        ps_line(48420, 48413, "09:40", "/usr/bin/python3 -m pytest tests/python"),
        ps_line(52736, 58828, "1-02:03:04", wrapper("sleep 600")),
        ps_line(52800, 58828, "00:01", "/bin/zsh -c eval 'no snapshot here' < /dev/null"),
        ps_line(52801, 58828, "00:01", "/bin/zsh -c source /x/.claude/shell-snapshots/s.sh"),
        ps_line(58828, 1, "2-00:00:00", "/opt/homebrew/bin/claude --resume"),
    ]
)


def test_background_jobs_lists_wrapper_lines_with_pid_etime_and_command() -> None:
    assert session_state.background_jobs(PS_TEXT, set()) == [
        "pid 48413, running 09:41: pytest tests/python",
        "pid 52736, running 1-02:03:04: sleep 600",
    ]


def test_background_jobs_skips_excluded_pids() -> None:
    assert session_state.background_jobs(PS_TEXT, {48413}) == [
        "pid 52736, running 1-02:03:04: sleep 600"
    ]


def test_background_jobs_cuts_a_long_command_to_100_characters() -> None:
    command = "echo " + "x" * 150
    text = ps_line(700, 1, "00:05", wrapper(command))
    assert session_state.background_jobs(text, set()) == [
        f"pid 700, running 00:05: {command[:100]}"
    ]


def test_background_jobs_takes_the_command_to_the_end_of_an_unterminated_line() -> None:
    text = ps_line(
        701, 1, "00:05", "/bin/zsh -c source /u/.claude/shell-snapshots/s.sh && eval 'sleep 5 && ls"
    )
    assert session_state.background_jobs(text, set()) == ["pid 701, running 00:05: sleep 5 && ls"]


def test_background_jobs_shows_four_then_how_many_more() -> None:
    text = "\n".join(ps_line(800 + n, 1, "00:05", wrapper(f"job {n}")) for n in range(6))
    assert session_state.background_jobs(text, set()) == [
        *(f"pid {800 + n}, running 00:05: job {n}" for n in range(4)),
        "... 2 more",
    ]


def test_background_jobs_of_an_empty_listing_is_empty() -> None:
    assert session_state.background_jobs("", set()) == []


HOOK_LAUNCHER = '/bin/sh -c python3 "$CLAUDE_PROJECT_DIR/tools/session_state.py"'
SELF = "/usr/bin/python3 /r/tools/session_state.py"


@pytest.mark.parametrize(
    ("ancestors", "verdict"),
    [
        ([SELF, wrapper("python3 x.py"), "/opt/homebrew/bin/claude"], "confirmed"),
        ([wrapper("x"), "claude --resume"], "confirmed"),
        # §8.7: a recognised wrapper confirms the format even with no Claude Code above it.
        ([SELF, wrapper("x"), "/sbin/launchd"], "confirmed"),
        ([SELF, "/bin/zsh -c run it some new way", "/usr/local/bin/claude"], "drift"),
        ([SELF, "bash -c x", "claude.exe"], "drift"),
        ([SELF, HOOK_LAUNCHER, "/opt/homebrew/bin/claude"], "unknown"),
        ([SELF, "/bin/zsh -l", "/Applications/iTerm.app/Contents/MacOS/iTerm2"], "unknown"),
        ([SELF, "/bin/zsh -c make", "/sbin/launchd"], "unknown"),
        # It stops at the first Claude Code process: a shell above it is not looked at.
        ([SELF, "/opt/homebrew/bin/claude", "/bin/zsh -c other"], "unknown"),
        ([], "unknown"),
    ],
)  # fmt: skip
def test_wrapper_check(ancestors: list[str], verdict: str) -> None:
    assert session_state.wrapper_check(ancestors) == verdict


def _fake_ps(tmp_path: Path, script: str) -> Path:
    """A directory holding an executable `ps` that runs `script` under /bin/sh."""
    bin_dir = tmp_path / "fake-bin"
    bin_dir.mkdir(exist_ok=True)
    ps = bin_dir / "ps"
    ps.write_text("#!/bin/sh\n" + script)
    ps.chmod(0o755)
    return bin_dir


def _recap_with(repo: Path, bin_dir: Path) -> list[str]:
    """The recap with `bin_dir` first on PATH, so `ps` is the fake."""
    env = harness_fixtures.clean_env()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    result = subprocess.run(
        [sys.executable, str(repo / "tools" / "session_state.py")],
        input="",
        capture_output=True,
        text=True,
        cwd=repo,
        env=env,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout.splitlines()


def test_the_recap_lists_a_running_wrapper_by_its_pid(harness_repo: Path, tmp_path: Path) -> None:
    # Real processes through the real `ps`, filtered to this test's own wrapper
    # so other sessions' jobs on the machine cannot crowd it past the cap.
    real_ps = shutil.which("ps")
    assert real_ps is not None
    marker = f"h8-test-{os.getpid()}-{time.monotonic_ns()}"
    bin_dir = _fake_ps(
        tmp_path,
        f'"{real_ps}" "$@" | awk -v m="{marker}" '
        "'index($0, \"/.claude/shell-snapshots/\") == 0 || index($0, m) > 0'\n",
    )
    job = subprocess.Popen(
        ["/bin/sh", "-c", f": /x/.claude/shell-snapshots/{marker}.sh"
         " && eval 'sleep 60' < /dev/null && pwd -P >| /dev/null"],
        start_new_session=True,
    )  # fmt: skip
    try:
        lines = _recap_with(harness_repo, bin_dir)
    finally:
        job.kill()
        job.wait()
    block = _block(lines, JOBS)
    assert len(block) == 1, block
    assert block[0].startswith(f"pid {job.pid}, running ")
    assert block[0].endswith(": sleep 60")


def test_the_recap_places_the_jobs_after_in_flight_and_before_waiting_on_ola(
    harness_repo: Path, tmp_path: Path
) -> None:
    lines = _recap_with(harness_repo, _fake_ps(tmp_path, "exit 0\n"))
    assert (
        lines.index("In flight:")
        < lines.index(JOBS)
        < lines.index("Waiting on Ola: (none recorded as ASK OLA)")
    )
    assert _block(lines, JOBS) == ["(none)"]


def test_the_recap_says_when_ps_fails_and_prints_the_rest(
    harness_repo: Path, tmp_path: Path
) -> None:
    lines = _recap_with(harness_repo, _fake_ps(tmp_path, "echo 'ps: planted' >&2\nexit 1\n"))
    [shown] = _block(lines, JOBS)
    assert shown.startswith("(could not list: ")
    assert "Next on ROADMAP.md:" in lines
    assert "== .claude/current-task/session.md ==" in lines


SELF_LINE = "/usr/bin/python3 tools/session_state.py"


def _printf(pid: int | str, ppid: int, etime: str, command: str) -> str:
    """A shell line printing one `ps` line; pid `$PPID` is the recap's own pid."""
    return f'printf "%5d %5d %11s %s\\n" {pid} {ppid} {etime} {shlex.quote(command)}\n'


def test_the_recap_says_the_wrapper_format_drifted_instead_of_none(
    harness_repo: Path, tmp_path: Path
) -> None:
    # The recap's parent is a `zsh -c` of another shape under `claude`: drift.
    script = (
        _printf('"$PPID"', 900002, "00:01", SELF_LINE)
        + _printf(900002, 900003, "00:05", "/bin/zsh -c run it some new way")
        + _printf(900003, 1, "01:00", "/opt/homebrew/bin/claude")
    )
    lines = _recap_with(harness_repo, _fake_ps(tmp_path, script))
    assert _block(lines, JOBS) == [NOT_RECOGNISED]


def test_the_recap_under_the_hook_launcher_says_none(harness_repo: Path, tmp_path: Path) -> None:
    script = (
        _printf('"$PPID"', 900002, "00:01", SELF_LINE)
        + _printf(900002, 900003, "00:01", HOOK_LAUNCHER)
        + _printf(900003, 1, "01:00", "/opt/homebrew/bin/claude")
    )
    lines = _recap_with(harness_repo, _fake_ps(tmp_path, script))
    assert _block(lines, JOBS) == ["(none)"]


def test_the_recap_does_not_list_its_own_wrapper(harness_repo: Path, tmp_path: Path) -> None:
    # The recap's parent is a recognised wrapper under `claude`: confirmed, and
    # the parent is an ancestor, so it is excluded; another session's job is listed.
    script = (
        _printf('"$PPID"', 900002, "00:01", SELF_LINE)
        + _printf(900002, 900003, "00:01", wrapper("python3 tools/session_state.py"))
        + _printf(900004, 900003, "00:09", wrapper("sleep 600"))
        + _printf(900003, 1, "01:00", "/opt/homebrew/bin/claude")
    )
    lines = _recap_with(harness_repo, _fake_ps(tmp_path, script))
    assert _block(lines, JOBS) == ["pid 900004, running 00:09: sleep 600"]


# ---------------------------------------------------------------- test 9, §3.6


def _sized_repo(repo: Path) -> None:
    (repo / ".claude" / "agents").mkdir(parents=True)
    (repo / ".claude" / "agents" / "tester.md").write_text("test " * 8)
    retro = repo / "docs" / "retrospectives" / "2026-10-01-first.md"
    retro.parent.mkdir(parents=True)
    retro.write_text("# Retro\n")
    harness_fixtures.git(repo, "add", "-A")
    harness_fixtures.git(repo, "commit", "-q", "-m", "rule files and a retrospective")
    (repo / "CLAUDE.md").write_text("# rules\n" + "more " * 12)
    harness_fixtures.git(repo, "commit", "-q", "-am", "grow CLAUDE.md")


def test_the_recap_prints_the_size_table_last(harness_repo: Path) -> None:
    assert (harness_repo / "tools" / "rule_sizes.py").exists(), "rule_sizes.py is not copied"
    _sized_repo(harness_repo)
    lines = _recap(harness_repo)
    heading = next(i for i, line in enumerate(lines) if line.startswith("Rule text in words"))
    assert heading > lines.index("Next on ROADMAP.md:")
    assert heading < lines.index("== .claude/current-task/session.md ==")
    assert lines[heading].startswith("Rule text in words, change since 2026-10-01-first.md (")
    assert f"  {14:5} {'+12':>6} CLAUDE.md" in lines
    assert f"  {8:5} {'0':>6} .claude/agents/tester.md" in lines


def test_the_recap_with_a_broken_size_module_prints_one_line_and_the_rest(
    harness_repo: Path,
) -> None:
    (harness_repo / "tools" / "rule_sizes.py").write_text('raise RuntimeError("planted")\n')
    lines = _recap(harness_repo)
    unavailable = [line for line in lines if "(size table unavailable: " in line]
    assert len(unavailable) == 1
    assert "planted" in unavailable[0]
    assert not any(line.startswith("Rule text in words") for line in lines)
    assert lines.index("Next on ROADMAP.md:") < lines.index(unavailable[0])
    assert "== .claude/current-task/session.md ==" in lines


# ---------------------------------------------------------------- test 11, §3.6
#
# The four sections h8 adds or widens, at their caps, fit 3,000 characters
# together. Each block runs from its heading to the line before the next
# section's heading, newlines included, so the fixture's transcripts and
# ROADMAP.md do not enter the sum.

BUDGET = 3000
RULE_GLOBS = (".claude/agents/*.md", ".claude/skills/**/*.md")


def _copy_real_rule_files(repo: Path) -> int:
    """The real checkout's rule files at their paths; how many."""
    real = harness_fixtures.REAL
    paths = [real / p for p in ("CLAUDE.md", ".claude/REQUIRED-READING.md",
                                "docs/increments/README.md", "docs/PRINCIPLES.md")]  # fmt: skip
    for pattern in RULE_GLOBS:
        paths += sorted(real.glob(pattern))
    for source in paths:
        target = repo / source.relative_to(real)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return len(paths)


def _blocks(lines: list[str]) -> dict[str, str]:
    """The text of each budgeted section, heading to the next heading."""
    starts = {
        "jobs": lines.index(JOBS),
        "waiting": next(i for i, line in enumerate(lines) if line.startswith("Waiting on Ola:")),
        "format": lines.index("session.md format:"),
        "sizes": next(i for i, line in enumerate(lines) if line.startswith("Rule text in words")),
    }
    ends = {
        "jobs": starts["waiting"],
        "waiting": starts["format"],
        "format": next(
            i for i, line in enumerate(lines)
            if i > starts["format"]
            and (line.startswith("Queued while") or line == "Next on ROADMAP.md:")
        ),
        "sizes": lines.index("== .claude/current-task/session.md =="),
    }  # fmt: skip
    return {k: "".join(line + "\n" for line in lines[starts[k] : ends[k]]) for k in starts}


def test_the_four_sections_fit_the_budget_at_their_caps(tmp_path: Path) -> None:
    root = tmp_path.resolve()
    main = harness_fixtures.make_repo(root / "main")
    files = _copy_real_rule_files(main)
    retro = main / "docs" / "retrospectives" / "2026-10-01-first.md"
    retro.parent.mkdir(parents=True, exist_ok=True)
    retro.write_text("# Retro\n")
    harness_fixtures.git(main, "add", "-A")
    harness_fixtures.git(main, "commit", "-q", "-m", "real rule files and a retrospective")
    # 30 ASK OLA lines of 300 characters, ten in each of three worktrees.
    for n in range(3):
        tree = harness_fixtures.add_worktree(
            main, main / ".claude" / "worktrees" / f"w{n}", f"w{n}"
        )
        asks = "".join(f"ASK OLA: {n}-{k} " + "a" * (300 - 14) + "\n" for k in range(10))
        (_tasks(tree) / "tester-101010.md").write_text(asks)
    # 15 format faults in the main session.md, each a long "not a" line.
    faults = "".join(f"Ruling {k:02d}: " + "r" * 120 + "\n" for k in range(15))
    (_tasks(main) / "session.md").write_text(f"NOW: n\nQUEUE: q\n{faults}")
    # 15 wrapper processes whose eval'd commands are 300 characters long.
    command = "sleep 600; : " + "c" * (300 - 13)
    script = "".join(_printf(100000 + k, 1, "10-02:03:04", wrapper(command)) for k in range(15))
    lines = _recap_with(main, _fake_ps(tmp_path, script))

    blocks = _blocks(lines)

    assert _block(lines, JOBS)[-1] == "... 11 more"
    assert len(_block(lines, JOBS)) == 5
    waiting = _block(lines, next(line for line in lines if line.startswith("Waiting on Ola:")))
    assert len(waiting) == 6 and waiting[-1] == "... 25 more"
    assert all(line.endswith("...") and len(line) <= 163 for line in waiting[:5])
    assert _block(lines, "session.md format:")[-1] == "... 11 more"
    assert len(_block(lines, "session.md format:")) == 5
    assert len(blocks["sizes"].splitlines()) >= files + 2
    total = sum(len(text) for text in blocks.values())
    assert total <= BUDGET, {k: len(v) for k, v in blocks.items()}


def test_a_long_ask_ola_line_is_shown_cut_to_160_characters(harness_repo: Path) -> None:
    line = "ASK OLA: " + "b" * 291  # 300 characters
    (_tasks(harness_repo) / "session.md").write_text(f"NOW: n\nQUEUE: q\n{line}\n")
    [shown] = _block(_recap(harness_repo), "Waiting on Ola:")
    # §8.11: the first 157 characters, then `...`, 160 in all.
    assert len(shown) == 160
    assert shown == f"session.md: {line}"[:157] + "..."
    # The full line is still printed further down, with session.md.
    assert line in _recap(harness_repo)
