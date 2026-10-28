"""`tools/work_pipes.py`: named work lists, the order a free writer takes, and resurfacing.

The spec is `docs/increments/h31-pipes.md`, *11. Tests (`tests/python/test_work_pipes.py`)*.
The pure functions take `now`, `unattended` and the busy branches as
arguments, so no test sleeps or patches the clock; the fixture time is
2026-10-25T19:00:00Z and every age is measured from it. Commands run as the
harness runs them, as a copied script in a temporary repository
(`harness_fixtures.make_repo`), so the real live lists are never read.

Records are built as the dictionaries the live list holds, one JSON object per
line, and folded by `pipes.fold`; a pipe's priority comes from a `Pipe`.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import subprocess
from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import harness_fixtures
from harness_fixtures import REAL, clean_env, git, make_repo, run_script

# A harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

pipes = harness_fixtures.Tool("work_pipes")

NOW = datetime(2026, 10, 25, 19, 0, tzinfo=UTC)
HOUR, DAY = timedelta(hours=1), timedelta(days=1)
IDLE = "IDLE: no runnable item in any pipe"

CORE_TOML = 'name = "core"\npurpose = "the public repository: product, harness and process"\n'
LSM_TOML = 'name = "lsm"\npurpose = "the private repository\'s land-surface work"\n'


# ---------------------------------------------------------------- helpers


def rec(item: str, state: str, pipe: str = "core", age: timedelta = HOUR, **extra: Any
        ) -> dict[str, Any]:  # fmt: skip
    """One live-list record for `item`, its state changed `age` before NOW."""
    record: dict[str, Any] = {
        "id": item, "title": f"title of {item}", "state": state, "pipe": pipe,
        "next": "", "branch": f"worktree-{item}", "pr": None, "question": None,
        "default": None, "plain": False, "fallback": False,
        "since": (NOW - age).isoformat(timespec="seconds"),
    }  # fmt: skip
    record.update(extra)
    return record


def two_pipes() -> list[Any]:
    return [
        pipes.Pipe(name="core", purpose="p", priority=1, source="tree"),
        pipes.Pipe(name="lsm", purpose="p", priority=2, source="tree"),
    ]


def choose(records: list[dict[str, Any]], busy: set[str] | None = None,
           unattended: bool = False, pipe_list: list[Any] | None = None) -> Any:  # fmt: skip
    return pipes.order(
        pipe_list if pipe_list is not None else two_pipes(),
        pipes.fold(records), NOW, unattended, busy or set(),
    )  # fmt: skip


def label(choice: Any) -> str:
    """`<pipe>/<id>` of the chosen item, or `IDLE` when none is chosen."""
    if choice.item is None:
        return "IDLE"
    return f"{choice.item.pipe}/{choice.item.id}"


def common_dir(repo: Path) -> Path:
    found = Path(git(repo, "rev-parse", "--git-common-dir").strip())
    return found if found.is_absolute() else (repo / found).resolve()


def live_list(repo: Path, pipe: str) -> Path:
    return common_dir(repo) / "pipes" / f"{pipe}.jsonl"


def append(repo: Path, *records: dict[str, Any]) -> None:
    """Append records to their pipes' live lists, as `tools/work_pipes.py` would."""
    for record in records:
        path = live_list(repo, record["pipe"])
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a") as out:
            out.write(json.dumps(record) + "\n")


def note(repo: Path, name: str, persona: str = "developer") -> Path:
    """A brief's note file for worktree `name`, as `tools/brief.py` names it."""
    path = repo / ".claude" / "current-task" / f"{persona}-{name}-190000.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ask\npersona\nfile\n")
    return path


def cli(repo: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return run_script(repo, "tools/work_pipes.py", "", *args)


def ok(repo: Path, *args: str) -> str:
    result = cli(repo, *args)
    assert result.returncode == 0, f"work_pipes.py {' '.join(args)}: {result.stderr}"
    return result.stdout


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """The fixture repository with `core` and `lsm` defined, `core` first by priority."""
    root = make_repo(tmp_path / "repo")
    assert (root / "tools" / "work_pipes.py").exists(), "work_pipes.py missing from the copy"
    (root / ".claude" / "pipes").mkdir(parents=True)
    (root / ".claude" / "pipes" / "core.toml").write_text(CORE_TOML)
    (root / ".claude" / "pipes" / "lsm.toml").write_text(LSM_TOML)
    ok(root, "priority", "core=1", "lsm=2", "--words", "core first, then the private work")
    return root


def after(lines: list[str], header: str) -> list[str]:
    """The lines after the first one that starts with `header`."""
    start = next(i for i, line in enumerate(lines) if line.startswith(header))
    return lines[start + 1 :]


def first_index(lines: list[str], prefix: str) -> int:
    found = [i for i, line in enumerate(lines) if line.startswith(prefix)]
    assert found, f"no line starts with {prefix!r}: {lines}"
    return found[0]


COLUMNS = ("PIPE", "PRI", "ITEM", "STEP", "WAITING ON", "SINCE", "NEXT")


def table(lines: list[str]) -> list[dict[str, str]]:
    """The table's rows, each cell cut at the column heading's offset, `mark` first.

    A row is any line after the column heading row and before `Questions for
    Ola:` or the closing `next:` line. Each column must start where its heading
    starts, so a cell that spills into the next column fails here.
    """
    at = [i for i, line in enumerate(lines) if all(name in line for name in COLUMNS)]
    assert len(at) == 1, f"one column heading row expected: {lines}"
    heading = lines[at[0]]
    offsets = [heading.index(name) for name in COLUMNS]
    assert offsets == sorted(offsets), heading
    rows = []
    for line in lines[at[0] + 1 :]:
        if line.startswith(("Questions for Ola:", "next:")):
            break
        for name, offset in zip(COLUMNS, offsets, strict=True):
            if name == "NEXT" and len(line) <= offset:
                continue  # an empty NEXT cell, its trailing blanks stripped
            assert line[offset - 1] == " ", f"{name} does not start at {offset}: {line!r}"
            if name != "NEXT":
                assert line[offset] != " ", f"{name} cell empty or shifted: {line!r}"
        bounds = [0, *offsets, len(line)]
        cells = [line[a:b].strip() for a, b in itertools.pairwise(bounds)]
        rows.append(dict(zip(("mark", *COLUMNS), cells, strict=True)))
    return rows


# ---------------------------------------------------------------- 1. fold


def test_fold_keeps_the_newest_record_per_item() -> None:
    items = pipes.fold([rec("h31", "queued", age=2 * DAY), rec("h31", "designing", age=DAY),
                        rec("h31", "red", age=HOUR)])  # fmt: skip
    assert list(items) == ["h31"]
    assert items["h31"].state == "red"
    assert items["h31"].since == NOW - HOUR


def test_an_unreadable_line_is_counted_not_raised(repo: Path) -> None:
    path = live_list(repo, "core")
    path.parent.mkdir(parents=True, exist_ok=True)
    good = [json.dumps(rec("a", "queued")), json.dumps(rec("b", "queued"))]
    path.write_text(f"{good[0]}\n{{not json\n{good[1]}\n")
    items = pipes.load_items(path)
    assert {"a", "b"} <= set(items)
    shown = ok(repo, "show")
    assert "1 unreadable" in shown


# ---------------------------------------------------------------- 2. step 1


def test_the_in_flight_item_nearest_to_merging_goes_first() -> None:
    choice = choose([rec("a", "review", pipe="lsm"), rec("b", "green"), rec("c", "queued")])
    assert label(choice) == "lsm/a"
    assert "nearest to merging" in choice.reason


def test_step_1_ties_go_to_the_higher_priority_pipe() -> None:
    choice = choose([rec("a", "review", pipe="lsm", age=5 * HOUR), rec("b", "review")])
    assert label(choice) == "core/b"


def test_step_1_ties_in_one_pipe_go_to_the_older_since() -> None:
    choice = choose([rec("young", "green", age=2 * HOUR), rec("old", "green", age=5 * HOUR)])
    assert label(choice) == "core/old"


# ---------------------------------------------------------------- 3. caps and busy


CAPPED = [rec("r", "red"), rec("g", "green"), rec("cq", "queued", age=2 * HOUR),
          rec("lq", "queued", pipe="lsm", age=HOUR)]  # fmt: skip


def test_a_pipe_with_two_in_flight_yields_to_another_pipe() -> None:
    choice = choose(CAPPED, busy={"worktree-r", "worktree-g"})
    assert pipes.in_flight(pipes.fold(CAPPED)) == 2
    assert label(choice) == "lsm/lq"


def test_a_not_busy_in_flight_item_is_taken_first() -> None:
    choice = choose(CAPPED, busy={"worktree-g"})
    assert label(choice) == "core/r"
    assert "nearest to merging" in choice.reason


def test_three_in_flight_all_busy_is_idle_never_queued_or_fallback() -> None:
    records = [*CAPPED, rec("d", "designing", pipe="lsm"),
               rec("fb", "queued", pipe="lsm", plain=True, fallback=True)]  # fmt: skip
    choice = choose(records, busy={"worktree-r", "worktree-g", "worktree-d"})
    assert pipes.in_flight(pipes.fold(records)) == 3
    assert label(choice) == "IDLE"


def test_the_only_pipe_with_runnable_items_may_hold_a_third() -> None:
    records = [rec("r", "red"), rec("g", "green"), rec("cq", "queued")]
    choice = choose(records, busy={"worktree-r", "worktree-g"})
    assert label(choice) == "core/cq"


def test_busy_branches_reads_the_note_names(tmp_path: Path) -> None:
    tasks = tmp_path / "current-task"
    tasks.mkdir()
    for name in ("developer-r-190000.md", "tester-tri-budget-rms-064037.md",
                 "reviewer-h31-x-190000.md", "session.md"):  # fmt: skip
        (tasks / name).write_text("x\n")
    busy = pipes.busy_branches(tasks)
    assert busy == {"worktree-r", "worktree-tri-budget-rms", "worktree-h31-x"}


def test_a_longer_worktree_name_does_not_make_its_prefix_busy(tmp_path: Path) -> None:
    tasks = tmp_path / "current-task"
    tasks.mkdir()
    (tasks / "developer-h31-x-190000.md").write_text("x\n")
    assert "worktree-h31" not in pipes.busy_branches(tasks)


# ---------------------------------------------------------------- 4. waiting


def test_waiting_items_are_not_in_flight() -> None:
    records = [rec("w1", "waiting"), rec("w2", "waiting"), rec("w3", "waiting", pipe="lsm"),
               rec("r", "red"), rec("q", "queued", pipe="lsm")]  # fmt: skip
    assert pipes.in_flight(pipes.fold(records)) == 1
    assert label(choose(records, busy={"worktree-r"})) == "lsm/q"


def test_a_waiting_item_is_never_chosen() -> None:
    assert label(choose([rec("w", "waiting", age=9 * DAY)])) == "IDLE"


# ---------------------------------------------------------------- 5. unattended


def test_unattended_skips_items_that_are_not_plain() -> None:
    records = [rec("r", "red"), rec("governed", "queued", age=2 * HOUR),
               rec("plain", "queued", pipe="lsm", plain=True)]  # fmt: skip
    assert label(choose(records, busy={"worktree-r"})) == "core/governed"
    assert label(choose(records, busy={"worktree-r"}, unattended=True)) == "lsm/plain"


def test_unattended_with_nothing_plain_is_idle(repo: Path) -> None:
    records = [rec("r", "red"), rec("governed", "queued")]
    assert label(choose(records, busy={"worktree-r"}, unattended=True)) == "IDLE"
    append(repo, *records)
    note(repo, "r")
    harness_fixtures.flag_on(repo)
    lines = ok(repo, "next").splitlines()
    assert IDLE in lines
    assert not any(line.startswith("next:") for line in lines)


# ---------------------------------------------------------------- 6. fallback


def test_the_fallback_is_taken_when_everything_else_waits() -> None:
    records = [rec("w1", "waiting"), rec("w2", "waiting", pipe="lsm"),
               rec("fb", "queued", pipe="lsm", plain=True, fallback=True)]  # fmt: skip
    choice = choose(records)
    assert label(choice) == "lsm/fb"
    assert "fallback" in choice.reason


def test_a_plain_queued_item_beats_the_fallback() -> None:
    records = [rec("w1", "waiting"), rec("w2", "waiting", pipe="lsm"),
               rec("fb", "queued", pipe="lsm", plain=True, fallback=True, age=9 * HOUR),
               rec("q", "queued", plain=True)]  # fmt: skip
    assert label(choose(records)) == "core/q"


def test_steps_2_and_3_skip_a_fallback_even_in_the_first_pipe() -> None:
    records = [rec("fb", "queued", plain=True, fallback=True),
               rec("q", "queued", pipe="lsm")]  # fmt: skip
    assert label(choose(records)) == "lsm/q"


def test_add_refuses_fallback_with_governed(repo: Path) -> None:
    result = cli(repo, "add", "lsm", "fb", "a fallback", "--fallback", "--governed")
    assert result.returncode != 0
    assert "--fallback" in result.stderr and "--governed" in result.stderr
    assert not any(
        line.get("id") == "fb" for line in harness_fixtures.jsonl(live_list(repo, "lsm"))
    )


# ---------------------------------------------------------------- 7. private definitions


def recording_run(repo: Path, calls: list[list[str]]) -> Any:
    """A `run` that records each git argument list, then runs it in `repo`."""

    def run(args: Any) -> tuple[int, str]:
        argv = [str(a) for a in args]
        calls.append(argv)
        command = argv if argv[0] == "git" else ["git", *argv]
        done = subprocess.run(command, cwd=repo, capture_output=True, text=True,
                              env=clean_env(), check=False)  # fmt: skip
        return done.returncode, done.stdout

    return run


@pytest.fixture
def with_private(repo: Path, tmp_path: Path) -> Path:
    """`repo` with a remote `private` whose master holds `lsm.toml` and `core.toml`."""
    source = tmp_path / "private-src"
    source.mkdir()
    git(source, "init", "-q", "-b", "master")
    git(source, "config", "user.email", "t@example.invalid")
    git(source, "config", "user.name", "T")
    (source / ".claude" / "pipes").mkdir(parents=True)
    (source / ".claude" / "pipes" / "core.toml").write_text(CORE_TOML)
    (source / ".claude" / "pipes" / "lsm.toml").write_text(LSM_TOML)
    git(source, "add", "-A")
    stamp = "2026-10-25T12:00:00+00:00"
    env = {**clean_env(), "GIT_COMMITTER_DATE": stamp, "GIT_AUTHOR_DATE": stamp}
    subprocess.run(["git", "-C", str(source), "commit", "-q", "-m", "pipes"], check=True,
                   capture_output=True, env=env)  # fmt: skip
    bare = tmp_path / "private.git"
    git(tmp_path, "clone", "-q", "--bare", str(source), str(bare))
    (repo / ".claude" / "pipes" / "lsm.toml").unlink()  # lsm is defined only on private
    git(repo, "remote", "add", "private", str(bare))
    return repo


def test_private_definitions_are_read_from_the_fetched_ref(with_private: Path) -> None:
    repo = with_private
    git(repo, "fetch", "-q", "private")
    sha7 = git(repo, "rev-parse", "--short", "private/master").strip()
    calls: list[list[str]] = []
    found, warnings = pipes.definitions(repo, run=recording_run(repo, calls))
    by_name = {p.name: p for p in found}
    assert sorted(by_name) == ["core", "lsm"]
    assert by_name["core"].source == "tree"
    assert by_name["lsm"].source == f"private/master {sha7}"
    assert warnings == []
    assert calls, "definitions reached git without the injected run"
    assert not any("fetch" in argv for argv in calls), calls
    header = pipes.show_lines(repo, NOW)[0]
    assert f"private/master {sha7}" in header
    assert "committed 2026-10-25" in header


def test_a_never_fetched_remote_warns_once_and_reads_the_tree(with_private: Path) -> None:
    repo = with_private
    calls: list[list[str]] = []
    found, warnings = pipes.definitions(repo, run=recording_run(repo, calls))
    assert [p.name for p in found] == ["core"]
    assert len(warnings) == 1
    assert "no private/master locally" in warnings[0]
    assert "git fetch private" in warnings[0]
    assert not any("fetch" in argv for argv in calls), calls
    result = cli(repo, "show")
    assert result.returncode == 0, result.stderr
    assert "no private/master locally" in result.stdout + result.stderr


def test_an_unknown_key_in_a_definition_is_warned(repo: Path) -> None:
    (repo / ".claude" / "pipes" / "core.toml").write_text(CORE_TOML + "checks = []\n")
    found, warnings = pipes.definitions(repo)
    assert "core" in [p.name for p in found]
    assert any("checks" in w for w in warnings), warnings


# ---------------------------------------------------------------- 8. stuck


STUCK_SET = [rec("old", "queued", pipe="lsm", age=4 * DAY), rec("fresh", "queued", age=DAY),
             rec("q", "waiting", age=5 * DAY, question="drop it?", default="keep")]  # fmt: skip


def test_the_stuck_threshold_is_three_days() -> None:
    assert timedelta(days=3) == pipes.STUCK_AFTER


def test_stuck_lists_queued_and_waiting_items_oldest_first() -> None:
    assert [i.id for i in pipes.stuck(pipes.fold(STUCK_SET), NOW)] == ["q", "old"]


@pytest.mark.parametrize(
    ("age", "is_stuck"),
    [(3 * DAY - timedelta(seconds=1), False), (3 * DAY, False),
     (3 * DAY + timedelta(seconds=1), True)],
)  # fmt: skip
def test_the_threshold_is_strict(age: timedelta, is_stuck: bool) -> None:
    items = pipes.fold([rec("x", "queued", age=age)])
    assert bool(pipes.stuck(items, NOW)) is is_stuck


def test_an_in_flight_item_is_never_stuck() -> None:
    items = pipes.fold([rec("x", "red", age=9 * DAY)])
    assert pipes.stuck(items, NOW) == []


def test_a_stuck_item_in_a_lower_pipe_goes_before_a_fresh_one() -> None:
    choice = choose(STUCK_SET)
    assert label(choice) == "lsm/old"
    assert "stuck 4 d" in choice.reason


def test_the_recap_block_lists_stuck_items_first_with_their_age(repo: Path) -> None:
    append(repo, *STUCK_SET)
    lines = pipes.show_lines(repo, NOW)
    assert lines[0].startswith("Pipes (in flight 0 of 3")
    assert "Stuck:" not in lines  # one table: no `Stuck:` heading, no per-pipe headings
    assert not any(line.startswith(("core (priority", "lsm (priority")) for line in lines)
    rows = table(lines)
    first, second = rows[0], rows[1]
    key = ("PIPE", "SINCE")
    assert [first[k] for k in key] == ["core", "5 d stuck"] and first["ITEM"].startswith("q ")
    assert [second[k] for k in key] == ["lsm", "4 d stuck"] and second["ITEM"].startswith("old ")
    assert (first["mark"], second["mark"]) == ("", "next")
    fresh = [row for row in rows[2:] if row["ITEM"].split()[0] == "fresh"]
    assert len(fresh) == 1 and fresh[0]["SINCE"] == "1 d", rows


# ---------------------------------------------------------------- 9. questions


def test_questions_number_the_stuck_one_first_each_with_its_default(repo: Path) -> None:
    append(
        repo,
        rec("new", "waiting", age=DAY, question="ship it?", default="wait"),
        rec("old", "waiting", pipe="lsm", age=5 * DAY, question="drop it?", default="keep"),
    )
    lines = after(pipes.show_lines(repo, NOW), "Questions for Ola:")
    first, second = lines[0].strip(), lines[1].strip()
    assert first.startswith("1.") and "lsm/old" in first and "drop it?" in first
    assert first.endswith("keep")
    assert second.startswith("2.") and "core/new" in second and second.endswith("wait")


def test_ask_sets_waiting_and_answer_restores_the_state(repo: Path) -> None:
    ok(repo, "add", "core", "q1", "a question item")
    ok(repo, "set", "q1", "designing", "--branch", "worktree-q1")
    ok(repo, "ask", "q1", "split it in two?", "--default", "keep it whole")
    shown = ok(repo, "show")
    asked = [line for line in after(shown.splitlines(), "Questions for Ola:") if "core/q1" in line]
    assert asked and asked[0].endswith("keep it whole")
    assert pipes.load_items(live_list(repo, "core"))["q1"].state == "waiting"
    ok(repo, "answer", "q1", "keep it whole")
    assert pipes.load_items(live_list(repo, "core"))["q1"].state == "designing"
    assert "split it in two?" not in ok(repo, "show")


def test_answer_records_the_ruling_on_the_item_record(repo: Path) -> None:
    ok(repo, "add", "core", "r1", "an item with a question")
    ok(repo, "ask", "r1", "split it in two?", "--default", "keep it whole")
    ok(repo, "answer", "r1", "split it, the reader part first")
    records = [json.loads(line) for line in live_list(repo, "core").read_text().splitlines()]
    newest = [r for r in records if r["id"] == "r1"][-1]
    assert newest.get("ruling") == "split it, the reader part first", newest


def test_answer_without_an_open_question_is_refused_and_changes_nothing(repo: Path) -> None:
    ok(repo, "add", "core", "h1", "an item with no question")
    ok(repo, "set", "h1", "red", "--branch", "worktree-h1")
    lines_before = live_list(repo, "core").read_text()
    result = cli(repo, "answer", "h1", "stray")
    assert result.returncode != 0, result.stdout
    assert "h1" in result.stderr and "no open question" in result.stderr, result.stderr
    assert live_list(repo, "core").read_text() == lines_before
    assert pipes.load_items(live_list(repo, "core"))["h1"].state == "red"


def test_a_question_on_a_pipe_changes_no_item(repo: Path) -> None:
    ok(repo, "add", "core", "k1", "an item")
    ok(repo, "ask", "core", "move the queue by hand?", "--default", "yes")
    asked = after(ok(repo, "show").splitlines(), "Questions for Ola:")
    assert any("move the queue by hand?" in line and line.endswith("yes") for line in asked)
    assert pipes.load_items(live_list(repo, "core"))["k1"].state == "queued"


# ---------------------------------------------------------------- the commands


def test_next_prints_one_parseable_line_and_a_busy_line_per_skipped_item(repo: Path) -> None:
    ok(repo, "add", "core", "a", "first")
    ok(repo, "set", "a", "red", "--branch", "worktree-a")
    ok(repo, "add", "core", "b", "second")
    held = note(repo, "a")
    lines = ok(repo, "next").splitlines()
    chosen = [line for line in lines if line.startswith("next: ")]
    assert len(chosen) == 1, lines
    assert chosen[0].startswith("next: core/b queued (") and chosen[0].endswith(")")
    busy = [line for line in lines if line.startswith("busy: ")]
    assert len(busy) == 1 and busy[0].startswith("busy: core/a (") and held.name in busy[0]
    assert len(lines) == len(chosen) + len(busy), lines


def test_top_moves_an_item_to_the_head_of_its_pipe(repo: Path) -> None:
    ok(repo, "add", "core", "x", "added first")
    ok(repo, "add", "core", "y", "added second")
    assert ok(repo, "next").splitlines()[-1].startswith("next: core/x ")
    ok(repo, "top", "y")
    assert ok(repo, "next").splitlines()[-1].startswith("next: core/y ")


def test_an_id_is_unique_across_pipes(repo: Path) -> None:
    ok(repo, "add", "core", "dup1", "in core")
    result = cli(repo, "add", "lsm", "dup1", "in lsm too")
    assert result.returncode != 0
    assert "dup1" in result.stderr
    assert not any(
        line.get("id") == "dup1" for line in harness_fixtures.jsonl(live_list(repo, "lsm"))
    )


def test_the_block_is_capped_at_40_lines_of_120_characters(repo: Path) -> None:
    append(repo, *[rec(f"i{n:02d}", "queued", title="t" * 200) for n in range(60)])
    lines = pipes.show_lines(repo, NOW)
    assert len(lines) <= 41
    assert lines[-1].strip().startswith("...") and lines[-1].strip().endswith("more")
    assert max(len(line) for line in lines) <= 120


# ---------------------------------------------------------------- 10. the recap


def test_the_recap_prints_the_pipes_after_waiting_on_ola(repo: Path) -> None:
    ok(repo, "add", "core", "r", "in flight")
    ok(repo, "set", "r", "red", "--branch", "worktree-r")
    ok(repo, "add", "core", "q", "queued")
    result = run_script(repo, "tools/session_state.py", "")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    waiting = first_index(lines, "Waiting on Ola:")
    block = first_index(lines, "Pipes (in flight 1 of 3")
    assert waiting < block < first_index(lines, "Next on ROADMAP.md:")


# ---------------------------------------------------------------- 11. the guard


def load_guard() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "h31_guard_governance", REAL / ".claude" / "hooks" / "guard_governance.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    "command",
    ['python3 tools/work_pipes.py set h31 red',
     'python3 tools/work_pipes.py add core x "x" --plain',
     'python3 tools/work_pipes.py ask h31 "q?" --default "yes"',
     "echo '{}' >> .git/pipes/core.jsonl"],
)  # fmt: skip
def test_the_commands_pass_the_governance_guard(command: str) -> None:
    assert load_guard().judge_bash(command) is None


def test_work_pipes_py_is_governed() -> None:
    assert load_guard().governed("tools/work_pipes.py") is True


# ---------------------------------------------------------------- 12. the table


def live(item: str, state: str, age: timedelta, **extra: Any) -> dict[str, Any]:
    """`rec` aged from the real clock, for commands that read it.

    Ages sit 10 minutes past a whole hour, so neither rounding nor flooring to
    hours moves across a boundary while the test runs.
    """
    record = rec(item, state, **extra)
    since = datetime.now(UTC).replace(microsecond=0) - age - timedelta(minutes=10)
    record["since"] = since.isoformat(timespec="seconds")
    return record


def table_set() -> list[dict[str, Any]]:
    """Test 12's items: four open in `core`, one in `lsm`, and a merged one."""
    return [
        live("r", "red", 2 * HOUR), live("g", "green", 3 * HOUR),
        live("w", "waiting", 5 * HOUR, question="ship it?", default="wait"),
        live("q", "queued", 2 * HOUR, title="t" * 200),
        live("l", "queued", DAY + 2 * HOUR, pipe="lsm"), live("m", "merged", HOUR, pr=7),
    ]  # fmt: skip


def test_the_table_has_one_aligned_row_per_open_item(repo: Path) -> None:
    append(repo, *table_set())
    note(repo, "r")  # `r` is busy: its worktree's note is in current-task
    shown = ok(repo, "show")
    assert ok(repo) == shown, "work_pipes.py with no arguments must print `show`"
    lines = shown.splitlines()
    assert max(len(line) for line in lines) <= 120
    rows = table(lines)
    by_id = {row["ITEM"].split()[0]: row for row in rows}
    assert [row["ITEM"].split()[0] for row in rows] == ["g", "r", "w", "q", "l"]
    assert len(by_id) == len(rows) and "m" not in by_id
    chosen = [line for line in ok(repo, "next").splitlines() if line.startswith("next: ")]
    assert chosen and chosen[0].startswith("next: core/g ")
    assert [row["ITEM"].split()[0] for row in rows if row["mark"] == "next"] == ["g"]
    assert {i: by_id[i]["WAITING ON"] for i in "rgwql"} == {
        "r": "its agent", "g": "free agent", "w": "Ola", "q": "free agent", "l": "free agent"
    }  # fmt: skip
    assert {i: (by_id[i]["PIPE"], by_id[i]["PRI"], by_id[i]["STEP"]) for i in "rwl"} == {
        "r": ("core", "1", "red"), "w": ("core", "1", "waiting"), "l": ("lsm", "2", "queued")
    }  # fmt: skip
    assert by_id["w"]["NEXT"].startswith("question 1 below")
    assert (by_id["r"]["SINCE"], by_id["l"]["SINCE"]) == ("2 h", "1 d")
    assert not any("stuck" in row["SINCE"] for row in rows)
    assert by_id["q"]["ITEM"].endswith("...") and len(by_id["q"]["ITEM"]) <= 30


def test_unattended_an_item_that_is_not_plain_waits_on_ola_present(repo: Path) -> None:
    append(repo, live("gov", "queued", HOUR), live("pl", "queued", HOUR, plain=True),
           live("w", "waiting", HOUR, question="ship it?", default="wait"))  # fmt: skip
    harness_fixtures.flag_on(repo)
    rows = table(ok(repo).splitlines())
    waiting_on = {row["ITEM"].split()[0]: row["WAITING ON"] for row in rows}
    assert waiting_on == {"gov": "Ola present", "pl": "free agent", "w": "Ola"}
