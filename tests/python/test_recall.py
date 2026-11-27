"""`tools/recall.py`: what is already known on a topic, and the gate `questions.py add` keeps.

Every test runs `command()` on a git repository, a memory folder and a search
record in `tmp_path`, so neither the private repository, the main session's
memories nor the live record under the git common directory is read or written.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

import harness_fixtures

pytestmark = pytest.mark.harness

recall = harness_fixtures.Tool("recall")
findings = harness_fixtures.Tool("findings")
questions = harness_fixtures.Tool("questions")

NOW = datetime(2026, 11, 27, 6, 0, tzinfo=UTC)
#: Question mode began ten hours before `NOW`: no clock, only the mode's start.
BEGAN = NOW - timedelta(hours=10)

LEDGER = [
    {"text": "Sliver triangles are 4 % of the mesh.", "commit": "a" * 40,
     "date": "2026-11-20", "topics": ["sliver"], "ola": "Yes"},
    {"text": "The breakline adds nothing to the refine time.", "commit": "b" * 40,
     "date": "2026-11-21", "topics": ["breakline", "mesh"], "ola": "Yes"},
    {"text": "The clip step is too slow.", "commit": "c" * 40,
     "date": "2026-11-22", "topics": ["clip"], "ola": "Yes"},
]  # fmt: skip

SLIVER_MEMORY = """---
name: breakline-sliver-idea
description: "Breaklines so slivers split (Ola 2026-11-23)"
metadata:
  modified: 2026-11-26T05:58:30.170Z
---

Ola, 2026-11-23: "A SLIVER must never stay in the mesh." Then more.
Ola 2026-11-26, on the breakline: "Keep it straight." A second line.
"""

OTHER_MEMORY = """---
name: plain-writing
description: "Ola finds the mannered register hard to parse"
metadata:
  modified: 2026-10-01T00:00:00Z
---

Nothing about mesh here.
"""

BENCHMARKS = {
    "docs/benchmarks/2026-11-12/thin-a/README.md":
        "# Sliver probe\n\nSome text about mesh slivers.\n\n**Verdict: ACCEPTED.** It moves.\n",
    "docs/benchmarks/2026-11-13/thin-b/README.md":
        "# Second sliver probe\n\nMesh text.\n\n## Results\n\nSliver adds 2 %.\n",
    "docs/benchmarks/2026-11-14/thin-c/README.md":
        "# Third\n\n## Verdict: REJECTED\n\nmesh sliver\n",
    "docs/benchmarks/2026-11-14/thin-c/deeper/README.md": "**mesh sliver deeper**\n",
    "docs/benchmarks/2026-11-15/other/README.md": "**Clip only.**\n",
}  # fmt: skip


def git(where: Path, *argv: str) -> str:
    return subprocess.run(["git", "-C", str(where), "-c", "user.name=t", "-c", "user.email=t@t",
                           "-c", "commit.gpgsign=false", *argv], env=harness_fixtures.clean_env(),
                          capture_output=True, text=True, check=True).stdout  # fmt: skip


@pytest.fixture
def where(tmp_path: Path) -> dict[str, Path]:
    """A repository whose `master` holds the ledger and the benchmark records, a memory
    folder, and the search record, not yet written."""
    repo = tmp_path / "private"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "master")
    (repo / "docs").mkdir()
    (repo / "docs" / "findings.jsonl").write_text("".join(json.dumps(e) + "\n" for e in LEDGER))
    for path, text in BENCHMARKS.items():
        (repo / path).parent.mkdir(parents=True, exist_ok=True)
        (repo / path).write_text(text)
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "records")
    memory = tmp_path / "memory"
    memory.mkdir()
    (memory / "breakline-sliver-idea.md").write_text(SLIVER_MEMORY)
    (memory / "plain-writing.md").write_text(OTHER_MEMORY)
    searches = tmp_path / "common" / "findings-searches.jsonl"
    harness_fixtures.question_mode(searches.parent, BEGAN)
    return {"repo": repo, "memory": memory, "searches": searches}


def run(where: dict[str, Path], *words: str, limit: int = 10, ledger: Path | None = None,
        now: datetime = NOW) -> list[str]:  # fmt: skip
    return recall.command(list(words), repo=where["repo"], ref="master", ledger=ledger,
                          memory=where["memory"], searches=where["searches"], now=now,
                          limit=limit)  # fmt: skip


def section(lines: list[str], head: str) -> list[str]:
    """The lines of the section whose heading starts with `head`, heading excluded."""
    start = next(n for n, line in enumerate(lines) if line.startswith(head))
    rest = lines[start + 1 :]
    end = next((n for n, line in enumerate(rest) if not line.startswith(" ")), len(rest))
    return rest[:end]


def test_the_ledger_needs_every_word_in_any_case_and_keeps_its_numbers(where: dict) -> None:
    lines = section(run(where, "MESH", "sliver"), "Ledger")
    assert len(lines) == 1
    assert lines[0].strip().startswith("1. Sliver triangles")
    assert "2026-11-20" in lines[0]
    topical = section(run(where, "mesh", "breakline"), "Ledger")
    assert len(topical) == 1  # matched by its topic
    assert topical[0].strip().startswith("2. The breakline")


def test_findings_ledger_names_a_ledger_file(where: dict, tmp_path: Path) -> None:
    own = tmp_path / "own.jsonl"
    own.write_text(json.dumps({**LEDGER[2], "text": "Own clip line."}) + "\n")
    assert "Own clip line." in "\n".join(section(run(where, "clip", ledger=own), "Ledger"))


def test_a_memory_prints_its_path_date_description_and_dated_ola_lines(where: dict) -> None:
    lines = section(run(where, "sliver", "mesh"), "Memories")
    path = str(where["memory"] / "breakline-sliver-idea.md")
    assert lines[0].strip() == (
        f"1. {path}, 2026-11-26: breakline-sliver-idea: Breaklines so slivers split "
        "(Ola 2026-11-23)"
    )
    assert lines[1].strip() == 'Ola, 2026-11-23: "A SLIVER must never stay in the mesh."'
    assert lines[2].strip() == 'Ola, 2026-11-26: "Keep it straight."'
    assert len(lines) == 3  # plain-writing has "mesh" but not "sliver"


def test_benchmark_verdicts_take_the_first_bold_or_verdict_line(where: dict) -> None:
    lines = [line.strip() for line in section(run(where, "mesh", "sliver"), "Benchmark")]
    assert lines == [
        "1. docs/benchmarks/2026-11-14/thin-c/README.md, 2026-11-14: Verdict: REJECTED",
        "2. docs/benchmarks/2026-11-13/thin-b/README.md, 2026-11-13: Sliver adds 2 %.",
        "3. docs/benchmarks/2026-11-12/thin-a/README.md, 2026-11-12: Verdict: ACCEPTED. It moves.",
    ]  # newest first; a README deeper than <date>/<name>/ is no record


def test_a_cap_prints_how_many_more_matched(where: dict) -> None:
    lines = [line.strip() for line in section(run(where, "sliver", limit=1), "Benchmark")]
    assert len(lines) == 2 and lines[1] == "+2 more: rerun with --limit 3"


def test_no_match_is_said_per_source(where: dict) -> None:
    lines = run(where, "quadtree")
    for head in ("Ledger", "Memories", "Benchmark"):
        assert [line.strip() for line in section(lines, head)] == ["no match"]


def test_a_run_is_recorded_and_counts_as_a_ledger_search(where: dict) -> None:
    run(where, "Sliver", "breakline")
    record = json.loads(where["searches"].read_text().splitlines()[-1])
    assert record == {"at": NOW.isoformat(), "words": ["Sliver", "breakline"], "tool": "recall"}
    assert findings.unsearched(where["searches"], NOW) is None


# ---------------------------------------------------------------- the gate


def ask(store: Path, text: str, *extra: str) -> list[str]:
    args = questions.parser().parse_args(["add", text, "--default", "wait", *extra])
    return list(questions.command(args, store, NOW))


def test_questions_add_is_refused_until_recall_ran(tmp_path: Path) -> None:
    store = tmp_path / "common" / "questions.jsonl"
    harness_fixtures.question_mode(store.parent, BEGAN)
    harness_fixtures.searched(store.parent, NOW)  # a ledger search alone is not enough
    with pytest.raises(questions.RefusedError) as caught:
        ask(store, "Shall a sliver leave the domain?")
    assert str(caught.value) == (
        "recall what is already known on this topic first: "
        "python3 tools/recall.py <a word of the question or ledger line>"
    )
    assert not store.exists()


def test_outside_question_mode_no_recall_counts(where: dict) -> None:
    store = where["searches"].with_name("questions.jsonl")
    (where["searches"].parent / "harness" / "question-mode.json").unlink()
    run(where, "sliver")
    with pytest.raises(questions.RefusedError) as caught:
        ask(store, "Shall a sliver leave the domain?")
    assert str(caught.value) == "enter question mode first: python3 tools/questions.py mode enter"
    assert not store.exists()


def test_recall_on_other_words_does_not_open_the_gate(where: dict, tmp_path: Path) -> None:
    store = where["searches"].with_name("questions.jsonl")
    run(where, "clip")
    with pytest.raises(questions.RefusedError):
        ask(store, "Shall a sliver leave the domain?")


def test_a_recall_counts_from_the_start_of_the_mode_however_old(where: dict) -> None:
    store = where["searches"].with_name("questions.jsonl")
    run(where, "clip", now=BEGAN + timedelta(minutes=1))
    run(where, "Sliver", "hook", now=BEGAN + timedelta(minutes=1))  # nine hours before NOW
    assert ask(store, "Shall slivers leave the mesh?") == ["stored as q1"]
    where["searches"].unlink()
    run(where, "sliver", now=BEGAN - timedelta(minutes=1))  # before the mode began
    with pytest.raises(questions.RefusedError):
        ask(store, "Shall a sliver leave the breakline?")


def test_a_word_in_the_context_counts_and_small_words_do_not(where: dict) -> None:
    store = where["searches"].with_name("questions.jsonl")
    run(where, "the", "breakline")
    stored = ask(store, "Shall it go on?", "--context", "The breakline lets slivers split.")
    assert stored == ["stored as q1"]
    where["searches"].unlink()
    run(where, "the", "shall")
    with pytest.raises(questions.RefusedError):
        ask(store, "Shall the sliver go?")
