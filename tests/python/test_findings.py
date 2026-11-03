"""`tools/findings.py`: the ledger of settled results, and the search gate.

Every test runs `command()` on a ledger and a search record in `tmp_path`, so
the repository's ledger and the live record under the git common directory are
never read or written. The commit a test names is this checkout's HEAD.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

import harness_fixtures

pytestmark = pytest.mark.harness

findings = harness_fixtures.Tool("findings")
questions = harness_fixtures.Tool("questions")

NOW = datetime(2026, 11, 3, 6, 0, tzinfo=UTC)
HEAD = subprocess.run(["git", "-C", str(harness_fixtures.REAL), "rev-parse", "HEAD"],
                      capture_output=True, text=True, check=True).stdout.strip()  # fmt: skip


@pytest.fixture
def paths(tmp_path: Path) -> tuple[Path, Path]:
    """The ledger and the search record, neither yet written."""
    return tmp_path / "ledger.jsonl", tmp_path / "common" / "findings-searches.jsonl"


def run(paths: tuple[Path, Path], *argv: str, now: datetime = NOW) -> list[str]:
    args = findings.parser().parse_args(list(argv))
    return list(findings.command(args, *paths, now))


def refused(paths: tuple[Path, Path], *argv: str) -> str:
    with pytest.raises(findings.RefusedError) as caught:
        run(paths, *argv)
    assert not paths[0].exists()
    return str(caught.value)


def add(paths: tuple[Path, Path], text: str, *topics: str) -> None:
    run(paths, "add", text, "--commit", HEAD, *(f"--topic={t}" for t in topics))


def test_add_stores_the_result_with_its_full_commit_and_topics(paths: tuple[Path, Path]) -> None:
    add(paths, "Two delays in series slow the response twice over.", "delay", "response")
    [entry] = [json.loads(line) for line in paths[0].read_text().splitlines()]
    assert entry["text"] == "Two delays in series slow the response twice over."
    assert entry["commit"] == HEAD and entry["topics"] == ["delay", "response"]


def test_add_refuses_a_commit_git_does_not_have(paths: tuple[Path, Path]) -> None:
    assert "no commit" in refused(paths, "add", "A result.", "--commit", "0" * 40, "--topic=x")


def test_add_needs_a_topic(paths: tuple[Path, Path]) -> None:
    with pytest.raises(SystemExit):
        run(paths, "add", "A result.", "--commit", HEAD)


def test_add_refuses_text_that_is_not_plain_words(paths: tuple[Path, Path]) -> None:
    said = refused(paths, "add", "See tools/brief.py for it.", "--commit", HEAD, "--topic=x")
    assert "plain words" in said


def test_add_refuses_a_text_over_the_cap(paths: tuple[Path, Path]) -> None:
    said = refused(paths, "add", "word " * 60, "--commit", HEAD, "--topic=x")
    assert "at most" in said


def test_search_prints_matches_numbered_with_commit_and_date(paths: tuple[Path, Path]) -> None:
    add(paths, "The slope term doubles the delay.", "delay")
    add(paths, "Coarse cells hide the ridges.", "cells")
    add(paths, "A second delay sits in the routing.", "routing")
    lines = run(paths, "search", "delay")
    assert len(lines) == 2
    assert lines[0].startswith("1. The slope term doubles the delay.")
    assert lines[1].startswith("2. A second delay sits")
    assert all(HEAD[:8] in line and "20" in line for line in lines)


def test_search_matches_a_topic_and_ignores_case(paths: tuple[Path, Path]) -> None:
    add(paths, "Coarse cells hide the ridges.", "Resolution")
    assert len(run(paths, "search", "resolution")) == 1


def test_search_with_no_ledger_says_nothing_matched_and_records_the_search(
    paths: tuple[Path, Path],
) -> None:
    [line] = run(paths, "search", "slope", "delay")
    assert "no settled result" in line.lower()
    record = json.loads(paths[1].read_text().splitlines()[-1])
    assert record["words"] == ["slope", "delay"] and record["at"] == NOW.isoformat()


def test_list_prints_every_entry(paths: tuple[Path, Path]) -> None:
    add(paths, "First result.", "a")
    add(paths, "Second result.", "b")
    assert [line[:3] for line in run(paths, "list")] == ["1. ", "2. "]


def test_the_gate_refuses_until_a_search_and_names_the_command(paths: tuple[Path, Path]) -> None:
    said = findings.unsearched(paths[1], NOW)
    assert said is not None and "python3 tools/findings.py search" in said
    run(paths, "search", "anything")
    assert findings.unsearched(paths[1], NOW + timedelta(minutes=59)) is None
    assert findings.unsearched(paths[1], NOW + timedelta(minutes=61)) is not None


def test_questions_add_is_refused_without_a_recent_search(tmp_path: Path) -> None:
    store = tmp_path / "questions.jsonl"
    args = questions.parser().parse_args(["add", "Shall we go on?", "--default", "go on"])
    with pytest.raises(questions.RefusedError) as caught:
        questions.command(args, store, NOW)
    assert "findings.py search" in str(caught.value)
    assert not store.exists()


def test_questions_act_and_taken_need_no_search(tmp_path: Path) -> None:
    store = tmp_path / "questions.jsonl"
    for argv in (["add", "May I push?", "--default", "no", "--act"],
                 ["taken", "Kept the old name.", "--default", "keep"]):  # fmt: skip
        questions.command(questions.parser().parse_args(argv), store, NOW)
