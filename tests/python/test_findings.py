"""`tools/findings.py`: the ledger of settled results, and the search gate.

Every test runs `command()` on a ledger and a search record in `tmp_path`, so
the repository's ledger and the live record under the git common directory are
never read or written. The commit a test names is this checkout's HEAD, which
the ledger's own repository (a `git init` in `tmp_path`) borrows through git's
alternates file, as a private repository holding the same commits would.
"""

from __future__ import annotations

import importlib
import json
import subprocess
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

import harness_fixtures
from brief_fixtures import human

pytestmark = pytest.mark.harness

findings = harness_fixtures.Tool("findings")
questions = harness_fixtures.Tool("questions")

NOW = datetime(2026, 11, 3, 6, 0, tzinfo=UTC)
HEAD = subprocess.run(["git", "-C", str(harness_fixtures.REAL), "rev-parse", "HEAD"],
                      capture_output=True, text=True, check=True).stdout.strip()  # fmt: skip


AGREED = "Yes, say it that way."
STAMP = "2026-11-03T05:50:00.000Z"
SESSION = "0f1e2d3c-4b5a-6978-8796-a5b4c3d2e1f0"


@pytest.fixture
def paths(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """The ledger and the search record, neither yet written, and a session transcript
    in which Ola said `AGREED`, read by `brief.py`'s quotation check."""
    assert findings.LIMIT  # loads the tool, which puts tools/ on the path
    transcripts = tmp_path / "transcripts"
    transcripts.mkdir()
    (transcripts / f"{SESSION}.jsonl").write_text(json.dumps(human(AGREED, STAMP)) + "\n")
    monkeypatch.setattr(importlib.import_module("session_state"), "TRANSCRIPTS", transcripts)
    monkeypatch.setenv("CLAUDE_CODE_SESSION_ID", SESSION)
    git(tmp_path, "init", "-q")
    common = git(harness_fixtures.REAL, "rev-parse", "--path-format=absolute", "--git-common-dir")
    (tmp_path / ".git" / "objects" / "info" / "alternates").write_text(f"{common}/objects\n")
    return tmp_path / "ledger.jsonl", tmp_path / "common" / "findings-searches.jsonl"


def git(where: Path, *argv: str) -> str:
    return subprocess.run(["git", "-C", str(where), "-c", "user.name=t", "-c", "user.email=t@t",
                           "-c", "commit.gpgsign=false", *argv],
                          capture_output=True, text=True, check=True).stdout.strip()  # fmt: skip


def run(paths: tuple[Path, Path], *argv: str, now: datetime = NOW) -> list[str]:
    args = findings.parser().parse_args(list(argv))
    return list(findings.command(args, *paths, now))


def refused(paths: tuple[Path, Path], *argv: str) -> str:
    with pytest.raises(findings.RefusedError) as caught:
        run(paths, *argv)
    assert not paths[0].exists()
    return str(caught.value)


def add(paths: tuple[Path, Path], text: str, *topics: str) -> None:
    run(paths, "add", text, "--commit", HEAD, "--ola", AGREED, *(f"--topic={t}" for t in topics))


def entries(paths: tuple[Path, Path]) -> list[dict[str, object]]:
    return [json.loads(line) for line in paths[0].read_text().splitlines()]


def test_add_stores_the_result_with_its_full_commit_and_topics(paths: tuple[Path, Path]) -> None:
    add(paths, "Two delays in series slow the response twice over.", "delay", "response")
    [entry] = entries(paths)
    assert entry["text"] == "Two delays in series slow the response twice over."
    assert entry["commit"] == HEAD and entry["topics"] == ["delay", "response"]
    assert entry["ola"] == f'"{AGREED}" ({STAMP})'


def test_add_needs_olas_words_from_this_sessions_transcript(paths: tuple[Path, Path]) -> None:
    with pytest.raises(SystemExit):
        run(paths, "add", "A result.", "--commit", HEAD, "--topic=x")
    said = refused(paths, "add", "A result.", "--commit", HEAD, "--topic=x", "--ola", "Sure.")
    assert "transcript" in said


def test_add_refuses_a_commit_git_does_not_have(paths: tuple[Path, Path]) -> None:
    assert "no commit" in refused(
        paths, "add", "A result.", "--commit", "0" * 40, "--topic=x", "--ola", AGREED
    )


def test_add_needs_a_topic(paths: tuple[Path, Path]) -> None:
    with pytest.raises(SystemExit):
        run(paths, "add", "A result.", "--commit", HEAD)


def test_add_refuses_text_that_is_not_plain_words(paths: tuple[Path, Path]) -> None:
    said = refused(
        paths, "add", "See tools/brief.py for it.", "--commit", HEAD, "--topic=x", "--ola", AGREED
    )
    assert "plain words" in said


def test_add_refuses_a_text_over_the_cap(paths: tuple[Path, Path]) -> None:
    said = refused(paths, "add", "word " * 60, "--commit", HEAD, "--topic=x", "--ola", AGREED)
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


def test_replace_rewrites_an_entry_in_place_with_olas_words(paths: tuple[Path, Path]) -> None:
    add(paths, "First result.", "a")
    add(paths, "Second result.", "b")
    with pytest.raises(findings.RefusedError, match="transcript"):
        run(paths, "replace", "1", "New first.", "--ola", "Sure.")
    run(paths, "replace", "1", "New first.", "--ola", AGREED)
    assert [e["text"] for e in entries(paths)] == ["New first.", "Second result."]
    assert entries(paths)[0]["topics"] == ["a"]


def test_drop_removes_an_entry_and_a_missing_number_is_refused(paths: tuple[Path, Path]) -> None:
    add(paths, "First result.", "a")
    add(paths, "Second result.", "b")
    run(paths, "drop", "1")
    assert [e["text"] for e in entries(paths)] == ["Second result."]
    with pytest.raises(findings.RefusedError):
        run(paths, "drop", "2")


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


def main(capsys: pytest.CaptureFixture[str], *argv: str) -> list[str]:
    assert findings.main(list(argv)) == 0
    return capsys.readouterr().out.splitlines()


@pytest.fixture
def live(paths: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """`main()` with its search record moved to `tmp_path`, and the default ledger too;
    patched in the loaded modules themselves, since `Tool` only forwards reads."""
    record = paths[1].with_name("q.jsonl")
    monkeypatch.setattr(importlib.import_module("questions"), "store_of", lambda _root: record)
    monkeypatch.setattr(importlib.import_module("u1_findings"), "LEDGER", paths[0].with_name("d"))
    monkeypatch.delenv("FINDINGS_LEDGER", raising=False)
    return paths


def test_the_variable_redirects_add_search_and_list(
    live: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("FINDINGS_LEDGER", str(live[0]))
    main(capsys, "add", "A private result.", "--commit", HEAD, "--topic=x", "--ola", AGREED)
    assert [e["text"] for e in entries(live)] == ["A private result."]
    assert main(capsys, "search", "private")[0].startswith("1. A private result.")
    assert main(capsys, "list")[0].startswith("1. A private result.")
    assert live[1].exists() and not findings.LEDGER.exists()


def test_unset_the_ledger_stays_the_default(
    live: tuple[Path, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    main(capsys, "add", "A public result.", "--commit", HEAD, "--topic=x", "--ola", AGREED)
    assert main(capsys, "list")[0].startswith("1. A public result.")
    assert findings.LEDGER.exists() and not live[0].exists()


def test_the_commit_is_checked_in_the_ledgers_own_repository(paths: tuple[Path, Path]) -> None:
    private = paths[0].parent / "private"
    git(private.parent, "init", "-q", str(private))
    git(private, "commit", "-q", "--allow-empty", "-m", "a private commit")
    own = git(private, "rev-parse", "HEAD")
    ledger = (private / "docs" / "findings.jsonl", paths[1])
    assert "no commit" in refused(
        ledger, "add", "A.", "--commit", HEAD, "--topic=x", "--ola", AGREED
    )
    run(ledger, "add", "A.", "--commit", own, "--topic=x", "--ola", AGREED)
    assert [e["commit"] for e in entries(ledger)] == [own]


LINE = "Two delays in series slow the response twice over."


def asked(paths: tuple[Path, Path], text: str, words: str | None) -> str:
    """Store a question in the store beside the search record, answered with `words`."""
    record = {"id": "q7", "kind": "question", "text": text, "default": "keep it",
              "choices": [], "item": None, "source": "main", "state": "open", "words": None,
              "at": STAMP}  # fmt: skip
    store = paths[1].with_name("questions.jsonl")
    lines = [record] + ([{**record, "state": "answered", "words": words}] if words else [])
    store.parent.mkdir(parents=True, exist_ok=True)
    store.write_text("".join(json.dumps(r) + "\n" for r in lines))
    return "q7"


def by_question(paths: tuple[Path, Path], qid: str, text: str = LINE) -> list[str]:
    return run(paths, "add", text, "--commit", HEAD, "--topic=delay", "--question", qid)


def test_an_answered_agreeing_question_lets_the_line_in(paths: tuple[Path, Path]) -> None:
    qid = asked(paths, f'May the ledger say: "{LINE}"?', "Yes")
    by_question(paths, qid)
    [entry] = entries(paths)
    assert entry["text"] == LINE and "q7" in str(entry["ola"]) and "Yes" in str(entry["ola"])


def test_an_unanswered_question_is_refused(paths: tuple[Path, Path]) -> None:
    qid = asked(paths, f'May the ledger say: "{LINE}"?', None)
    with pytest.raises(findings.RefusedError, match="not answered"):
        by_question(paths, qid)
    assert not paths[0].exists()


@pytest.mark.parametrize("words", ["Change it", "No", "Yes, but say three delays"])
def test_an_answer_other_than_agreement_is_refused(paths: tuple[Path, Path], words: str) -> None:
    qid = asked(paths, f'May the ledger say: "{LINE}"?', words)
    with pytest.raises(findings.RefusedError, match="agree"):
        by_question(paths, qid)
    assert not paths[0].exists()


def test_a_line_the_question_does_not_hold_is_refused(paths: tuple[Path, Path]) -> None:
    qid = asked(paths, f'May the ledger say: "{LINE}"?', "Yes")
    with pytest.raises(findings.RefusedError, match="does not hold"):
        by_question(paths, qid, "Two delays in series slow the response.")
    assert not paths[0].exists()


def test_an_unknown_question_is_refused(paths: tuple[Path, Path]) -> None:
    with pytest.raises(findings.RefusedError, match="no question"):
        by_question(paths, "q9")


def test_replace_takes_an_answered_question_too(paths: tuple[Path, Path]) -> None:
    add(paths, "First result.", "a")
    qid = asked(paths, f'May the ledger say: "{LINE}"?', "Agreed.")
    run(paths, "replace", "1", LINE, "--question", qid)
    assert [e["text"] for e in entries(paths)] == [LINE]


def test_ola_and_question_together_are_refused(paths: tuple[Path, Path]) -> None:
    with pytest.raises(SystemExit):
        run(paths, "add", LINE, "--commit", HEAD, "--topic=x", "--ola", AGREED, "--question", "q7")
