"""`.claude/hooks/record_review.py`: the verdict copied from `@reviewer`'s handback (h20 §4.5).

The spec is `docs/increments/h20-away-push-grants.md` §3.3. Each event is a
`PreToolUse` event for `SubagentHandback` with `tool_input.message` set,
unless a test says otherwise, run through a copy of the hook in a temporary
repository, so `reviews.jsonl` is `<tmp>/.git/harness/reviews.jsonl`. Every
case prints nothing and exits 0; the expected result is the lines of
`reviews.jsonl`.

Pinned here beyond the design's words: a recorded line's `event` is the
event's `hook_event_name`.
"""

from __future__ import annotations

import os
import stat
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import (
    RECORD_REVIEW,
    jsonl,
    make_repo,
    reviews_path,
    run_script,
    set_mode,
    state_dir,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

SHA = "0123456789abcdef0123456789abcdef01234567"
OTHER = "fedcba9876543210fedcba9876543210fedcba98"
HANDBACK = f"Result: two findings fixed.\n\nQuestions for Ola: none\n\nReview: APPROVED {SHA}"


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


def handback_event(repo: Path, message: str, **who: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": "SubagentHandback",
        "tool_input": {"message": message},
        "session_id": "s1",
        "cwd": str(repo),
        **who,
    }


def stop_event(repo: Path, message: str, **who: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "SubagentStop",
        "last_assistant_message": message,
        "session_id": "s1",
        "cwd": str(repo),
        **who,
    }


REVIEWER = {"agent_type": "reviewer", "agent_id": "r1"}


def recorded(repo: Path, event: dict[str, Any]) -> list[dict[str, Any]]:
    """Run the hook on `event`; it must print nothing and exit 0. The lines of reviews.jsonl."""
    result = run_script(repo, RECORD_REVIEW, event)
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    assert result.stdout == "", result.stdout
    assert result.stderr == "", result.stderr
    return jsonl(reviews_path(repo))


# ---------------------------------------------------------------- test 33


@pytest.mark.parametrize(
    "who", [{"agent_type": "developer", "agent_id": "d1"}, {"agent_id": "x1"}, {}],
    ids=["developer", "no-agent-type", "main-session"],
)  # fmt: skip
def test_33_a_handback_from_anyone_but_reviewer_records_nothing(
    repo: Path, who: dict[str, str]
) -> None:
    assert recorded(repo, handback_event(repo, HANDBACK, **who)) == []


@pytest.mark.parametrize("mode", ["off", "on", "broken"])
def test_33_a_reviewer_handback_records_one_line_in_every_mode(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    [line] = recorded(repo, handback_event(repo, HANDBACK, **REVIEWER))
    assert set(line) == {"at", "sha", "verdict", "agent_id", "session_id", "event"}
    assert (line["sha"], line["verdict"]) == (SHA, "APPROVED")
    assert (line["agent_id"], line["session_id"], line["event"]) == ("r1", "s1", "PreToolUse")


def test_33_a_subagent_stop_from_reviewer_records_its_last_message(repo: Path) -> None:
    [line] = recorded(repo, stop_event(repo, HANDBACK, **REVIEWER))
    assert (line["sha"], line["verdict"], line["event"]) == (SHA, "APPROVED", "SubagentStop")


def test_33_a_subagent_stop_from_another_persona_records_nothing(repo: Path) -> None:
    assert recorded(repo, stop_event(repo, HANDBACK, agent_type="tester", agent_id="t1")) == []


@pytest.mark.parametrize("verdict", ["APPROVED WITH WORDINGS", "CHANGES REQUESTED"])
def test_each_verdict_is_recorded_as_written(repo: Path, verdict: str) -> None:
    [line] = recorded(repo, handback_event(repo, f"Result.\nReview: {verdict} {SHA}", **REVIEWER))
    assert (line["sha"], line["verdict"]) == (SHA, verdict)


def test_records_append_newest_last(repo: Path) -> None:
    recorded(repo, handback_event(repo, f"Review: CHANGES REQUESTED {OTHER}", **REVIEWER))
    lines = recorded(repo, handback_event(repo, HANDBACK, **REVIEWER))
    assert [(r["sha"], r["verdict"]) for r in lines] == [
        (OTHER, "CHANGES REQUESTED"),
        (SHA, "APPROVED"),
    ]


# ---------------------------------------------------------------- test 34


@pytest.mark.parametrize(
    "line",
    [
        f"Review: APPROVED {SHA[:39]}",
        f"Review: APPROVED {SHA.upper()}",
        "Review: APPROVED",
        f"Review: LGTM {SHA}",
        f"Review: APPROVED {SHA} extra",
        f"  Review: APPROVED {SHA}",
    ],
    ids=["39-hex", "upper-case", "no-sha", "unknown-verdict", "trailing-word", "indented"],
)
def test_34_a_malformed_review_line_records_nothing(repo: Path, line: str) -> None:
    assert recorded(repo, handback_event(repo, f"Result: done.\n{line}", **REVIEWER)) == []


# ---------------------------------------------------------------- test 35


@pytest.mark.parametrize(
    "lines",
    [(f"Review: APPROVED {SHA}", f"Review: APPROVED {SHA}"),
     (f"Review: APPROVED {SHA}", f"Review: CHANGES REQUESTED {OTHER}"),
     (f"Review: CHANGES REQUESTED {OTHER}", f"Review: APPROVED {SHA}")],
    ids=["the-same-twice", "approved-first", "changes-first"],
)  # fmt: skip
def test_35_two_review_lines_record_nothing(repo: Path, lines: tuple[str, str]) -> None:
    message = f"Result.\n{lines[0]}\nmore text\n{lines[1]}"
    assert recorded(repo, handback_event(repo, message, **REVIEWER)) == []


# ---------------------------------------------------------------- test 36


@pytest.fixture
def unwritable_state(repo: Path) -> Iterator[Path]:
    if os.geteuid() == 0:  # pragma: no cover - root writes through any mode
        pytest.skip("root ignores directory permissions")
    state = state_dir(repo)
    state.mkdir(parents=True, exist_ok=True)
    state.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        yield state
    finally:
        state.chmod(stat.S_IRWXU)


def test_36_an_unwritable_state_directory_records_nothing_and_still_exits_0(
    repo: Path, unwritable_state: Path
) -> None:
    result = run_script(repo, RECORD_REVIEW, handback_event(repo, HANDBACK, **REVIEWER))
    assert (result.returncode, result.stdout, result.stderr) == (0, "", "")
    assert not reviews_path(repo).exists()
