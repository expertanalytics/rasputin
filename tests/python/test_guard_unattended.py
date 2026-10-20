"""`.claude/hooks/guard_unattended.py`: the five hooks (h3 U1, T11; h26; h28).

The spec is `docs/increments/h3-unattended-u1.md` §3.7 and §4. When the mode is
off the hook is silent for every event; when it is on or broken it queues the
act and refuses it in the output shape each event takes. The two non-PreToolUse
shapes were checked against the Claude Code `hooks` page (*PermissionRequest
decision control*: `hookSpecificOutput.decision.behavior` and `message`;
*ConfigChange decision control*: top-level `decision: "block"` and `reason`).

The fourth hook, `Stop`, is `docs/increments/h26-away-stop.md` §2 and §4: in
mode `on`, a main-session stop with nothing running is sent back once with a
reason, never queued; every other stop goes through.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import (
    GUARD_UNATTENDED,
    SUBAGENT,
    flag_granting,
    flag_on,
    jsonl,
    make_repo,
    queue_lines,
    run_script,
    set_mode,
    state_dir,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

QUESTION = "Which framework?"
COMMAND = "rm -rf node_modules"


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


def ask_event(repo: Path, question: str = QUESTION) -> dict[str, Any]:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": "AskUserQuestion",
        "tool_input": {
            "questions": [
                {"question": question, "header": "Q", "options": [{"label": "A"}]},
                {"question": "A second one?", "header": "Q2", "options": [{"label": "B"}]},
            ]
        },
        "cwd": str(repo),
    }


def permission_event(repo: Path) -> dict[str, Any]:
    return {
        "hook_event_name": "PermissionRequest",
        "tool_name": "Bash",
        "tool_input": {"command": COMMAND, "description": "Remove node_modules"},
        "cwd": str(repo),
    }


def config_event(repo: Path) -> dict[str, Any]:
    return {
        "hook_event_name": "ConfigChange",
        "source": "project_settings",
        "file_path": str(repo / ".claude" / "settings.json"),
        "cwd": str(repo),
    }


def output(repo: Path, event: dict[str, Any]) -> dict[str, Any] | None:
    result = run_script(repo, GUARD_UNATTENDED, event)
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    return json.loads(result.stdout) if result.stdout.strip() else None


def refusal_text(mode: str, text: str) -> None:
    """The queue reason (on) or the broken reason; both keep the ASK OLA instruction."""
    assert text.startswith("Refused: ")
    assert "ASK OLA:" in text
    if mode == "on":
        assert "unattended mode is on until" in text
    else:
        assert "cannot be read" in text


@pytest.mark.parametrize(
    "event", [ask_event, permission_event, config_event], ids=["ask", "permission", "config"]
)
def test_the_hook_is_silent_while_attended(repo: Path, event: Any) -> None:
    assert output(repo, event(repo)) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["on", "broken"])
def test_ask_user_question_is_denied_and_queued(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    found = output(repo, ask_event(repo))
    assert found is not None
    specific = found["hookSpecificOutput"]
    assert specific["hookEventName"] == "PreToolUse"
    assert specific["permissionDecision"] == "deny"
    refusal_text(mode, specific["permissionDecisionReason"])
    assert (
        "a question needs Ola; write it as the ASK OLA line"
        in (specific["permissionDecisionReason"])
    )
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_unattended", f"AskUserQuestion: {QUESTION}")


def test_ask_user_question_queues_only_the_first_200_characters(repo: Path) -> None:
    set_mode(repo, "on")
    question = "Q" * 500
    assert output(repo, ask_event(repo, question)) is not None
    [line] = queue_lines(repo)
    assert line["act"] == f"AskUserQuestion: {question[:200]}"


@pytest.mark.parametrize("mode", ["on", "broken"])
def test_a_permission_request_is_denied_with_a_message_and_queued(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    found = output(repo, permission_event(repo))
    assert found is not None
    assert set(found) == {"hookSpecificOutput"}
    specific = found["hookSpecificOutput"]
    assert specific["hookEventName"] == "PermissionRequest"
    assert specific["decision"]["behavior"] == "deny"
    refusal_text(mode, specific["decision"]["message"])
    assert "a permission prompt needs Ola" in specific["decision"]["message"]
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_unattended", f"Bash: {COMMAND}")


@pytest.mark.parametrize("mode", ["on", "broken"])
def test_a_config_change_is_blocked_and_queued(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    found = output(repo, config_event(repo))
    assert found == {
        "decision": "block",
        "reason": "unattended: configuration changes wait for Ola",
    }
    [line] = queue_lines(repo)
    settings = repo / ".claude" / "settings.json"
    assert (line["hook"], line["act"]) == (
        "guard_unattended",
        f"ConfigChange project_settings {settings}",
    )


# ------------------------------------------------ T19, T20: failure direction
#
# §3.9 (review 1): once stdin parses as a JSON object, any failure refuses in
# the shape of the event, in both modes, without `harness_mode` and without a
# queue line. A silent exit 1 would let the prompt through to nobody.

FAILED_PREFIX = "guard_unattended failed ("
FAILED_SUFFIX = (
    "), so this is refused rather than left for a prompt nobody may answer. It is not "
    "queued: record it as an ASK OLA line (main session: in session.md; subagent: in its "
    "handback) and continue."
)


def failure_reason(event_name: str, found: dict[str, Any] | None) -> str:
    """R from the §3.9 failure shape for `event_name`; fails on any other shape."""
    assert found is not None, "the hook printed nothing"
    if event_name == "ConfigChange":
        assert set(found) == {"decision", "reason"}
        assert found["decision"] == "block"
        reason: str = found["reason"]
    elif event_name == "PermissionRequest":
        assert set(found) == {"hookSpecificOutput"}
        specific = found["hookSpecificOutput"]
        assert specific["hookEventName"] == "PermissionRequest"
        assert specific["decision"]["behavior"] == "deny"
        reason = specific["decision"]["message"]
    else:
        assert set(found) == {"hookSpecificOutput"}
        specific = found["hookSpecificOutput"]
        assert specific["hookEventName"] == "PreToolUse"
        assert specific["permissionDecision"] == "deny"
        reason = specific["permissionDecisionReason"]
    assert reason.startswith(FAILED_PREFIX), reason
    assert reason.endswith(FAILED_SUFFIX), reason
    return reason


def with_list_input(event: dict[str, Any]) -> dict[str, Any]:
    return {**event, "tool_input": ["not", "an", "object"]}


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("event", [ask_event, permission_event], ids=["ask", "permission"])
def test_a_crash_after_parsing_refuses_in_the_events_shape_and_queues_nothing(
    repo: Path, mode: str, event: Any
) -> None:
    set_mode(repo, mode)
    malformed = with_list_input(event(repo))
    failure_reason(malformed["hook_event_name"], output(repo, malformed))
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("stdin", ["[1, 2]", "{not json"], ids=["json-list", "not-json"])
def test_stdin_that_is_not_a_json_object_prints_nothing(repo: Path, mode: str, stdin: str) -> None:
    set_mode(repo, mode)
    result = run_script(repo, GUARD_UNATTENDED, stdin)
    assert (result.returncode, result.stdout) == (0, ""), result.stderr
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize(
    "event", [ask_event, permission_event, config_event], ids=["ask", "permission", "config"]
)
def test_a_missing_harness_mode_refuses_in_the_events_shape(
    repo: Path, mode: str, event: Any
) -> None:
    set_mode(repo, mode)
    (repo / "tools" / "harness_mode.py").unlink()
    supplied = event(repo)
    reason = failure_reason(supplied["hook_event_name"], output(repo, supplied))
    assert reason.startswith(f"{FAILED_PREFIX}ModuleNotFoundError"), reason
    assert queue_lines(repo) == []


# ------------------------------------------------ h26: the Stop hook
#
# docs/increments/h26-away-stop.md §4. The idle event is the design's; each
# pass case changes one field of it. A stop is never queued.


#: A `stop_event` value that removes the key instead of setting it.
DROP = object()


def stop_event(repo: Path, **changed: Any) -> dict[str, Any]:
    """The design's idle Stop event, with `changed` applied; a value of DROP removes the key."""
    event: dict[str, Any] = {
        "hook_event_name": "Stop",
        "cwd": str(repo),
        "stop_hook_active": False,
        "last_assistant_message": "Done.",
        "background_tasks": [],
        "session_crons": [],
    }
    for key, value in changed.items():
        if value is DROP:
            event.pop(key)
        else:
            event[key] = value
    return event


def task(kind: str) -> dict[str, str]:
    return {"id": "t1", "type": kind, "status": "running", "description": "x"}


def test_an_idle_stop_while_away_is_blocked_once_with_the_reason(repo: Path) -> None:
    until = flag_on(repo)
    found = output(repo, stop_event(repo))
    assert found is not None, "the hook let an idle stop through in mode on"
    assert set(found) == {"decision", "reason"}
    assert found["decision"] == "block"
    reason = found["reason"]
    assert f"Unattended mode is on until {until:%Y-%m-%d %H:%M} UTC" in reason, reason
    for fragment in ("nothing is running", "ASK OLA:", "QUEUE:", "NOW: idle:"):
        assert fragment in reason, f"{fragment!r} missing from {reason!r}"
    assert queue_lines(repo) == []


PASSING_STOPS = {
    "already-sent-back": {"stop_hook_active": True},
    "stop-hook-active-absent": {"stop_hook_active": DROP},
    "subagent": {"agent_id": "a1"},
    "task-registry-unreachable": {"background_tasks": DROP},
    "shell-running": {"background_tasks": [task("shell")]},
    "subagent-running": {"background_tasks": [task("subagent")]},
    "monitor-running": {"background_tasks": [task("monitor")]},
    "cron-scheduled": {
        "session_crons": [{"id": "c1", "schedule": "*/5 * * * *", "prompt": "check"}]
    },
}


@pytest.mark.parametrize("changed", PASSING_STOPS.values(), ids=PASSING_STOPS.keys())
def test_a_stop_that_is_not_idle_or_already_sent_back_goes_through(
    repo: Path, changed: dict[str, Any]
) -> None:
    flag_on(repo)
    assert output(repo, stop_event(repo, **changed)) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "expired", "broken"])
def test_an_idle_stop_goes_through_unless_the_mode_is_on(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    assert output(repo, stop_event(repo)) is None
    assert queue_lines(repo) == []


def test_an_idle_stop_goes_through_when_harness_mode_is_missing(repo: Path) -> None:
    flag_on(repo)
    (repo / "tools" / "harness_mode.py").unlink()
    assert output(repo, stop_event(repo)) is None
    assert queue_lines(repo) == []


# ------------------------------------------------ h28: the prompt record and the window note
#
# docs/increments/h28-away-prompt.md §2 and §4, tests 1 to 4 of §6. While Ola
# is present, each permission prompt the hook lets open is written to
# `prompts.jsonl`; while away, nothing is (the prompt is refused). In mode on,
# the main session's first tool call per window gets a note naming the end of
# the window and the granted branches.

PROMPT_INPUT = {"command": "git push", "description": "Push the branch"}


def prompts_lines(repo: Path) -> list[dict[str, Any]]:
    return jsonl(state_dir(repo) / "prompts.jsonl")


def recorded_permission_event(repo: Path, **who: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "PermissionRequest",
        "tool_name": "Bash",
        "tool_input": PROMPT_INPUT,
        "session_id": "5e55102a-0000-4000-8000-000000000001",
        "transcript_path": str(repo / "transcripts" / "5e55102a.jsonl"),
        "cwd": str(repo),
        **who,
    }


def input_sha(supplied: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(supplied, sort_keys=True).encode()).hexdigest()


@pytest.mark.parametrize("who", [{}, SUBAGENT], ids=["main-session", "subagent"])
def test_an_attended_permission_prompt_is_recorded_once(repo: Path, who: dict[str, Any]) -> None:
    event = recorded_permission_event(repo, **who)
    assert output(repo, event) is None
    [line] = prompts_lines(repo)
    assert line["session_id"] == event["session_id"]
    assert line["transcript_path"] == event["transcript_path"]
    assert line["agent_id"] == who.get("agent_id")
    assert line["agent_type"] == who.get("agent_type")
    assert line["tool_name"] == "Bash"
    assert line["input_sha"] == input_sha(PROMPT_INPUT)
    assert line["act"] == "Bash: git push"
    at = datetime.fromisoformat(line["at"])
    assert at.tzinfo is not None, line["at"]
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["on", "broken"])
def test_a_refused_permission_prompt_is_not_recorded(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    found = output(repo, recorded_permission_event(repo))
    assert found is not None
    assert found["hookSpecificOutput"]["decision"]["behavior"] == "deny"
    assert not (state_dir(repo) / "prompts.jsonl").exists()
    [queued] = queue_lines(repo)
    assert queued["act"] == "Bash: git push"


GRANTED = ["grant-alpha", "grant-beta"]


def post_tool_event(repo: Path, session_id: str = "main-1", **extra: Any) -> dict[str, Any]:
    return {
        "hook_event_name": "PostToolUse",
        "tool_name": "Bash",
        "tool_input": {"command": "git status"},
        "tool_response": {"stdout": "", "stderr": "", "interrupted": False},
        "session_id": session_id,
        "transcript_path": str(repo / "transcripts" / f"{session_id}.jsonl"),
        "cwd": str(repo),
        **extra,
    }


def note_of(found: dict[str, Any] | None) -> str:
    assert found is not None, "the hook gave the main session no note"
    assert set(found) == {"hookSpecificOutput"}
    specific = found["hookSpecificOutput"]
    assert specific["hookEventName"] == "PostToolUse"
    note: str = specific["additionalContext"]
    return note


def assert_window_note(note: str, until: datetime) -> None:
    assert f"Unattended mode is on until {until:%Y-%m-%d %H:%M} UTC" in note, note
    for name in GRANTED:
        assert name in note, f"{name!r} missing from {note!r}"
    assert "ASK OLA:" in note, note


def test_the_main_session_hears_the_window_once_per_session_and_window(repo: Path) -> None:
    now = datetime.now(UTC).replace(microsecond=0)
    since = flag_granting(repo, GRANTED, now)
    until = since + timedelta(hours=1)
    assert_window_note(note_of(output(repo, post_tool_event(repo))), until)
    assert output(repo, post_tool_event(repo)) is None, "the same session heard it twice"
    assert_window_note(note_of(output(repo, post_tool_event(repo, "main-2"))), until)
    assert output(repo, post_tool_event(repo, "main-3", **SUBAGENT)) is None
    # A new window (a new since) is heard again by a session that heard the last one.
    later = flag_granting(repo, GRANTED, now - timedelta(minutes=1))
    assert_window_note(note_of(output(repo, post_tool_event(repo))), later + timedelta(hours=1))
    assert queue_lines(repo) == []


def test_no_grants_are_named_as_none(repo: Path) -> None:
    flag_granting(repo, [])
    note = note_of(output(repo, post_tool_event(repo)))
    assert "Granted branches: none" in note, note


@pytest.mark.parametrize("mode", ["off", "expired"])
def test_no_window_note_unless_the_mode_is_on(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    assert output(repo, post_tool_event(repo)) is None
    assert not (state_dir(repo) / "noticed.jsonl").exists()


def test_a_post_tool_use_without_harness_mode_prints_nothing(repo: Path) -> None:
    flag_granting(repo, GRANTED)
    (repo / "tools" / "harness_mode.py").unlink()
    result = run_script(repo, GUARD_UNATTENDED, post_tool_event(repo))
    assert (result.returncode, result.stdout) == (0, ""), result.stderr
