"""`.claude/hooks/guard_unattended.py`: the three catch-all hooks (h3 U1, T11).

The spec is `docs/increments/h3-unattended-u1.md` §3.7 and §4. When the mode is
off the hook is silent for every event; when it is on or broken it queues the
act and refuses it in the output shape each event takes. The two non-PreToolUse
shapes were checked against the Claude Code `hooks` page (*PermissionRequest
decision control*: `hookSpecificOutput.decision.behavior` and `message`;
*ConfigChange decision control*: top-level `decision: "block"` and `reason`).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import (
    GUARD_UNATTENDED,
    make_repo,
    queue_lines,
    run_script,
    set_mode,
)

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
