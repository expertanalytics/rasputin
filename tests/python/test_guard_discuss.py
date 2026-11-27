"""`.claude/hooks/guard_discuss.py`: no new question box while Ola formulates his point.

Ola, 2026-11-27: when he answers a box with anything other than one of its
options (chat about this, discuss, typed text, an interrupt), no new box may be
asked until his next chat message. The fixtures are short synthetic copies of
the transcript shapes Claude Code writes: an assistant `tool_use`, a user
`tool_result` (with `toolUseResult`), a human message (`origin.kind` human),
agent hand-backs (`isMeta`, `origin.kind` peer) and the interrupt marker.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import REAL, clean_env, question_mode

pytestmark = pytest.mark.harness

HOOK = REAL / ".claude" / "hooks" / "guard_discuss.py"
OPTIONS = [{"label": "Yes"}, {"label": "No"}, {"label": "Discuss first"}]
QUESTION = "Merge the fix?"

CLARIFY = (
    "The user doesn't want to proceed with this tool use. The tool use was rejected. "
    "To tell you how to proceed, the user said:\nThe user wants to clarify these questions."
)
OFF = "enter question mode first: python3 tools/questions.py mode enter"
STOPPED = (
    "The user doesn't want to proceed with this tool use. The tool use was rejected. "
    "STOP what you are doing and wait for the user to tell you how to proceed."
)


def ask(tool_id: str, question: str = QUESTION) -> dict[str, Any]:
    use = {
        "type": "tool_use",
        "id": tool_id,
        "name": "AskUserQuestion",
        "input": {"questions": [{"question": question, "header": "Q", "options": OPTIONS}]},
    }
    return {"type": "assistant", "message": {"role": "assistant", "content": [use]}}


def answered(tool_id: str, answer: str, question: str = QUESTION) -> dict[str, Any]:
    text = f'Your questions have been answered: "{question}"="{answer}". You can now continue.'
    block = {"type": "tool_result", "content": text, "tool_use_id": tool_id}
    result = {
        "questions": [{"question": question, "options": OPTIONS}],
        "answers": {question: answer},
        "annotations": {},
    }
    return {
        "type": "user",
        "message": {"role": "user", "content": [block]},
        "toolUseResult": result,
    }


def rejected(tool_id: str, text: str = CLARIFY) -> dict[str, Any]:
    block = {"type": "tool_result", "content": text, "is_error": True, "tool_use_id": tool_id}
    return {
        "type": "user",
        "message": {"role": "user", "content": [block]},
        "toolUseResult": "Error: " + text,
        "toolDenialKind": "user-rejected",
    }


def human(text: str) -> dict[str, Any]:
    return {
        "type": "user",
        "origin": {"kind": "human"},
        "message": {"role": "user", "content": text},
    }


def interrupted() -> dict[str, Any]:
    text = [{"type": "text", "text": "[Request interrupted by user for tool use]"}]
    return {"type": "user", "message": {"role": "user", "content": text}}


def handback() -> dict[str, Any]:
    text = 'Another Claude session sent a message:\n<agent-message from="a1">\nResult: done.'
    return {
        "type": "user",
        "isMeta": True,
        "origin": {"kind": "peer", "handback": True},
        "message": {"role": "user", "content": text},
    }


def reminder() -> dict[str, Any]:
    text = "<system-reminder>\nA note.\n</system-reminder>"
    return {"type": "user", "isMeta": True, "message": {"role": "user", "content": text}}


def tool_output(tool_id: str) -> dict[str, Any]:
    block = {"type": "tool_result", "content": "ok", "tool_use_id": tool_id}
    return {"type": "user", "message": {"role": "user", "content": [block]}}


def bash_use(tool_id: str) -> dict[str, Any]:
    use = {"type": "tool_use", "id": tool_id, "name": "Bash", "input": {"command": "ls"}}
    return {"type": "assistant", "message": {"role": "assistant", "content": [use]}}


def run(
    tmp_path: Path,
    entries: list[dict[str, Any]] | None,
    raw: str | None = None,
    path: Path | None = None,
    mode: bool = True,
) -> subprocess.CompletedProcess[str]:
    """The hook on a box asked from `tmp_path`, a git repository with question mode on
    (or off, `mode=False`)."""
    if not (tmp_path / ".git").exists():
        subprocess.run(["git", "init", "-q", str(tmp_path)], env=clean_env(), check=True)
    if mode:
        question_mode(tmp_path / ".git")
    transcript = path or tmp_path / "t.jsonl"
    if entries is not None:
        transcript.write_text("".join(json.dumps(e) + "\n" for e in entries) + (raw or ""))
    event = {
        "hook_event_name": "PreToolUse",
        "tool_name": "AskUserQuestion",
        "tool_input": {"questions": [{"question": "Next?", "options": OPTIONS}]},
        "transcript_path": str(transcript),
        "cwd": str(tmp_path),
    }
    return subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(event),
        capture_output=True,
        text=True,
        env=clean_env(),
        timeout=30,
        check=False,
    )


def refused(result: subprocess.CompletedProcess[str]) -> bool:
    if result.returncode == 2:
        assert result.stderr.strip(), "a refusal names its reason on stderr"
        assert len(result.stderr.strip().splitlines()) == 1
        assert "Ola" in result.stderr
        return True
    assert result.returncode == 0, result.stderr
    return False


# -- a box that opened a discussion, with nothing from Ola since ---------------


@pytest.mark.parametrize(
    "result",
    [
        pytest.param(rejected("t1"), id="chat-about-this"),
        pytest.param(rejected("t1", STOPPED), id="rejected-stop"),
        pytest.param(answered("t1", "Can we split this first?"), id="typed-text"),
        pytest.param(answered("t1", "Discuss first"), id="discuss-option"),
    ],
)
def test_refused_after_a_discussion_box(tmp_path: Path, result: dict[str, Any]) -> None:
    assert refused(run(tmp_path, [ask("t1"), result]))


def test_refused_after_an_interrupt(tmp_path: Path) -> None:
    assert refused(run(tmp_path, [ask("t1"), rejected("t1", STOPPED), interrupted()]))


def test_refused_when_only_harness_content_follows(tmp_path: Path) -> None:
    entries = [
        ask("t1"),
        rejected("t1"),
        handback(),
        reminder(),
        bash_use("b1"),
        tool_output("b1"),
        interrupted(),
    ]
    assert refused(run(tmp_path, entries))


def test_refused_with_the_new_box_already_in_the_transcript(tmp_path: Path) -> None:
    """The box being asked may already be written, with no result yet."""
    assert refused(run(tmp_path, [ask("t1"), answered("t1", "Discuss first"), ask("t2")]))


# -- passes ----------------------------------------------------------------------


def test_passes_after_a_listed_answer(tmp_path: Path) -> None:
    assert not refused(run(tmp_path, [ask("t1"), answered("t1", "Yes"), handback()]))


def test_passes_once_ola_has_written(tmp_path: Path) -> None:
    entries = [ask("t1"), rejected("t1"), handback(), human("Split it in two."), handback()]
    assert not refused(run(tmp_path, entries))


def test_passes_with_no_earlier_box(tmp_path: Path) -> None:
    assert not refused(run(tmp_path, [human("Good morning!"), bash_use("b1"), tool_output("b1")]))


def test_only_the_latest_box_counts(tmp_path: Path) -> None:
    entries = [ask("t1"), rejected("t1"), human("Go on."), ask("t2"), answered("t2", "No")]
    assert not refused(run(tmp_path, entries))


def test_a_later_discussion_box_refuses_again(tmp_path: Path) -> None:
    entries = [ask("t1"), answered("t1", "Yes"), ask("t2"), rejected("t2")]
    assert refused(run(tmp_path, entries))


def test_a_torn_last_line_is_skipped(tmp_path: Path) -> None:
    assert not refused(run(tmp_path, [ask("t1"), answered("t1", "Yes")], raw='{"type": "us'))


def test_missing_transcript_passes_with_a_note(tmp_path: Path) -> None:
    result = run(tmp_path, None, path=tmp_path / "absent.jsonl")
    assert result.returncode == 0
    assert "guard_discuss" in result.stderr


def test_no_transcript_path_passes(tmp_path: Path) -> None:
    subprocess.run(["git", "init", "-q", str(tmp_path)], env=clean_env(), check=True)
    question_mode(tmp_path / ".git")
    event = {"hook_event_name": "PreToolUse", "tool_name": "AskUserQuestion", "tool_input": {},
             "cwd": str(tmp_path)}  # fmt: skip
    result = subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(event),
        capture_output=True,
        text=True,
        env=clean_env(),
        timeout=30,
        check=False,
    )
    assert result.returncode == 0


# -- question mode -----------------------------------------------------------------


@pytest.mark.parametrize(
    "entries",
    [
        pytest.param([], id="no-earlier-box"),
        pytest.param([ask("t1"), answered("t1", "Yes"), human("Go on.")], id="after-a-yes"),
    ],
)
def test_refused_while_question_mode_is_off(tmp_path: Path, entries: list[dict[str, Any]]) -> None:
    result = run(tmp_path, entries, mode=False)
    assert result.returncode == 2
    assert result.stderr.strip() == OFF


def test_refused_outside_a_git_repository(tmp_path: Path) -> None:
    """No git common dir, so no mode state: the mode reads as off."""
    transcript = tmp_path / "t.jsonl"
    transcript.write_text("")
    event = {"tool_name": "AskUserQuestion", "transcript_path": str(transcript),
             "cwd": str(tmp_path / "absent")}  # fmt: skip
    result = subprocess.run([sys.executable, str(HOOK)], input=json.dumps(event),
                            capture_output=True, text=True, env=clean_env(), timeout=30,
                            check=False)  # fmt: skip
    assert result.returncode == 2
    assert result.stderr.strip() == OFF
