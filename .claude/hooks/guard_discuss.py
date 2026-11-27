#!/usr/bin/env python3
"""PreToolUse on AskUserQuestion: no new box while Ola formulates his point.

Ola, 2026-11-27: when he answers a box with anything other than one of its
options, no new box is asked until his next chat message. The latest answered
AskUserQuestion in the transcript opened a discussion if its result is an error
(he chose "Chat about this", or rejected or interrupted the box), or an answer
that is not one of the question's labels (typed text) or is a "Discuss" label.
Then the box is refused (exit 2) until a user entry from the person follows:
`origin.kind` human, or, without an origin, plain text that is not harness
content. Hand-backs, reminders, hook output and tool results do not count. A
missing or unreadable transcript passes with a note on stderr; a torn line is
skipped. Ola, 2026-11-27: no box outside question mode; while it is off
(`questions.py mode`, read from the git common dir of the event's `cwd`) every
box is refused. Tests: tests/python/test_guard_discuss.py. Standard library only.
"""

import json
import subprocess
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2] / "tools"))
import questions  # in tools/, standard library only

REASON = (
    "guard_discuss: Ola opened a discussion on the last question box (he chose to chat, "
    "discuss or typed his own answer); ask no new box until his next chat message."
)
OFF = "enter question mode first: python3 tools/questions.py mode enter"
HARNESS = ("<", "[Request interrupted", "Another Claude session", "Stop hook feedback")


def entries(path: Path) -> list[dict]:
    found = []
    for line in path.read_text().splitlines():
        try:
            entry = json.loads(line)
        except ValueError:
            continue
        if isinstance(entry, dict) and not entry.get("isSidechain"):
            found.append(entry)
    return found


def blocks(entry: dict) -> list[dict]:
    content = (entry.get("message") or {}).get("content")
    if isinstance(content, str):
        return [{"type": "text", "text": content}]
    return [b for b in content if isinstance(b, dict)] if isinstance(content, list) else []


def from_person(entry: dict) -> bool:
    """A user entry the person typed, not a tool result or harness content."""
    if entry.get("type") != "user" or entry.get("isMeta"):
        return False
    if "origin" in entry:
        return (entry.get("origin") or {}).get("kind") == "human"
    texts = [b.get("text", "") for b in blocks(entry) if b.get("type") == "text"]
    if any(b.get("type") == "tool_result" for b in blocks(entry)):
        return False
    return any(t.strip() and not t.lstrip().startswith(HARNESS) for t in texts)


def opened(result: dict, entry: dict) -> bool:
    """The box's result left the question open: an error, or an answer off its labels."""
    if result.get("is_error"):
        return True
    record = entry.get("toolUseResult")
    if not isinstance(record, dict) or not isinstance(record.get("answers"), dict):
        return True
    for question in record.get("questions") or []:
        labels = {o.get("label") for o in question.get("options") or []}
        answer = record["answers"].get(question.get("question"))
        if not isinstance(answer, str):
            return True
        parts = [answer] if answer in labels else answer.split(", ")  # multiSelect joins
        if any(p not in labels or p.lower().startswith("discuss") for p in parts):
            return True
    return False


def refuses(log: list[dict]) -> bool:
    asked = {
        b.get("id")
        for e in log
        if e.get("type") == "assistant"
        for b in blocks(e)
        if b.get("type") == "tool_use" and b.get("name") == "AskUserQuestion"
    }
    for at in range(len(log) - 1, -1, -1):
        entry = log[at]
        if entry.get("type") != "user":
            continue
        for block in blocks(entry):
            if block.get("type") == "tool_result" and block.get("tool_use_id") in asked:
                return opened(block, entry) and not any(map(from_person, log[at + 1 :]))
    return False


def mode_on(event: dict) -> bool:
    """Question mode is on in the repository of the event's `cwd` (else this hook's)."""
    cwd = event.get("cwd") or Path(__file__).resolve().parents[2]
    common = subprocess.run(
        ["git", "-C", str(cwd), "rev-parse", "--path-format=absolute", "--git-common-dir"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    return bool(common) and questions.mode_since(Path(common)) is not None


def main() -> int:
    try:
        event = json.load(sys.stdin)
        if not mode_on(event):
            print(OFF, file=sys.stderr)
            return 2
        path = Path(event["transcript_path"])
        log = entries(path)
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"guard_discuss: transcript not read ({error}); the box passes", file=sys.stderr)
        return 0
    if refuses(log):
        print(REASON, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
