#!/usr/bin/env python3
"""Print what the previous session was in the middle of.

Answers one question: when a session starts cold or resumes after a context
loss, what was actually being asked? Reads `.claude/current-task.md` (the
in-flight ask, if the previous session wrote one) and the predecessor session
transcript under ~/.claude/projects/, including prompts the user queued and the
harness absorbed mid-turn -- which is where a lost instruction hides, because
such a prompt never appears as a normal user turn.

Known gap: only user entries whose content is a plain string are captured, so a
prompt carrying an attachment or image arrives as a list and is dropped. No turn
in this project is currently lost that way, but the failure is silent, which is
the worst mode for a tool whose job is to surface a dropped turn.

Usage: python tools/session_state.py [--turns N]
"""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SLUG = re.sub(r"[^A-Za-z0-9]", "-", str(REPO))
TRANSCRIPTS = Path.home() / ".claude" / "projects" / SLUG


SKIP = (
    "<task-notification>",
    "<local-command",
    "<command-name>",
    "Caveat:",
    "Another Claude session sent a message:",
)


def human_turns(path: Path) -> list[tuple[str, str, str]]:
    """(timestamp, session, text) for human turns and absorbed queued prompts."""
    found: list[tuple[str, str, str]] = []
    for line in path.read_text(errors="replace").splitlines():
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue
        kind = entry.get("type")
        prefix = ""
        if kind == "user":
            content = entry.get("message", {}).get("content")
            text = content if isinstance(content, str) else None
        elif kind == "queue-operation" and entry.get("reason") == "absorbed_mid_turn":
            text = entry.get("content")
            prefix = "[queued, absorbed mid-turn] "
        else:
            continue
        if text and not text.lstrip().startswith(SKIP):
            text = prefix + text
            found.append((entry.get("timestamp", ""), path.stem[:8], " ".join(text.split())))
    return found


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--turns", type=int, default=5)
    args = parser.parse_args()

    current = REPO / ".claude" / "current-task.md"
    print("== .claude/current-task.md ==")
    if current.exists():
        print(current.read_text().rstrip())
    else:
        print("(absent -- no ask was recorded in flight)")

    # Claude Code exports CLAUDE_CODE_SESSION_ID; the older spelling is kept as a
    # fallback so the script still excludes the current session if that changes.
    here = os.environ.get("CLAUDE_CODE_SESSION_ID") or os.environ.get(
        "CLAUDE_SESSION_ID", ""
    )
    others = sorted(
        (p for p in TRANSCRIPTS.glob("*.jsonl") if p.stem != here),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not others:
        print(f"\n(no predecessor transcript under {TRANSCRIPTS})")
        return 0

    # Enough transcripts that --turns can always be satisfied, since a cold
    # session's own predecessor may itself be short.
    turns = [t for path in others[: max(3, args.turns)] for t in human_turns(path)]
    turns.sort(key=lambda t: t[0])
    print(f"\n== last {args.turns} human turns before this session ==")
    for stamp, session, text in turns[-args.turns :]:
        shown = text if len(text) <= 600 else text[:600] + " [...truncated]"
        print(f"- [{stamp[:19]} {session}] {shown}")
    print("\nA turn above with no answering commit or file is still pending. Ask before acting.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
