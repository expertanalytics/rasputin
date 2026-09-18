#!/usr/bin/env python3
"""Print what the previous session was in the middle of.

Answers one question: when a session starts cold or resumes after a context
loss, what was actually being asked? Reads `.claude/current-task/` (the in-flight
asks -- `session.md` first, each subagent's file after it as context) and the
predecessor session transcript under ~/.claude/projects/, including prompts the
user queued and the harness absorbed mid-turn -- which is where a lost
instruction hides, because such a prompt never appears as a normal user turn.

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
from datetime import datetime
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


def _read(path: Path, absent: str) -> str:
    """Contents, or a named placeholder. Never raises.

    A recovery tool that dies on one unreadable file takes the human-turn
    half -- the reason it exists -- down with it. Measured: a chmod 000 file
    raised PermissionError out of read_text() and main() never reached the
    transcript scan.
    """
    try:
        return path.read_text().rstrip()
    except FileNotFoundError:
        return absent
    except OSError:
        return "(unreadable)"


def print_current_task() -> None:
    """Print the session's ask first, then each subagent's as context.

    A cold session needs one of these promoted, not a flat dump: `session.md` is
    the record of what the round is for, and a subagent file is a fragment of it
    delegated. Both are printed because a subagent file may be the only trace of
    a step that died, but the order says which one to believe about the round.

    Ordering by mtime would be wrong here -- the newest file is whichever
    subagent wrote last, which is precisely not the thing to read first.
    """
    tasks = REPO / ".claude" / "current-task"
    session = tasks / "session.md"
    print("== .claude/current-task/session.md ==")
    print(_read(session, "(absent -- no session-level ask was recorded in flight)"))

    # Every entry, not just *.md: a subagent that names its file without an
    # extension must not become invisible to the tool whose job is surfacing
    # what would otherwise be lost.
    others = sorted(p for p in tasks.glob("*") if p.is_file() and p.name != "session.md")
    if not others:
        return
    # No staleness flag here. It was computed from session.md's mtime, but
    # session.md is overwritten as a round progresses, so an ordinary
    # write-spawn-update sequence made every LIVE subagent file predate it and
    # get marked for sweeping -- the one file a cold session must not lose.
    # Liveness is not visible from a directory listing, so the sweep rule in
    # .claude/REQUIRED-READING.md owns it and this prints what is there.
    print(f"\n== {len(others)} subagent ask(s), as context ==")
    for path in others:
        try:
            stamp = datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        except OSError:
            stamp = "unknown time"
        print(f"\n-- {path.name} ({stamp})")
        print(_read(path, "(unreadable)"))
    print(
        "\nDelete a subagent file once its handback is read; sweep any belonging"
        " to an agent you are no longer waiting on, after reading the turns below."
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--turns", type=int, default=5)
    args = parser.parse_args()

    print_current_task()

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
