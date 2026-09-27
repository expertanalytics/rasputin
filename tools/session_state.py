#!/usr/bin/env python3
"""Print the round recap, then what the previous session was in the middle of.

The recap (retrospective rule 5, Ola): the last thing landed, what is in
flight, decisions waiting on Ola, and the next three ROADMAP items -- all read
from files on disk and git, so it is the same whoever runs it.

Then it answers one question: when a session starts cold or resumes after a
context loss, what was actually being asked? Reads `.claude/current-task/` (the in-flight
asks -- `session.md` first, each subagent's file after it as context) and the
predecessor session transcript under ~/.claude/projects/, including prompts the
user queued and the harness absorbed mid-turn -- which is where a lost
instruction hides, because such a prompt never appears as a normal user turn.

Known gap: only user entries whose content is a plain string are captured, so a
prompt carrying an attachment or image arrives as a list and is dropped. No turn
in this project is currently lost that way, but the failure is silent, which is
the worst mode for a tool whose job is to surface a dropped turn.

Partly tested: tests/python/test_session_state.py covers the recap functions
only. The transcript and current-task readers are still untested, and every
defect in them so far was found by hand. Cover the input classes that actually
bit, since two of them nobody guessed: FIFO (open() blocks with no writer, and
an in-process timeout cannot see it), symlink loop, dangling symlink,
directory, non-UTF-8, chmod 000, no extension, missing directory.

Usage: python3 tools/session_state.py [--turns N]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import stat
import subprocess
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
    """Contents, or a named placeholder.

    A recovery tool that dies on one unreadable file takes the human-turn
    half -- the reason it exists -- down with it. Measured twice: a chmod 000
    file raised PermissionError, and once the glob widened past *.md, a
    non-UTF-8 file raised UnicodeDecodeError, which is a ValueError and so
    passed both excepts. Hence errors="replace", as human_turns already does.
    """
    try:
        # Only a regular file may be opened. Dropping is_file() to surface broken
        # symlinks also admitted FIFOs, and open() on a FIFO with no writer BLOCKS
        # -- a recovery tool that hangs is worse than one that crashes, because
        # nothing prints at all. Measured before this guard: with a FIFO in the
        # directory the run timed out; without it, it completed. Time it from
        # outside the process -- signal.alarm raises TimeoutError, which is an
        # OSError, which the except below would swallow.
        if not stat.S_ISREG(path.stat().st_mode):
            return "(not a regular file)"
        return path.read_text(errors="replace").rstrip()
    except FileNotFoundError:
        return absent
    except OSError:
        return "(unreadable)"


def _git(repo: Path, *args: str) -> str | None:
    """stdout of a git command, or None if it failed."""
    result = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def last_landed(repo: Path) -> str:
    """The newest merge commit reachable from HEAD: how work lands here."""
    found = _git(repo, "log", "-1", "--merges", "--format=%h %cs %s")
    return found or "(no merge commit reachable from HEAD)"


def default_base(repo: Path) -> str:
    """origin/master when it resolves, else master.

    A local master lags the remote after a PR merges on GitHub, and counting
    against it lists already-merged commits as in flight.
    """
    remote = _git(repo, "rev-parse", "--verify", "--quiet", "origin/master")
    return "origin/master" if remote else "master"


def in_flight(repo: Path, base: str = "master") -> list[str]:
    """The current branch and its commits ahead of base, newest first."""
    branch = _git(repo, "rev-parse", "--abbrev-ref", "HEAD") or "(unknown)"
    ahead = _git(repo, "log", "--format=%h %s", f"{base}..HEAD")
    if ahead is None:
        return [f"branch {branch} (base {base} not found)"]
    commits = ahead.splitlines()
    return [f"branch {branch}, {len(commits)} commit(s) ahead of {base}", *commits]


def pending_decisions(tasks: Path) -> list[str]:
    """Lines marked `ASK OLA` in any current-task file, any case."""
    if not tasks.is_dir():
        return []
    found: list[str] = []
    for path in sorted(p for p in tasks.glob("*") if not p.is_dir()):
        for line in _read(path, "").splitlines():
            if "ask ola" in line.lower():
                found.append(f"{path.name}: {line.strip()}")
    return found


# A table row: | # | What it is | Status | Record |
_ROW = re.compile(r"^\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|")
_DONE = ("shipped", "landed")


def roadmap_next(text: str, n: int = 3) -> list[str]:
    """The first n open rows of ROADMAP.md's first table, with their status.

    Open means the status neither starts with "shipped" or "landed" nor says
    "unscheduled". Table order is the reading order; the status is printed so
    the reader, not this function, judges what is next.
    """
    rows: list[str] = []
    in_table = False
    for line in text.splitlines():
        match = _ROW.match(line)
        if not match:
            if in_table:
                break
            continue
        in_table = True
        number, what, status = match.groups()
        if number == "#" or set(number) <= {"-"}:
            continue
        plain = status.replace("*", "")
        if plain.lower().startswith(_DONE) or "unscheduled" in plain.lower():
            continue
        rows.append(f"{number}: {what[:100]} [{plain}]")
    return rows[:n]


def print_recap() -> None:
    """Retrospective rule 5: the structured recap every round opens with."""
    print("== recap ==")
    print(f"Last landed: {last_landed(REPO)}")
    print("In flight:")
    flight = in_flight(REPO, default_base(REPO))
    for line in flight[:11]:
        print(f"  {line}")
    if len(flight) > 11:
        print(f"  ... {len(flight) - 11} older")
    decisions = pending_decisions(REPO / ".claude" / "current-task")
    print("Waiting on Ola:" + ("" if decisions else " (none recorded as ASK OLA)"))
    for line in decisions:
        print(f"  {line}")
    print("Next on ROADMAP.md:")
    for row in roadmap_next(_read(REPO / "ROADMAP.md", "")):
        print(f"  {row}")
    print()


def print_current_task() -> None:
    """Print the session's ask first, then each subagent's as context.

    A cold session needs one of these promoted, not a flat dump: `session.md` is
    the record of what the round is for, and a subagent file is a fragment of it
    delegated. Both are printed because a subagent file may be the only trace of
    a step that died, but the order says which one to believe about the round.

    Sorted by name for stable output. session.md is promoted above the rest
    explicitly, so the order among subagent files decides nothing.
    """
    tasks = REPO / ".claude" / "current-task"
    session = tasks / "session.md"
    print("== .claude/current-task/session.md ==")
    print(_read(session, "(absent -- no session-level ask was recorded in flight)"))

    # Every entry, not just *.md, and not p.is_file(): a subagent that names its
    # file without an extension, or leaves a broken symlink, must not become
    # invisible to the tool whose job is surfacing what would otherwise be lost.
    # is_file() is false for a dangling symlink, so it dropped one silently;
    # not is_dir() surfaces it and _read's FileNotFoundError names it.
    others = sorted(p for p in tasks.glob("*") if not p.is_dir() and p.name != "session.md")
    if not others:
        return
    # No staleness flag here. It was computed from session.md's mtime, but
    # session.md is overwritten as a round progresses, so an ordinary
    # write-spawn-update sequence made every LIVE subagent file predate it and
    # get marked for sweeping -- the one file a cold session must not lose.
    # Liveness is not visible from a directory listing, so the lifecycle rule
    # in .claude/REQUIRED-READING.md owns it and this prints what is there.
    print(f"\n== {len(others)} subagent ask(s), as context ==")
    for path in others:
        try:
            stamp = datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        except OSError:
            stamp = "unknown time"
        print(f"\n-- {path.name} ({stamp})")
        print(_read(path, "(unreadable)"))
    print(
        "\nThe spawner deletes a subagent file once its handback is read; a file"
        " whose writer died is deleted once the turns below have been read."
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--turns", type=int, default=5)
    args = parser.parse_args()

    print_recap()
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
