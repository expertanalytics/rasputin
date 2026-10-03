#!/usr/bin/env python3
"""Enter or leave unattended mode. Ola runs this, in his own terminal.

    python3 tools/away.py 8h        # or 90m, 1h30m; at most 72h
    python3 tools/away.py --back

Spec: docs/increments/h3-unattended-u1.md §3.2 and §3.3. Entering needs a `y`
typed on /dev/tty, which a pipe cannot supply; leaving needs no terminal,
because it loosens nothing. The guards deny agents both (§3.6).
"""

from __future__ import annotations

import argparse
import contextlib
import errno
import io
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
from collections.abc import Callable
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any, NamedTuple, TextIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness_mode as hm

DURATION = re.compile(r"^(?:(\d+)h)?(?:(\d+)m)?$")
USAGE = "usage: python3 tools/away.py <duration: 8h, 90m, 1h30m; at most 72h> | --back"


def stop_keep_awake(pid: object) -> None:
    """SIGTERM the keep-awake, only if `pid` is an int naming a running caffeinate."""
    if not isinstance(pid, int) or isinstance(pid, bool):
        return
    found = subprocess.run(
        ["ps", "-p", str(pid), "-o", "comm="], capture_output=True, text=True, check=False
    )
    if found.returncode == 0 and os.path.basename(found.stdout.strip()) == "caffeinate":
        with contextlib.suppress(OSError):
            os.kill(pid, signal.SIGTERM)


def minutes_of(duration: str | None) -> int | None:
    match = DURATION.match(duration or "")
    if match is None or not any(match.groups()):
        return None
    minutes = int(match[1] or 0) * 60 + int(match[2] or 0)
    return minutes if 0 < minutes <= hm.MAX_HOURS * 60 else None


def read_flag(state: Path) -> dict[str, Any] | None:
    """The flag as written, `{}` if present but unreadable, None if absent."""
    path = state / "unattended.json"
    if not os.path.lexists(path):
        return None
    try:
        flag = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    return flag if isinstance(flag, dict) else {}


def append(path: Path, entry: dict[str, Any]) -> None:
    with path.open("a") as out:
        out.write(json.dumps(entry) + "\n")


def when(moment: datetime) -> str:
    return f"{moment.astimezone():%Y-%m-%d %H:%M %Z} ({moment.isoformat(timespec='seconds')})"


class Quiet(NamedTuple):
    """The longest interval without a commit; `opener` is None at the window's start.

    A NamedTuple, not a frozen dataclass: both refuse assignment, but a dataclass
    with string annotations needs its module in sys.modules, and the pty tests
    load this file without registering it.
    """

    start: datetime
    end: datetime
    opener: str | None
    count: int  # type: ignore[assignment]  # shadows tuple.count; the name is pinned


def longest_quiet(since: datetime, end: datetime, commits: list[tuple[datetime, str]]) -> Quiet:
    """The longest gap between since, the commits inside [since, end], and end (h8 §3.4)."""
    inside = sorted(c for c in commits if since <= c[0] <= end)
    points: list[tuple[datetime, str | None]] = [(since, None), *inside, (end, None)]
    best = max(range(len(points) - 1), key=lambda i: (points[i + 1][0] - points[i][0], -i))
    return Quiet(points[best][0], points[best + 1][0], points[best][1], len(inside))


def commits_between(
    root: Path, since: datetime, end: datetime
) -> list[tuple[datetime, str]] | None:
    """(commit time, short hash) on any ref inside [since, end]; None if git failed.

    No --since/--until: git stops a date-limited walk early on out-of-order dates.
    """
    log = subprocess.run(
        ["git", "-C", str(root), "log", "--all", "--format=%h %cI"],
        capture_output=True,
        text=True,
        check=False,
    )
    if log.returncode:
        return None
    found = []
    for line in log.stdout.splitlines():
        short, _, stamp = line.partition(" ")
        try:  # an unparsable or naive date is a failed log: --back must still archive
            if since <= (moment := datetime.fromisoformat(stamp)) <= end:
                found.append((moment, short))
        except (ValueError, TypeError):
            return None
    return found


def _moment(value: object) -> datetime | None:
    """A timezone-aware time from the flag, or None."""
    try:
        moment = datetime.fromisoformat(str(value))
    except ValueError:
        return None
    return moment if moment.tzinfo else None


def stretch(root: Path, flag: dict[str, Any], now: datetime) -> str:
    """The line `--back` prints about the longest stretch without a commit."""
    unknown = "Longest stretch without a commit: unknown ({})."
    since = _moment(flag.get("since"))
    if since is None:
        return unknown.format("the flag has no readable since")
    end = min(now, _moment(flag.get("until")) or now)
    commits = commits_between(root, since, end)
    if commits is None:
        return unknown.format("git log failed")
    quiet = longest_quiet(since, end, commits)
    minutes = int((quiet.end - quiet.start).total_seconds() // 60)
    opener = f"after {quiet.opener}" if quiet.opener else "from the window's start"
    start, stop = (f"{t.astimezone(UTC):%Y-%m-%d %H:%M}" for t in (quiet.start, quiet.end))
    return (
        f"Longest stretch without a commit (any ref): {minutes // 60} h {minutes % 60:02d} min, "
        f"{start} to {stop} UTC, {opener}; {quiet.count} commit(s) in the window."
    )


def enter(
    duration: str,
    minutes: int,
    state: Path,
    now: datetime,
    open_tty: Callable[[], TextIO],
    spawn: Callable[..., Any],
    stop: Callable[[object], None],
) -> int:
    try:
        tty = open_tty()
    except OSError as error:
        # errno None: no OS error behind it, so name the exception, not a code.
        name = (
            f"{type(error).__name__}: {error}"
            if error.errno is None
            else errno.errorcode.get(error.errno, str(error))
        )
        print(
            f"away.py needs your own terminal: /dev/tty is not available ({name}). "
            "Run it in a terminal window, not through the agent.",
            file=sys.stderr,
        )
        return 3
    seconds = math.ceil(round(minutes * 60 * hm.BUFFER, 6))
    until = now + timedelta(seconds=seconds)
    old = read_flag(state)
    with tty:
        tty.write(
            f"Unattended mode for {duration} x {hm.BUFFER}, until {when(until)}. "
            "Guarded acts will be refused and queued.\n"
        )
        previous = hm.read_mode(state, now)
        if previous.state == "on":
            tty.write(f"This replaces the window ending {previous.until}.\n")
        tty.write("Type y to confirm: ")
        tty.flush()
        if tty.readline().strip().lower() not in ("y", "yes"):
            print("Not confirmed; nothing changed.")
            return 1
    try:
        battery = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True)
        if battery.returncode == 0 and "Battery Power" in battery.stdout:
            print(
                "On battery: caffeinate -s holds only on AC power. Plug in, or the Mac may sleep."
            )
    except OSError:
        pass
    if old is not None:
        stop(old.get("keep_awake_pid"))
    argv = [part.format(seconds=seconds) for part in hm.KEEP_AWAKE]
    try:
        pid = spawn(
            argv,
            start_new_session=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).pid
        awake = f"caffeinate pid {pid}"
    except OSError as error:
        pid, awake = None, f"not started ({error})"
    since_iso, until_iso = now.isoformat(timespec="seconds"), until.isoformat(timespec="seconds")
    flag = {"since": since_iso, "until": until_iso, "set_by": "away", "keep_awake_pid": pid}
    temporary = None
    try:  # one unit: on failure, undo step 5 and write nothing (§3.3 step 6)
        state.mkdir(mode=0o700, parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(dir=state)  # created 0o600
        with os.fdopen(fd, "w") as out:
            json.dump(flag, out)
        os.replace(temporary, state / "unattended.json")
    except OSError as error:
        if pid is not None:
            stop(pid)
        if temporary is not None:
            with contextlib.suppress(OSError):
                os.unlink(temporary)
        print(
            f"away.py could not write the flag ({type(error).__name__}: {error}); "
            "keep-awake stopped, nothing changed.",
            file=sys.stderr,
        )
        return 1
    append(state / "windows.jsonl", {"event": "enter", "since": since_iso, "until": until_iso})
    print(
        f"UNATTENDED until {when(until)}. Keep-awake: {awake}. Guarded acts are refused and "
        "queued. End early: python3 tools/away.py --back"
    )
    return 0


def back(root: Path, state: Path, now: datetime, stop: Callable[[object], None]) -> int:
    flag = read_flag(state)
    ended = now.isoformat(timespec="seconds")
    if flag is None:
        print("No unattended flag was set.")
    else:
        stop(flag.get("keep_awake_pid"))
        path = state / "unattended.json"
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)  # a directory flag is broken, and --back clears any state
        else:
            path.unlink()
        append(
            state / "windows.jsonl", {"event": "back", "since": flag.get("since"), "until": ended}
        )
        print(
            f"Back. Unattended since {flag.get('since')} until {flag.get('until')} (ended {ended})."
        )
        print(stretch(root, flag, now))
    queue = state / "queue.jsonl"
    entries = []
    if queue.exists():
        for line in queue.read_text(errors="replace").splitlines():
            try:
                entries.append(json.loads(line))
            except ValueError:
                continue
    print(f"Queued while away ({len(entries)}):")
    branches: dict[str, list[dict[str, Any]]] = {}
    for entry in entries:
        branches.setdefault(str(entry.get("branch") or "(no branch)"), []).append(entry)
    for branch, queued in branches.items():
        print(f"  {branch}")
        for e in queued:
            print(
                f"    {e.get('at')} {e.get('hook')} {e.get('agent_type') or 'main'}: {e.get('act')}"
            )
    if not entries:
        print("  (none)")
    try:  # the recap module is optional to --back; the flag and the queue are not (§3.3)
        from session_state import all_decisions

        decisions = all_decisions(root)
    except Exception as error:
        decisions = [f"(unavailable: {type(error).__name__}: {error})"]
    print("ASK OLA lines, every worktree's .claude/current-task/:")
    for line in decisions or ["(none)"]:
        print(f"  {line}")
    if queue.exists():
        archive = state / f"queue-{now.astimezone(UTC):%Y-%m-%d}.jsonl"
        if archive.exists():
            with archive.open("a") as out:
                out.write(queue.read_text())
            queue.unlink()
        else:
            queue.rename(archive)
    return 0


def main(
    argv: list[str] | None = None,
    *,
    root: Path | None = None,
    now: datetime | None = None,
    open_tty: Callable[[], TextIO] | None = None,
    spawn: Callable[..., Any] | None = None,
    stop: Callable[[object], None] | None = None,
) -> int:
    parser = argparse.ArgumentParser(usage=USAGE)
    parser.add_argument("duration", nargs="?")
    parser.add_argument("--back", action="store_true")
    args = parser.parse_args(argv)
    root = root or Path(__file__).resolve().parents[1]
    now = now or datetime.now(UTC)
    stop = stop or stop_keep_awake
    state = hm.state_dir(root)
    minutes = minutes_of(args.duration)
    if not args.back and minutes is None:
        print(USAGE, file=sys.stderr)
        return 2
    if state is None:
        print("away.py: git common dir not found.", file=sys.stderr)
        return 1
    if args.back or minutes is None:
        return back(root, state, now, stop)
    return enter(
        args.duration,
        minutes,
        state,
        now,
        # A tty cannot seek, so no BufferedRandom: one FileIO per direction.
        # The ignore: typeshed's wrapped-buffer protocol wants a `name` that
        # BufferedRWPair lacks; TextIOWrapper reads it only on `.name`.
        open_tty
        or (
            lambda: io.TextIOWrapper(  # type: ignore[type-var]
                io.BufferedRWPair(io.FileIO("/dev/tty", "r"), io.FileIO("/dev/tty", "w")),
                line_buffering=True,
            )
        ),
        spawn or subprocess.Popen,
        stop,
    )


if __name__ == "__main__":
    sys.exit(main())
