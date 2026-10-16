#!/usr/bin/env python3
"""Print the round recap, then what the previous session was in the middle of.

The recap (retrospective rule 5, Ola): the last thing landed, what is in
flight, running Bash-tool jobs, decisions waiting on Ola (from every worktree),
format warnings for `session.md`, the next three ROADMAP items and the rule-text
size table -- all read from files on disk, git and `ps`, so it is the same
whoever runs it (docs/increments/h8-window-and-recap.md §3). The task files are
the main checkout's, found from the common git dir, so a session started inside
a worktree sees them too.

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
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Literal

REPO = Path(__file__).resolve().parent.parent
sys.path.append(str(REPO / "tools"))  # appended: the stdlib wins (h16 G4)
if TYPE_CHECKING:  # imported in print_recap, where a harness fault costs one line (§3.8)
    import harness_mode
UNKNOWN = (
    "UNATTENDED STATE UNKNOWN ({}): the recap could not read the flag or the queue. "
    "Guarded acts may be refused; record each refusal as an ASK OLA line and continue."
)

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


def main_checkout(repo: Path) -> Path:
    """The main checkout: the parent of the common git dir when that is `.git`."""
    common = _git(repo, "rev-parse", "--path-format=absolute", "--git-common-dir")
    return Path(common).parent if common and Path(common).name == ".git" else repo


def checkouts(repo: Path) -> list[Path]:
    """The main checkout, then every worktree whose directory still exists."""
    listed = _git(repo, "worktree", "list", "--porcelain")
    if listed is None:
        return [repo]
    paths = [Path(line[9:]) for line in listed.splitlines() if line.startswith("worktree ")]
    return [path for path in paths if path.is_dir()]


# h8 §3.2: at the start of a line, after an optional bullet; case-sensitive.
ASK = re.compile(r"^[ \t]*(?:[-*+][ \t]+)?ASK OLA:(.*)$")


def pending_decisions(tasks: Path) -> list[str]:
    """`ASK OLA:` lines in any current-task file; an empty one as a warning."""
    if not tasks.is_dir():
        return []
    found: list[str] = []
    for path in sorted(p for p in tasks.glob("*") if not p.is_dir()):
        for number, line in enumerate(_read(path, "").splitlines(), 1):
            if match := ASK.match(line):
                empty = f"{path.name}:{number}: WARNING: empty ASK OLA line"
                found.append(f"{path.name}: {line.strip()}" if match[1].strip() else empty)
    return found


def all_decisions(repo: Path) -> list[str]:
    """pending_decisions of every checkout; a worktree's lines carry its path."""
    trees = checkouts(repo)
    found: list[str] = []
    for tree in trees:
        name = tree.relative_to(trees[0]) if tree.is_relative_to(trees[0]) else tree
        prefix = "" if tree == trees[0] else f"{name}: "
        found += [prefix + line for line in pending_decisions(tree / ".claude" / "current-task")]
    return found


KIND = re.compile(r"^[ \t]*(?:[-*+][ \t]+)?(NOW|QUEUE|ASK OLA):(.*)$")
MAX_LINE = 300


def session_format(text: str, name: str = "session.md") -> list[str]:
    """Warnings on `session.md`: one NOW, one QUEUE, ASK OLA lines, nothing else."""
    counts, empty, other, long = {"NOW": 0, "QUEUE": 0}, [], [], []
    for number, line in enumerate(text.splitlines(), 1):
        line = line.rstrip()
        if not line:
            continue
        if len(line) > MAX_LINE:
            long.append(f"{name}:{number}: {len(line)} characters; at most {MAX_LINE}")
        match = KIND.match(line)
        if match is None:
            other.append(f"{name}:{number}: not a NOW, QUEUE or ASK OLA line: {line[:60]}")
        elif match[1] != "ASK OLA":  # an empty ASK OLA line is warned by pending_decisions
            counts[match[1]] += 1
            if not match[2].strip():
                empty.append(f"{name}:{number}: empty {match[1]} line")
    wrong = [
        f"{name}: {n} {kind} lines; exactly one expected" for kind, n in counts.items() if n != 1
    ]
    return wrong + empty + other + long


def _capped(lines: list[str], cap: int) -> list[str]:
    return lines[:cap] + ([f"... {len(lines) - cap} more"] if len(lines) > cap else [])


def _is_wrapper(command: str) -> bool:
    """Claude Code's Bash-tool shell: sources a shell snapshot and evals the command."""
    return "/.claude/shell-snapshots/" in command and "eval '" in command


def background_jobs(ps_text: str, exclude: set[int]) -> list[str]:
    """The Bash-tool wrappers in a `ps` listing, at most 4, then how many more."""
    jobs = []
    for line in ps_text.splitlines():
        parts = line.split(None, 3)
        if len(parts) < 4 or not parts[0].isdigit() or int(parts[0]) in exclude:
            continue
        if _is_wrapper(parts[3]):
            command = parts[3].split("eval '", 1)[1].split("' < /dev/null", 1)[0]
            jobs.append(f"pid {parts[0]}, running {parts[2]}: {' '.join(command.split())[:100]}")
    return _capped(jobs, 4)


def wrapper_check(ancestors: list[str]) -> Literal["confirmed", "drift", "unknown"]:
    """Whether the ancestors, parent upward, show Claude Code's wrapper format (§3.5)."""
    drift = False
    for command in ancestors:
        tokens = command.split()
        if not tokens:
            continue
        if os.path.basename(tokens[0]) in ("claude", "claude.exe"):
            return "drift" if drift else "unknown"
        if _is_wrapper(command):
            return "confirmed"
        shell = os.path.basename(tokens[0]) in ("sh", "bash", "zsh") and tokens[1:2] == ["-c"]
        hook = (
            tokens[2:3] != []
            and tokens[2].startswith("python")
            and "tools/session_state.py" in command
        )
        drift = drift or (shell and not hook)
    return "unknown"


def running_jobs() -> list[str]:
    """The jobs section's lines: from `ps`, without this process and its ancestors."""
    try:
        ps = subprocess.run(
            ["ps", "-e", "-ww", "-o", "pid=,ppid=,etime=,command="],
            capture_output=True,
            text=True,
            check=False,
        )
        if ps.returncode:
            raise OSError(ps.stderr.strip() or f"ps exited {ps.returncode}")
        table = {}
        for line in ps.stdout.splitlines():
            parts = line.split(None, 3)
            if len(parts) == 4 and parts[0].isdigit() and parts[1].isdigit():
                table[int(parts[0])] = (int(parts[1]), parts[3])
        pid, exclude, chain = os.getpid(), {os.getpid()}, []
        while pid in table and table[pid][0] not in exclude:
            pid = table[pid][0]
            exclude.add(pid)
            chain += [table[pid][1]] if pid in table else []
        jobs = background_jobs(ps.stdout, exclude)
        drift = ["(wrapper format not recognised)"] if wrapper_check(chain) == "drift" else []
        return jobs + drift or ["(none)"]
    except Exception as error:
        return [f"(could not list: {str(error)[:100]})"]


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


def unattended_header(mode: harness_mode.Mode) -> str | None:
    """The line above the recap while the flag is not off (h3 U1 §3.8)."""
    if mode.state == "broken":
        return (
            f"UNATTENDED FLAG UNREADABLE ({mode.detail}): every guarded act is refused. "
            "Ola: python3 tools/away.py --back."
        )
    if mode.until is None:
        return None
    until = f"{mode.until.astimezone(UTC):%Y-%m-%d %H:%M} UTC"
    if mode.state == "on":
        return (
            f"UNATTENDED until {until}. Guarded acts are refused and queued: "
            "record each as an ASK OLA line and continue."
        )
    return (
        f"Unattended mode ended at {until}. "
        "Ola: python3 tools/away.py --back prints and archives the queue."
    )


def queued(state: Path | None) -> list[str]:
    """The newest 10 queue lines, then how many are older or unreadable."""
    path = state / "queue.jsonl" if state else None
    if path is None or not path.is_file():
        return []
    entries, bad = [], 0
    for line in _read(path, "").splitlines():
        try:
            entry = json.loads(line)
            entries.append(
                f"  {entry['at']} {entry.get('branch') or '(no branch)'} {entry['hook']} "
                f"{entry.get('agent_type') or 'main'}: {str(entry['act'])[:100]}"
            )
        except (ValueError, KeyError, TypeError):
            bad += 1
    lines = [f"Queued while unattended ({len(entries)}):", *entries[-10:]]
    if len(entries) > 10:
        lines.append(f"  ... {len(entries) - 10} older")
    if bad:
        lines.append(f"  ({bad} unreadable lines)")
    return lines


def print_recap(tasks: Path) -> None:
    """Retrospective rule 5: the structured recap every round opens with."""
    try:
        import harness_mode

        state = harness_mode.state_dir(REPO)
        header = unattended_header(harness_mode.read_mode(state, datetime.now(UTC)))
        queue = queued(state)
    except Exception as error:
        header, queue = UNKNOWN.format(f"{type(error).__name__}: {error}"[:200]), []
    if header:
        print(header)
    print("== recap ==")
    print(f"Last landed: {last_landed(REPO)}")
    print("In flight:")
    flight = in_flight(REPO, default_base(REPO))
    for line in flight[:11]:
        print(f"  {line}")
    if len(flight) > 11:
        print(f"  ... {len(flight) - 11} older")
    print("Running Bash-tool jobs, any session (Claude Code shell wrappers):")
    for line in running_jobs():
        print(f"  {line}")
    decisions = all_decisions(REPO)
    print("Waiting on Ola:" + ("" if decisions else " (none recorded as ASK OLA)"))
    for line in _capped([d if len(d) <= 160 else d[:157] + "..." for d in decisions], 5):
        print(f"  {line}")
    session = tasks / "session.md"
    warnings = session_format(_read(session, "")) if os.path.lexists(session) else []
    if warnings:
        print("session.md format:")
    for line in _capped(warnings, 4):
        print(f"  {line}")
    for line in queue:
        print(line)
    print("Next on ROADMAP.md:")
    for row in roadmap_next(_read(REPO / "ROADMAP.md", "")):
        print(f"  {row}")
    try:  # a fault in the size table costs one line, like the harness (h8 §3)
        import rule_sizes

        sizes = rule_sizes.report(REPO)
    except Exception as error:
        sizes = [f"(size table unavailable: {f'{type(error).__name__}: {error}'[:150]})"]
    print("\n".join(sizes))
    print()


def print_current_task(tasks: Path) -> None:
    """Print the session's ask first, then each subagent's as context.

    A cold session needs one of these promoted, not a flat dump: `session.md` is
    the record of what the round is for, and a subagent file is a fragment of it
    delegated. Both are printed because a subagent file may be the only trace of
    a step that died, but the order says which one to believe about the round.

    Sorted by name for stable output. session.md is promoted above the rest
    explicitly, so the order among subagent files decides nothing.
    """
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

    tasks = main_checkout(REPO) / ".claude" / "current-task"
    print_recap(tasks)
    print_current_task(tasks)

    # Claude Code exports CLAUDE_CODE_SESSION_ID; the older spelling is kept as a
    # fallback so the script still excludes the current session if that changes.
    here = os.environ.get("CLAUDE_CODE_SESSION_ID") or os.environ.get("CLAUDE_SESSION_ID", "")
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
