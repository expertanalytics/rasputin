#!/usr/bin/env python3
"""List, and with `--kill` stop, processes this project left behind.

    python3 tools/reap.py [--kill] [--worktree <path>]

A process is listed, numbered with pid, age and command line, when its working
folder is (a) inside a worktree that was removed, or one whose branch is merged
into origin/master (its tip reachable from it, but not on its first-parent
line, so a fresh branch cut from master is not "merged"), or (b) inside the
repository or any of its worktrees and it has run for over 6 hours. Never
listed: this process, its ancestors, `claude` sessions, and anything on a
terminal (a person's shell; Bash-tool jobs have none). `--worktree`
limits both to that one worktree. `--kill` sends SIGTERM to exactly the listed
processes, SIGKILL to any still alive 5 s later, and prints what it stopped.
Working folders come from /proc on Linux and `lsof -d cwd` elsewhere (macOS).

Why: on 2026-10-30 one process waited 23 h on input after its branch merged
(#400), and four ran 4 days after their worktree was removed.
Tests: tests/python/test_reap.py. Standard library only.
"""

from __future__ import annotations

import argparse
import contextlib
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import NamedTuple

LIMIT = 6 * 3600
GRACE = 5.0


class Proc(NamedTuple):
    pid: int
    ppid: int
    age: int  # seconds
    cwd: str
    command: str
    tty: str = ""  # `ps`'s terminal; "?" or "??" for none


class Tree(NamedTuple):
    path: str
    merged: bool


def _inside(path: str, folder: str) -> bool:
    folder = folder.rstrip("/")
    return path == folder or path.startswith(folder + "/")


def select(
    procs: list[Proc], trees: list[Tree], repo: str, *, me: int, only: str | None = None
) -> list[tuple[Proc, str]]:
    """The processes to list, each with why; pure over the given tables."""
    parent = {p.pid: p.ppid for p in procs}
    spared, pid = {me}, me
    while pid in parent and parent[pid] not in spared:
        pid = parent[pid]
        spared.add(pid)
    found = []
    for p in procs:
        cwd = p.cwd.removesuffix(" (deleted)")
        human = p.tty.strip("?") or os.path.basename((p.command.split() or [""])[0]) == "claude"
        if p.pid in spared or human:
            continue
        if only is not None and not _inside(cwd, only):
            continue
        home = [t for t in trees if _inside(cwd, t.path) and t.path != repo]
        if p.cwd.endswith(" (deleted)") or (not home and _inside(cwd, f"{repo}/.claude/worktrees")):
            found.append((p, "worktree removed"))
        elif home and max(home, key=lambda t: len(t.path)).merged:
            found.append((p, "branch merged"))
        elif p.age > LIMIT and any(_inside(cwd, t.path) for t in [Tree(repo, False), *trees]):
            found.append((p, "over 6 h"))
    return found


def seconds(etime: str) -> int:
    """`ps`'s elapsed time, `[[dd-]hh:]mm:ss`, in seconds."""
    days, _, clock = etime.rpartition("-")
    total = 0
    for part in clock.split(":"):
        total = total * 60 + int(part)
    return total + int(days or 0) * 86400


def listing(found: list[tuple[Proc, str]]) -> list[str]:
    lines = []
    for n, (p, why) in enumerate(found, 1):
        days, rest = divmod(p.age, 86400)
        age = f"{days}d {rest // 3600:02}:{rest % 3600 // 60:02}"
        lines.append(f"{n}. pid {p.pid}, {age}, {why}, in {p.cwd}: {p.command[:200]}")
    return lines or ["none"]


def _git(repo: Path, *args: str) -> str:
    done = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    return done.stdout if done.returncode == 0 else ""


def trees(repo: Path) -> list[Tree]:
    """The live worktrees (prunable ones are gone), each with whether its branch merged."""
    base = (
        "origin/master" if _git(repo, "rev-parse", "-q", "--verify", "origin/master") else "master"
    )
    reached = set(_git(repo, "rev-list", base).split())
    trunk = set(_git(repo, "rev-list", "--first-parent", base).split())
    found = []
    for block in _git(repo, "worktree", "list", "--porcelain").split("\n\n"):
        fields = dict(line.partition(" ")[::2] for line in block.splitlines())
        if "worktree" in fields and "prunable" not in fields:
            head = fields.get("HEAD", "")
            found.append(Tree(fields["worktree"], head in reached and head not in trunk))
    return found


def processes() -> list[Proc]:
    ps = subprocess.run(["ps", "-e", "-ww", "-o", "pid=,ppid=,etime=,tty=,command="],
                        capture_output=True, text=True, check=True).stdout  # fmt: skip
    cwd: dict[int, str] = {}
    if Path("/proc/self/cwd").exists():
        for entry in os.listdir("/proc"):
            with contextlib.suppress(OSError, ValueError):
                cwd[int(entry)] = os.readlink(f"/proc/{entry}/cwd")
    else:
        lsof = subprocess.run(["lsof", "-d", "cwd", "-Fn"], capture_output=True, text=True).stdout
        owner = 0
        for line in lsof.splitlines():
            if line[:1] == "p":
                owner = int(line[1:])
            elif line[:1] == "n":
                cwd[owner] = line[1:]
    procs = []
    for line in ps.splitlines():
        parts = line.split(None, 4)
        if len(parts) == 5 and int(parts[0]) in cwd:
            pid, ppid, age, tty, command = parts
            procs.append(Proc(int(pid), int(ppid), seconds(age), cwd[int(pid)], command, tty))
    return procs


def found(only: str | None = None) -> list[tuple[Proc, str]]:
    """The listing for this machine now; the recap prints `listing(found())`."""
    here = Path(__file__).resolve().parent.parent
    common = _git(here, "rev-parse", "--path-format=absolute", "--git-common-dir").strip()
    repo = Path(common).parent if common.endswith("/.git") else here
    only = str(Path(only).resolve()) if only else None
    return select(processes(), trees(repo), str(repo), me=os.getpid(), only=only)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--kill", action="store_true")
    parser.add_argument("--worktree")
    args = parser.parse_args(argv)
    listed = found(args.worktree)
    print("\n".join(listing(listed)))
    if args.kill and listed:
        stop([p.pid for p, _ in listed])
    return 0


def _signal(pid: int, sig: int) -> bool:
    """Send `sig` (0 only probes); False when the process is gone."""
    try:
        os.kill(pid, sig)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass
    return True


def stop(pids: list[int]) -> None:
    """SIGTERM, then SIGKILL for whatever is still alive after GRACE seconds."""
    alive = [pid for pid in pids if _signal(pid, signal.SIGTERM)]
    deadline = time.monotonic() + GRACE
    while alive and time.monotonic() < deadline:
        time.sleep(0.2)
        alive = [pid for pid in alive if _signal(pid, 0)]
    killed = [pid for pid in alive if _signal(pid, signal.SIGKILL)]
    print(f"stopped by SIGTERM: {' '.join(str(p) for p in pids if p not in alive) or 'none'}")
    if killed:
        print(f"stopped by SIGKILL after {GRACE:g} s: {' '.join(map(str, killed))}")


if __name__ == "__main__":
    sys.exit(main())
