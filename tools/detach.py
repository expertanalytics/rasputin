#!/usr/bin/env python3
"""Launch a heavy run as a detached job and re-attach to it from disk.

The increment is docs/increments/h40-detached-runs.md. A heavy run -- a `@perf`
timed acceptance / scaling / profile run, or an lsm model run -- must never live
on an agent's monitored stream: a laptop suspend or an agent kill silences the
stream and the Anthropic harness watchdog kills the run. The fix is standard
Unix job detachment: launch the command in its own session (`setsid`, here
`subprocess.Popen(start_new_session=True)`), redirect its stdio to files, write
a small manifest, and return at once. The launching persona hands back; a later
persona re-attaches by reading the manifest and the job's own completion
artifact. Liveness no longer depends on the launcher's process, group,
controlling terminal or stream.

Two verbs:

* `start --run-dir DIR --name NAME [--done-marker FILE] [--force] -- cmd...`
  launches the job and writes `DIR/NAME.detached.json`.
* `status --run-dir DIR [--name NAME]` reads each `*.detached.json` manifest and
  reports, per job, one of four states:
  - `done`    -- the manifest's done-marker file exists and parses as JSON;
  - `running` -- no done-marker, the pid is alive, and its start time matches
                 the manifest within START_TOLERANCE_S (the pid-reuse guard);
  - `failed`  -- no done-marker and the pid is gone, or a reused pid whose start
                 time does not match (a foreign process);
  - `unknown` -- the manifest was written on another host, so liveness cannot be
                 judged here; the host guard wins over any local kill(pid, 0).

The done-marker is the job's own final artifact, written last and atomically
(write `<name>.tmp` in the same directory, then `os.replace`). A `bench.py run`
writes `run.json` only on success, so that is its natural done-marker; a
truncated or `.tmp` marker never reads as done -- the marker must parse.

`status` is a pure read of disk plus one `kill(pid, 0)` probe; it starts no
heavy work and is safe from any persona or the main session. Standard library
only; no file decoding, no CRS, nothing that crosses the C++ I/O boundary.

Usage:
    python3 tools/detach.py start --run-dir DIR --name NAME -- cmd [args...]
    python3 tools/detach.py status --run-dir DIR [--name NAME]
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

State = Literal["running", "done", "failed", "unknown"]

# The pid-reuse clock tolerance (seconds). Start times are compared normalized
# to UTC and must agree within this bound. It is >= 1 s because `ps -o lstart=`
# is 1 s granular (macOS) and `/proc` start time is clock-tick granular (Linux);
# a tighter bound would flap on those sources. Fixed at 2 s (@tester, T6): one
# tick of slack over the 1 s granularity floor, far below any pid-reuse gap.
START_TOLERANCE_S = 2.0

_MARKER = ".detached.json"


def this_host() -> str:
    """The host name a manifest records so another machine reads `unknown`."""
    return platform.node()


def atomic_write_json(path: Path, obj: object) -> None:
    """Write `obj` as JSON to `path` atomically: a `.tmp` sibling, then replace.

    The temp and final names share a directory so `os.replace` is atomic
    (cross-filesystem rename is not); an interrupted write leaves only the
    `.tmp`, which no reader ever treats as the final file.
    """
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1) + "\n")
    os.replace(tmp, path)


def pid_start_utc(pid: int) -> str | None:
    """The process's start time as an ISO-8601 UTC string, or None if unreadable.

    Linux reads `/proc/<pid>/stat` field 22 (clock ticks since boot) against the
    boot time in `/proc/stat`. macOS (and other Unix) parses `ps -o lstart=`,
    which prints local time at 1 s granularity. Both are normalized to UTC so a
    manifest written on this machine and a live pid compare on one clock.
    """
    if sys.platform.startswith("linux"):
        return _linux_start_utc(pid)
    return _ps_start_utc(pid)


def _linux_start_utc(pid: int) -> str | None:
    try:
        boot = None
        for line in Path("/proc/stat").read_text().splitlines():
            if line.startswith("btime "):
                boot = int(line.split()[1])
        stat = Path(f"/proc/{pid}/stat").read_text()
        # Field 22 (1-based) is starttime in clock ticks; fields after comm,
        # which is parenthesized and may contain spaces, so split past the ')'.
        fields = stat[stat.rindex(")") + 2 :].split()
        ticks = int(fields[19])
        hz = os.sysconf("SC_CLK_TCK")
        if boot is None or hz <= 0:
            return None
        epoch = boot + ticks / hz
        return datetime.fromtimestamp(epoch, UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    except (OSError, ValueError, IndexError):
        return None


def _ps_start_utc(pid: int) -> str | None:
    try:
        out = subprocess.run(
            ["ps", "-o", "lstart=", "-p", str(pid)],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    text = out.stdout.strip()
    if out.returncode or not text:
        return None
    # `ps -o lstart=` prints e.g. "Mon Nov  2 11:20:11 2026" in local time.
    local = datetime.strptime(re.sub(r"\s+", " ", text), "%a %b %d %H:%M:%S %Y")
    return local.astimezone().astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _parse_utc(stamp: str) -> datetime | None:
    try:
        return datetime.strptime(stamp, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=UTC)
    except ValueError:
        return None


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _marker_done(run_dir: Path, manifest: dict[str, object]) -> bool:
    """True when the manifest's done-marker exists and parses as JSON.

    A truncated or `.tmp` marker fails to parse (or is a different name) and so
    reads not-done; pid-liveness then decides running vs failed.
    """
    name = manifest.get("done_marker")
    if not isinstance(name, str):
        return False
    marker = run_dir / name
    try:
        json.loads(marker.read_text())
    except (OSError, ValueError):
        return False
    return True


def _state(run_dir: Path, manifest: dict[str, object]) -> State:
    if _marker_done(run_dir, manifest):
        return "done"
    if manifest.get("host") != this_host():
        return "unknown"  # another machine: the host guard wins over kill(pid,0)
    pid = manifest.get("pid")
    if not isinstance(pid, int) or not _alive(pid):
        return "failed"
    # Pid-reuse guard: the live pid must be the original job, not a reused pid.
    started = _parse_utc(str(manifest.get("started", "")))
    actual = pid_start_utc(pid)
    now = _parse_utc(actual) if actual else None
    if started is None or now is None:
        return "running"  # cannot read a start time to compare; pid is alive
    if abs((now - started).total_seconds()) > START_TOLERANCE_S:
        return "failed"  # a foreign process reused the pid
    return "running"


def status(run_dir: Path, name: str | None = None) -> dict[str, State]:
    """Each job's state under `run_dir`, keyed by name. Missing dir -> empty."""
    states: dict[str, State] = {}
    if not run_dir.is_dir():
        return states
    for path in sorted(run_dir.glob(f"*{_MARKER}")):
        job = path.name[: -len(_MARKER)]
        if name is not None and job != name:
            continue
        try:
            manifest = json.loads(path.read_text())
        except (OSError, ValueError):
            continue  # a malformed manifest is not a job; status() does not crash
        states[job] = _state(run_dir, manifest)
    return states


def start(
    run_dir: Path,
    name: str,
    cmd: list[str],
    done_marker: str,
    *,
    force: bool = False,
) -> int:
    """Launch `cmd` detached, write its manifest, and return the child pid.

    The child lands in its own session (start_new_session=True) with stdin from
    /dev/null (so it never blocks on input) and stdout/stderr to files under
    `run_dir`. The manifest is written atomically after the launch. A second
    `start` with the same name refuses with FileExistsError unless `force`.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = run_dir / f"{name}{_MARKER}"
    if manifest_path.exists() and not force:
        raise FileExistsError(f"detached job {name!r} already exists: {manifest_path}")
    stdout = run_dir / f"{name}.stdout"
    stderr = run_dir / f"{name}.stderr"
    with open(os.devnull, "rb") as devnull, stdout.open("wb") as out, stderr.open("wb") as err:
        proc = subprocess.Popen(
            cmd,
            start_new_session=True,
            stdin=devnull,
            stdout=out,
            stderr=err,
            cwd=os.getcwd(),  # the launcher's repo/worktree root (design §"start")
        )
    manifest = {
        "name": name,
        "pid": proc.pid,
        "pgid": os.getpgid(proc.pid),
        "cmd": list(cmd),
        "cwd": str(Path.cwd()),
        "started": pid_start_utc(proc.pid) or datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "host": this_host(),
        "stdout": stdout.name,
        "stderr": stderr.name,
        "done_marker": done_marker,
        "state": "running",
    }
    atomic_write_json(manifest_path, manifest)
    return proc.pid


def _print_status(run_dir: Path, states: dict[str, State]) -> None:
    for job, state in states.items():
        if state == "done":
            path = run_dir / json.loads((run_dir / f"{job}{_MARKER}").read_text())["done_marker"]
            print(f"done {job} {path}")
        elif state == "running":
            pid = json.loads((run_dir / f"{job}{_MARKER}").read_text())["pid"]
            print(f"running {job} pid {pid}")
        elif state == "failed":
            print(f"failed {job} (no done-marker, process gone)")
        else:
            host = json.loads((run_dir / f"{job}{_MARKER}").read_text())["host"]
            print(f"unknown {job} (host {host}, not this machine)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Launch and re-attach detached heavy runs.")
    sub = parser.add_subparsers(dest="verb", required=True)
    s = sub.add_parser("start", help="launch a detached job")
    s.add_argument("--run-dir", required=True, type=Path)
    s.add_argument("--name", required=True)
    s.add_argument("--done-marker", default="run.json")
    s.add_argument("--force", action="store_true")
    s.add_argument("cmd", nargs=argparse.REMAINDER)
    t = sub.add_parser("status", help="report each detached job's state")
    t.add_argument("--run-dir", required=True, type=Path)
    t.add_argument("--name", default=None)
    args = parser.parse_args(argv)

    if args.verb == "start":
        cmd = args.cmd[1:] if args.cmd and args.cmd[0] == "--" else args.cmd
        if not cmd:
            parser.error("no command given after --")
        pid = start(args.run_dir, args.name, cmd, args.done_marker, force=args.force)
        print(args.run_dir / f"{args.name}{_MARKER}")
        print(f"pid {pid}")
        return 0
    states = status(args.run_dir, args.name)
    if not states:
        print("none")
    else:
        _print_status(args.run_dir, states)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
