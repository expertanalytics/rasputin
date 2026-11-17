#!/usr/bin/env python3
"""PreToolUse on Bash: a kill that can hit the wrong process is refused.

On 2026-11-17 a hand-written `kill <pid>` hit the session's own Claude Code
process, taken for a stale one; on 2026-11-07 a `pkill -f` hit Ola's
applications and other agents' runs. Refused, as a `deny` naming the reason:
every `pkill` and `killall`; a `kill` of a process group (`0` or a negative
pid); a `kill` of this hook's own ancestor (the Claude Code process and its
shells) or of a `claude` process. Any other pid passes, and so does one given
through a variable or `$(...)`, which the hook cannot read. `tools/reap.py`
excludes the same processes. Each simple command is read through
tools/shell_scan.py, so `kill` in a chain, in `$(...)`, after `xargs` (pids
from an `echo` before the pipe) or in `sh -c` is judged; a string is not.
Tests: tests/python/test_guard_kill.py. Standard library only.
"""

import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2] / "tools"))  # appended (h16 G4)
try:
    import shell_scan  # type: ignore[import-not-found,unused-ignore]
except ImportError:  # each line is then split on its separators
    shell_scan = None

CLAUDE = re.compile(r"(^|/)claude(\.exe)?(\s|$)")
INSTEAD = "use `python3 tools/reap.py --kill`, or `kill <pid>` of a pid you started"


def name(word: str) -> str:
    return word.rsplit("/", 1)[-1]


def targets(args: list[str]) -> list[str]:
    """The pid words of `kill args`: the first option is the signal, later ones are pids."""
    found: list[str] = []
    signal, at = False, 0
    while at < len(args):
        word = args[at]
        if word == "--":
            return found + args[at + 1 :]
        if word in ("-l", "-L"):
            return []
        if word.startswith("-") and not signal:
            signal, at = True, at + (2 if word in ("-s", "-n") else 1)
            continue
        found.append(word)
        at += 1
    return found


def refusal(pids: list[str]) -> str | None:
    table = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,comm="], capture_output=True,
                           text=True, check=True, timeout=10).stdout  # fmt: skip
    procs = {int(p): (int(pp), c) for p, pp, c in (ln.split(None, 2) for ln in table.splitlines()
                                                   if len(ln.split(None, 2)) == 3)}  # fmt: skip
    ancestors, pid = set(), os.getpid()
    while pid in procs and pid not in ancestors and pid > 1:
        ancestors.add(pid)
        pid = procs[pid][0]
    for word in pids:
        if not re.fullmatch(r"-?\d+", word):
            continue  # `%1`, `$PID`, `$(...)`: not readable here
        if int(word) <= 0:
            return f"`kill {word}` signals a whole process group"
        if int(word) in ancestors:
            return f"pid {word} is this session's own ancestor ({procs[int(word)][1]})"
        if int(word) in procs and CLAUDE.search(procs[int(word)][1]):
            return f"pid {word} is a claude process ({procs[int(word)][1]})"
    return None


def why(command: str) -> str | None:
    simples = shell_scan.parse(command) if shell_scan else None
    pairs: list[tuple[list[str], list[str]]] = []
    if simples is None:  # unreadable: split on separators; a quote error fails closed in main
        pairs = [(shlex.split(s), []) for s in re.split(r"&&|\|\||;|\||\n", command)]
    for at, s in enumerate(simples or []):
        fed = at > 0 and s.piped and "xargs" in map(name, s.written)  # `echo 12 | xargs kill`
        pairs.append((s.argv, simples[at - 1].argv[1:] if simples and fed else []))
    pids: list[str] = []
    for argv, fed in pairs:
        if argv and name(argv[0]) in ("pkill", "killall"):
            return f"`{name(argv[0])}` matches by name and can hit any process; {INSTEAD}"
        if argv and name(argv[0]) == "kill":
            pids += targets(argv[1:]) + fed
    return refusal(pids) if pids else None


def main() -> int:
    try:
        event = json.load(sys.stdin)
        command = (event.get("tool_input") or {}).get("command", "")
        reason = why(command) if event.get("tool_name") == "Bash" else None
    except Exception as error:  # fail closed on a line that kills, open on any other
        text = locals().get("command", "")
        reason = (f"guard_kill failed ({type(error).__name__}: {error})"
                  if re.search(r"\b(p?kill|killall)\b", str(text)) else None)  # fmt: skip
    if reason:
        print(json.dumps({"hookSpecificOutput": {
            "hookEventName": "PreToolUse", "permissionDecision": "deny",
            "permissionDecisionReason": f"Refused: {reason}. {INSTEAD}."}}))  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
