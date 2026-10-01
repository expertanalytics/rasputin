"""The unattended mode: read the flag, decide a guard's verdict, append the queue.

Spec: docs/increments/h3-unattended-u1.md §3.1, §3.2 and §3.4. The flag only
ever restricts: while it is on, or cannot be read, an `ask` becomes a `deny`
that is recorded in the queue for Ola. Guards compute a verdict from the event
alone and call `guard`; `away.py` is the only writer of the flag.
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any, Literal

#: Rasputin's values of the design's §7 profile keys, until profile.toml exists.
OWNER = "Ola"
MAX_HOURS = 72
BUFFER = 1.2
KEEP_AWAKE = ("caffeinate", "-is", "-t", "{seconds}")
FORGE_HOST = "github.com"

ROOT = Path(__file__).resolve().parents[1]
MAX_WINDOW = timedelta(hours=MAX_HOURS * BUFFER, seconds=60)
ACT_LIMIT = 1000

Verdict = Literal["pass", "ask", "deny"]

RECORDED = "The refusal is already recorded in the queue."
QUEUE_REASON = (
    "Refused: unattended mode is on until {until}; {act} waits for Ola ({why}). "
    + RECORDED
    + " Do not retry it and do not work around it. Main session: add one line "
    "`ASK OLA: {act}: <what you wanted and why>` to `.claude/current-task/session.md`, "
    "then continue with other work. Subagent: put that `ASK OLA:` line in your handback "
    "and carry on with the rest of your brief. If the brief cannot proceed without this "
    "act, hand back now."
)
BROKEN_FIRST = (
    "Refused: the unattended flag cannot be read ({detail}), so every guarded act is "
    "refused until Ola runs python3 tools/away.py --back; {act} waits for Ola ({why})."
)


@dataclass(frozen=True)
class Mode:
    state: Literal["off", "on", "broken"]
    until: datetime | None
    detail: str


@dataclass(frozen=True)
class Decision:
    kind: Verdict
    reason: str
    queue: bool  # the caller appends a queue entry


def state_dir(root: Path) -> Path | None:
    """`<git-common-dir>/harness` of the repository at `root`, or None."""
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "--path-format=absolute", "--git-common-dir"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return None
    return Path(result.stdout.strip()) / "harness"


def _moment(flag: dict[str, Any], key: str) -> datetime | str:
    """The flag's time at `key`, or the text naming why it is not one."""
    value = flag.get(key)
    if not isinstance(value, str):
        return f"{key} missing or not a string"
    try:
        moment = datetime.fromisoformat(value)
    except ValueError:
        return f"{key} is not an ISO time"
    return moment if moment.utcoffset() is not None else f"{key} has no UTC offset"


def read_mode(state: Path | None, now: datetime) -> Mode:
    """The mode, per the table of §3.2: first failure wins, and a failure is `broken`."""
    if state is None:
        return Mode("broken", None, "git common dir not found")
    path = state / "unattended.json"
    try:
        # lstat for absence, so a dangling symlink is broken rather than off.
        path.lstat()
        if not stat.S_ISREG(path.stat().st_mode):
            return Mode("broken", None, "the flag is not a regular file")
        flag = json.loads(path.read_text())
    except FileNotFoundError:
        return Mode("off", None, "")
    except (OSError, ValueError) as error:
        return Mode("broken", None, f"unreadable: {type(error).__name__}")
    if not isinstance(flag, dict):
        return Mode("broken", None, "not a JSON object")
    since, until = _moment(flag, "since"), _moment(flag, "until")
    for moment in (since, until):
        if isinstance(moment, str):
            return Mode("broken", None, moment)
    assert isinstance(since, datetime) and isinstance(until, datetime)
    if until <= since:
        return Mode("broken", None, "until is not after since")
    if until - since > MAX_WINDOW:
        return Mode("broken", None, f"window longer than {MAX_HOURS} h x {BUFFER}")
    if now >= until:
        return Mode("off", until, "expired")
    return Mode("on", until, "")


def decide(verdict: Verdict, reason: str, mode: Mode, act: str, why: str) -> Decision:
    """Map a guard's verdict and the mode to a decision (§3.4's table). Pure."""
    if verdict != "ask" or mode.state == "off":
        return Decision(verdict, reason, False)
    rest = QUEUE_REASON.format(until="", act=act, why=why).split(RECORDED, 1)[1]
    if mode.state == "broken":
        first = BROKEN_FIRST.format(detail=mode.detail, act=act, why=why)
        return Decision("deny", f"{first} {RECORDED}{rest}", True)
    assert mode.until is not None
    until = f"{mode.until.astimezone(UTC):%Y-%m-%d %H:%M} UTC"
    return Decision("deny", QUEUE_REASON.format(until=until, act=act, why=why), True)


def queue_entry(
    event: dict[str, Any], hook: str, act: str, why: str, now: datetime
) -> dict[str, Any]:
    """One queue line: where and by whom the act was refused, and why."""
    cwd = event.get("cwd") or os.getcwd()
    result = subprocess.run(
        ["git", "-C", str(cwd), "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "at": now.isoformat(timespec="seconds"),
        "branch": result.stdout.strip() or None if result.returncode == 0 else None,
        "cwd": cwd,
        "agent_type": event.get("agent_type"),
        "agent_id": event.get("agent_id"),
        "hook": hook,
        "act": act[:ACT_LIMIT],
        "why": why,
    }


def append_queue(state: Path | None, entry: dict[str, Any]) -> str | None:
    """Append one line with one write (under PIPE_BUF); the error text, or None."""
    if state is None:
        return "git common dir not found"
    try:
        state.mkdir(mode=0o700, parents=True, exist_ok=True)
        fd = os.open(state / "queue.jsonl", os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
        try:
            os.write(fd, (json.dumps(entry) + "\n").encode())
        finally:
            os.close(fd)
    except OSError as error:
        return f"{type(error).__name__}: {error}"
    return None


def settle(
    event: dict[str, Any], hook: str, verdict: Verdict, reason: str, act: str, why: str
) -> Decision:
    """Resolve the mode, decide, and queue when told to. A refusal never waits on the queue."""
    now = datetime.now(UTC)
    state = state_dir(ROOT)
    decision = decide(verdict, reason, read_mode(state, now), act, why)
    if not decision.queue:
        return decision
    error = append_queue(state, queue_entry(event, hook, act, why, now))
    if error is None:
        return decision
    lost = f"The refusal could not be recorded ({error}); your ASK OLA line is the only record."
    return Decision(decision.kind, decision.reason.replace(RECORDED, lost), True)


def pretool(kind: str, reason: str) -> dict[str, Any]:
    """A PreToolUse hook's output for `ask` or `deny`."""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": kind,
            "permissionDecisionReason": reason,
        }
    }


def guard(
    event: dict[str, Any], hook: str, verdict: Verdict, reason: str, act: str, why: str
) -> dict[str, Any] | None:
    """The one call a PreToolUse guard makes: its output, or None for a pass."""
    if verdict == "pass":
        return None
    decision = settle(event, hook, verdict, reason, act, why)
    return pretool(decision.kind, decision.reason)
