#!/usr/bin/env python3
"""AskUserQuestion, PermissionRequest, ConfigChange, Stop and PostToolUse while Ola is away.

Spec: docs/increments/h3-unattended-u1.md §3.7. While Ola is at the keyboard
this hook prints nothing, except a refusal when it crashes or the
harness_mode module is missing. In unattended mode (on, or a flag that cannot be
read) a question, a permission prompt or a configuration reload would wait for
him and stall the run, so each is refused and queued instead.

Stop (docs/increments/h26-away-stop.md §2): only in mode on, and only when
nothing is running, a main-session stop is refused once per turn and never
queued. A crash, a missing harness_mode or a broken flag lets a stop through.

PermissionRequest, while Ola is present: the prompt that opens is recorded in
`prompts.jsonl`, so `tools/away.py` can refuse to start while it is unanswered.
PostToolUse: in mode on, the main session's first tool call per window gets a
note that Ola is away and which branches are granted. Neither ever refuses.
"""

import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path


def act_of(event: dict) -> tuple[str, str] | None:
    """(act, why) for an event this hook refuses at night, or None."""
    name = event.get("hook_event_name")
    supplied = event.get("tool_input") or {}
    if name == "PreToolUse" and event.get("tool_name") == "AskUserQuestion":
        first = (supplied.get("questions") or [{}])[0].get("question", "")
        why = "a question needs Ola; write it as the ASK OLA line"
        return f"AskUserQuestion: {first[:200]}", why
    if name == "PermissionRequest":
        target = supplied.get("command") or supplied.get("file_path") or ""
        return f"{event.get('tool_name')}: {target[:200]}", "a permission prompt needs Ola"
    if name == "ConfigChange":
        act = f"ConfigChange {event.get('source')} {event.get('file_path')}"
        return act, "a configuration change needs Ola"
    return None


FAILED = (
    "guard_unattended failed ({error}), so this is refused rather than left for a prompt "
    "nobody may answer. It is not queued: record it as an ASK OLA line (main session: in "
    "session.md; subagent: in its handback) and continue."
)


def refusal(name: str, reason: str) -> dict | None:
    """The refusal in `name`'s shape, built without harness_mode (§3.7, §3.9)."""
    if name == "PreToolUse":
        return {"hookSpecificOutput": {"hookEventName": "PreToolUse",
                                       "permissionDecision": "deny",
                                       "permissionDecisionReason": reason}}
    if name == "PermissionRequest":
        return {"hookSpecificOutput": {"hookEventName": "PermissionRequest",
                                       "decision": {"behavior": "deny", "message": reason}}}
    if name == "ConfigChange":
        return {"decision": "block", "reason": reason}
    return None


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0
    if not isinstance(event, dict):
        return 0
    try:
        output = settle(event)
    except Exception as error:  # after parsing, a crash refuses rather than waits (§3.9)
        output = refusal(event.get("hook_event_name", ""), FAILED.format(
            error=f"{type(error).__name__}: {error}"))
    if output is not None:
        print(json.dumps(output))
    return 0


STOP_REASON = (
    "Unattended mode is on until {until} and nothing is running: no agent, no background "
    "command, no monitor, no scheduled wake-up. If this turn ends now, the session stays "
    "idle until Ola returns. Before ending it: write any open question as an `ASK OLA:` "
    "line in `.claude/current-task/session.md` and do not wait for the answer; start the "
    "next item that needs no ruling (the `QUEUE:` line), or watch work that will finish, "
    "for example `gh pr checks <pr> --watch` run in the background. If no such work "
    "exists, set the `NOW:` line to `NOW: idle: <why>` and end the turn; this hook does "
    "not refuse it a second time."
)


def stop_output(event: dict, harness_mode) -> dict | None:
    """Refuse one idle stop of the main session while unattended mode is on."""
    if event.get("stop_hook_active") is not False or event.get("agent_id"):
        return None
    if "background_tasks" not in event or event["background_tasks"] or event.get("session_crons"):
        return None
    mode = harness_mode.read_mode(harness_mode.state_dir(harness_mode.ROOT), datetime.now(UTC))
    if mode.state != "on":
        return None
    until = f"{mode.until.astimezone(UTC):%Y-%m-%d %H:%M} UTC"
    return {"decision": "block", "reason": STOP_REASON.format(until=until)}


NOTE = (
    "Unattended mode is on until {until}: Ola is away. Granted branches: {names}. Each may "
    "be pushed, opened as a pull request and enqueued without asking, once `@reviewer`'s "
    "recorded verdict approves its head; do not wait for Ola's yes on these. Every other "
    "question is an `ASK OLA:` line."
)


def window_note(event: dict, harness_mode) -> dict | None:
    """Tell the main session once per window that unattended mode is on."""
    if event.get("agent_id"):
        return None
    state = harness_mode.state_dir(harness_mode.ROOT)
    mode = harness_mode.read_mode(state, datetime.now(UTC))
    if mode.state != "on":
        return None
    try:
        flag = json.loads((state / "unattended.json").read_text())
    except (OSError, ValueError):
        return None
    import push_grant

    seen = {"session_id": event.get("session_id"), "since": flag.get("since")}
    noticed = state / "noticed.jsonl"
    if noticed.exists() and seen in push_grant.jsonl(noticed):
        return None
    harness_mode.append_queue(state, seen, "noticed.jsonl")
    names = ", ".join(push_grant.grants(flag)) or "none"
    until = f"{mode.until.astimezone(UTC):%Y-%m-%d %H:%M} UTC"
    return {"hookSpecificOutput": {"hookEventName": "PostToolUse",
                                   "additionalContext": NOTE.format(until=until, names=names)}}


def record_prompt(event: dict, act: str, harness_mode) -> None:
    """Append the prompt about to open to `prompts.jsonl`; a failed append is silent."""
    supplied = json.dumps(event.get("tool_input") or {}, sort_keys=True)
    entry = {"at": datetime.now(UTC).isoformat(timespec="seconds"),
             "session_id": event.get("session_id"), "agent_id": event.get("agent_id"),
             "agent_type": event.get("agent_type"),
             "transcript_path": event.get("transcript_path"), "tool_name": event.get("tool_name"),
             "input_sha": hashlib.sha256(supplied.encode()).hexdigest(), "act": act[:200]}
    harness_mode.append_queue(harness_mode.state_dir(harness_mode.ROOT), entry, "prompts.jsonl")


def settle(event: dict) -> dict | None:
    """The output for a parsed event, or None to let it through."""
    sys.path.append(str(Path(__file__).resolve().parents[2] / "tools"))  # the stdlib wins (h16 G4)
    import harness_mode

    if event.get("hook_event_name") == "Stop":
        return stop_output(event, harness_mode)
    if event.get("hook_event_name") == "PostToolUse":
        return window_note(event, harness_mode)
    found = act_of(event)
    if found is None:
        return None
    act, why = found
    decision = harness_mode.settle(event, "guard_unattended", "ask", "", act, why)
    if decision.kind == "ask":  # attended: the normal question, prompt or reload
        if event["hook_event_name"] == "PermissionRequest":
            record_prompt(event, act, harness_mode)
        return None
    name = event["hook_event_name"]
    blocked = "unattended: configuration changes wait for Ola"
    return refusal(name, blocked if name == "ConfigChange" else decision.reason)


if __name__ == "__main__":
    sys.exit(main())
