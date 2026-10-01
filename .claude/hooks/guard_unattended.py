#!/usr/bin/env python3
"""AskUserQuestion, PermissionRequest and ConfigChange: nobody is there to answer.

Spec: docs/increments/h3-unattended-u1.md §3.7. While Ola is at the keyboard
this hook prints nothing, except a refusal when it crashes or the
harness_mode module is missing. In unattended mode (on, or a flag that cannot be
read) a question, a permission prompt or a configuration reload would wait for
him and stall the run, so each is refused and queued instead.
"""

import json
import sys
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


def settle(event: dict) -> dict | None:
    """The output for a parsed event, or None to let it through."""
    found = act_of(event)
    if found is None:
        return None
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
    import harness_mode

    act, why = found
    decision = harness_mode.settle(event, "guard_unattended", "ask", "", act, why)
    if decision.kind == "ask":  # attended: the normal question, prompt or reload
        return None
    name = event["hook_event_name"]
    blocked = "unattended: configuration changes wait for Ola"
    return refusal(name, blocked if name == "ConfigChange" else decision.reason)


if __name__ == "__main__":
    sys.exit(main())
