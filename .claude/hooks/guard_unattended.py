#!/usr/bin/env python3
"""AskUserQuestion, PermissionRequest and ConfigChange: nobody is there to answer.

Spec: docs/increments/h3-unattended-u1.md §3.7. While Ola is at the keyboard
this hook prints nothing. In unattended mode (on, or a flag that cannot be
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


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0
    found = act_of(event)
    if found is None:
        return 0
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
    import harness_mode

    act, why = found
    decision = harness_mode.settle(event, "guard_unattended", "ask", "", act, why)
    if decision.kind == "ask":  # attended: the normal question, prompt or reload
        return 0
    name = event["hook_event_name"]
    if name == "PreToolUse":
        output = harness_mode.pretool("deny", decision.reason)
    elif name == "PermissionRequest":
        output = {"hookSpecificOutput": {"hookEventName": "PermissionRequest",
                                         "decision": {"behavior": "deny",
                                                      "message": decision.reason}}}
    else:
        output = {"decision": "block", "reason": "unattended: configuration changes wait for Ola"}
    print(json.dumps(output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
