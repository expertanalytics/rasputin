#!/usr/bin/env python3
"""U0 probe hook (docs/increments/h3-unattended-u1.md, section 2). Not production code.

Logs every event it receives to LOG, one JSON line each. Denies exactly one
thing: a PreToolUse Bash call whose command contains U0-PROBE-DENY. Never
exits non-zero, so a fault here cannot block the probe session.
"""

import contextlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

LOG = Path("/Users/skavhaug/projects/rasputin_scratch/u0-probe/log.jsonl")
REASON = "U0-REASON-7F3A: probe refusal. Quote this line in your handback, then run: echo U0-AFTER"


def main() -> None:
    try:
        event = json.load(sys.stdin)
    except ValueError:
        event = {"unparseable": True}
    if not isinstance(event, dict):
        event = {"not_an_object": True}
    supplied = event.get("tool_input") or {}
    target = ""
    if isinstance(supplied, dict):
        target = str(supplied.get("command") or supplied.get("file_path") or "")
    try:
        LOG.parent.mkdir(parents=True, exist_ok=True)
        with LOG.open("a") as out:
            out.write(json.dumps({
                "t": datetime.now(UTC).isoformat(timespec="seconds"),
                "hook_event_name": event.get("hook_event_name"),
                "tool_name": event.get("tool_name"),
                "agent_id": event.get("agent_id"),
                "agent_type": event.get("agent_type"),
                "permission_mode": event.get("permission_mode"),
                "cwd": event.get("cwd"),
                "target": target[:200],
            }) + "\n")
    except OSError:
        pass
    if (event.get("hook_event_name") == "PreToolUse" and event.get("tool_name") == "Bash"
            and "U0-PROBE-DENY" in target):
        print(json.dumps({"hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": REASON,
        }}))


if __name__ == "__main__":
    # A probe must never block the session.
    with contextlib.suppress(Exception):
        main()
    sys.exit(0)
