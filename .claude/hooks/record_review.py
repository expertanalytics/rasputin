#!/usr/bin/env python3
"""PreToolUse (SubagentHandback) and SubagentStop: record `@reviewer`'s verdict.

Spec: docs/increments/h20-away-push-grants.md §3.3. `@reviewer` ends its
handback with one line `Review: <verdict> <40-hex sha>`; this hook appends it
to `<git-common-dir>/harness/reviews.jsonl`, which the away-mode push grant
(tools/push_grant.py) reads. Only an event whose `agent_type` is `reviewer`
counts; Claude Code sets that field, no agent can. It records in every mode,
never blocks, prints nothing and exits 0, also on an error.
"""

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

# Appended, not put first: a tools/ file named after a stdlib module must not
# replace it here (h16 §2 G4).
sys.path.append(str(Path(__file__).resolve().parents[2] / "tools"))


def main() -> int:
    try:
        import harness_mode
        import push_grant

        event = json.load(sys.stdin)
        if event.get("agent_type") != "reviewer":
            return 0
        if event.get("hook_event_name") == "SubagentStop":
            text = event.get("last_assistant_message")
        else:
            text = (event.get("tool_input") or {}).get("message")
        found = push_grant.parse_review(text) if isinstance(text, str) else None
        if found is None:
            return 0
        line = {"at": datetime.now(UTC).isoformat(timespec="seconds"), "sha": found[1],
                "verdict": found[0], "agent_id": event.get("agent_id"),
                "session_id": event.get("session_id"),
                "event": event.get("hook_event_name")}  # fmt: skip
        state = harness_mode.state_dir(harness_mode.ROOT)
        harness_mode.append_queue(state, line, "reviews.jsonl")
    except Exception:  # never block a handback; nothing is recorded
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
