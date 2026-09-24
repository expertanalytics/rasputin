#!/usr/bin/env python3
"""PreToolUse: the working tree is the agent's, the remote is the user's.

PRINCIPLES.md E1 states this and the user stated it directly -- "ensure one agent
working at the time, and alert me before pushing". It held because the agent
remembered to ask. That is not a mechanism; it is a habit, and the same session
has forgotten comparable ones.

`gh pr create`, `gh pr merge` and `git push` all publish. So does a `git commit`
carrying `--no-verify`, which is not a publish but is the disabling of somebody
else's guard, and belongs to the user for the same reason.

`ask`, not `deny`: the user says yes constantly. The point is that they say it.
"""

import json
import re
import sys

PUBLISHES = (
    (re.compile(r"\bgit\s+push\b"), "git push writes to the remote"),
    (re.compile(r"\bgh\s+pr\s+(create|merge|ready|edit)\b"), "gh pr changes a pull request"),
    (re.compile(r"\bgh\s+(release|repo\s+(create|delete|edit))\b"), "gh publishes or alters the repo"),
    (re.compile(r"--no-verify\b"), "--no-verify disables git's own hooks"),
    (re.compile(r"\bgit\s+(rebase|reset\s+--hard|filter-branch)\b|\bgit\s+commit\b.*--amend"),
     "this rewrites history, which is destructive once anything is published"),
    (re.compile(r"\bgit\s+push\b.*(--force|-f)\b"), "a force push can discard the user's commits"),
)


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    if event.get("tool_name") != "Bash":
        return 0
    command = (event.get("tool_input", {}) or {}).get("command", "")

    reasons = [why for pattern, why in PUBLISHES if pattern.search(command)]
    if not reasons:
        return 0

    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "ask",
                    "permissionDecisionReason": (
                        "This reaches beyond the working tree: "
                        + "; ".join(reasons)
                        + ".\nPRINCIPLES.md E1 -- the working tree is the agent's, the remote is "
                        "the user's. Approval of an earlier push does not carry to this one."
                    ),
                }
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
