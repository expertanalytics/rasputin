#!/usr/bin/env python3
"""PreToolUse on Agent and SendMessage: a persona spawn carries an intact, current brief.

Spec: docs/increments/h9-spawn-briefs.md §3.5 to §3.8. Always exits 0. A refusal
is a `deny` with its reason; a pass is silent, or carries `additionalContext`
when part of the check could not run. Without a working `tools/brief.py` the
spawn goes through with a notice (§3.7): refusing would also refuse the
@developer who repairs it. Unattended mode changes nothing, and nothing is queued.
"""

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
UNAVAILABLE = (
    "guard_spawn: brief check unavailable ({error}); this spawn was not checked. "
    "Tell Ola; repairing tools/brief.py is a @developer task."
)


def emit(**fields: str) -> int:
    print(json.dumps({"hookSpecificOutput": {"hookEventName": "PreToolUse", **fields}}))
    return 0


def refusal(brief, event: dict) -> str | None:
    """§3.5 rules 4 to 9, through brief.py's reader."""
    tool, given = event["tool_name"], event.get("tool_input") or {}
    text = given.get("prompt" if tool == "Agent" else "message") or ""
    try:
        blocks = brief.find_blocks(text)
    except brief.MalformedError as exc:
        return f"brief: malformed block ({exc}); paste the output of tools/brief.py unchanged"
    if len(blocks) > 1:
        return "brief: two blocks; a spawn carries one, from tools/brief.py"
    kind = given.get("subagent_type") if tool == "Agent" else None
    if not blocks:
        if kind in brief.PERSONAS:
            return (f"brief: no block for @{kind}; run: python3 tools/brief.py {kind} "
                    "--worktree <path> --beside <none|persona:path> [--increment <file>], "
                    "and paste its output before the task")  # fmt: skip
        return None
    block = blocks[0]
    reason = brief.check(block, ROOT)
    if reason and reason.startswith("brief: edited"):
        return reason
    if tool == "Agent" and kind != block.persona:
        return f"brief: persona: the block is for @{block.persona}, the spawn is @{kind}"
    if reason:
        return reason
    if tool == "Agent" and given.get("isolation"):
        return f"brief: isolation: @{kind} works in {block.worktree}, not in a new worktree"
    return None


def main() -> int:
    try:
        event = json.loads(sys.stdin.read())
    except ValueError:
        return 0
    if (not isinstance(event, dict) or event.get("tool_name") not in ("Agent", "SendMessage")
            or "agent_id" in event):  # fmt: skip
        return 0
    notes = []
    cwd, project = event.get("cwd"), os.environ.get("CLAUDE_PROJECT_DIR")
    if cwd and project:
        try:
            same = os.path.samefile(cwd, project)
        except OSError:
            same = False
        if not same:
            reason = f"cwd: this session is in {cwd}, not in {project}; run: cd {project}"
            return emit(permissionDecision="deny", permissionDecisionReason=reason)
    else:
        notes.append("guard_spawn: the working directory was not compared (no cwd in the "
                     "event, or CLAUDE_PROJECT_DIR unset)")  # fmt: skip
    try:
        sys.path.append(str(ROOT / "tools"))  # appended: the stdlib wins (h16 G4)
        import brief

        reason = refusal(brief, event)
    except Exception as exc:  # §3.7: unavailable, so unchecked, and said
        notes.append(UNAVAILABLE.format(error=f"{type(exc).__name__}: {exc}"))
        reason = None
    if reason:
        return emit(permissionDecision="deny", permissionDecisionReason=reason)
    if notes:
        return emit(additionalContext=" ".join(notes))
    return 0


if __name__ == "__main__":
    sys.exit(main())
