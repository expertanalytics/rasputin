#!/usr/bin/env python3
"""PreToolUse: a rule needs an author, so a rule-bearing file needs a human.

Why this exists: on 2026-09-17 commit d1db597 -- whose stated purpose was
*deleting* false governance text -- added `Boost.Geometry` to CLAUDE.md section 2
and a bare `boost` key to tools/check_prohibited_deps.py, prohibiting every Boost
library in existence. Nobody had asked for it. The agent had read
project_structure.md's "Boost.Geometry, **if no longer used** after
vector_simplify is in-tree" -- a conditional note about dropping an unnecessary
dependency -- and hardened it into a machine-enforced ban. Six days later the
same agent cited that ban back as settled policy. The user: "I've prohibited
CGAL and GDAL."

The identical failure had already happened on the CRS rule, and PRINCIPLES.md
was in force both times. Prose did not stop it. This does.

WHAT IT GUARDS. Not "important files" -- files that STATE RULES, which is a much
smaller set. A rule in one of these binds every future session, so it is the one
kind of edit whose blast radius exceeds the branch it is made on.

WHY IT SCANS Bash TOO. An agent told to prefer shell tooling writes files with
heredocs and `python3 -`, not with Edit. A guard that matched only Edit and Write
would have watched d1db597 go past. The Bash arm is a substring scan and is
therefore both over- and under-inclusive: it cannot parse shell, so a path built
by variable expansion escapes it. It is a tripwire, not a sandbox.

DECISION IS `ask`, NEVER `deny`. The user may well want the edit. What must not
happen is the edit arriving inside an unrelated commit with nobody having said
yes. Answering the prompt IS the authorship record.
"""

import json
import re
import sys
from fnmatch import fnmatch

#: Files whose content is normative. A change here changes what is allowed.
GOVERNED = (
    "CLAUDE.md",
    "docs/PRINCIPLES.md",
    ".claude/REQUIRED-READING.md",
    "docs/increments/README.md",
)

#: Matched as a glob, not a literal: .claude/settings.local.json is the
#: permission allow-list and .claude/settings.json.pending-... is where the
#: settings currently sit. A literal ".claude/settings.json" guards neither.
GOVERNED_GLOBS = (".claude/settings*.json*",)

#: Directory prefixes where every file is a gate: these turn prose into refusals.
GOVERNED_PREFIXES = ("tools/check_", ".claude/agents/", ".claude/hooks/")

#: Shell constructs that write. Anything else containing a governed path is a
#: read -- `sed -n`, `grep`, `cat`, `git show` -- and must pass silently, or the
#: guard becomes noise and gets disabled, which is how guards die.
WRITES = re.compile(
    r"(>>?\s|<<\s*['\"]?\w*EOF|\bsed\s+-i|\btee\b|\bcp\b|\bmv\b|\bdd\b|\btruncate\b"
    r"|\bgit\s+(checkout|restore|apply|revert)\b|\bpatch\b|\bpython3?\s+-\b)",
)


def governed(path: str) -> bool:
    # removeprefix, not lstrip: lstrip("./") strips CHARACTERS, so it eats the
    # leading dot of ".claude/..." and every dotfile path stops matching.
    norm = path.replace("\\", "/").removeprefix("./")
    tail = norm.split("/")[-1]
    if any(norm.endswith(g) or tail == g for g in GOVERNED):
        return True
    if any(fnmatch(norm, f"*{pattern}") or fnmatch(tail, pattern.split("/")[-1])
           for pattern in GOVERNED_GLOBS):
        return True
    return any(prefix in norm for prefix in GOVERNED_PREFIXES)


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # A guard that crashes the session is worse than no guard.

    tool = event.get("tool_name", "")
    supplied = event.get("tool_input", {}) or {}
    hits: list[str] = []

    if tool in ("Edit", "Write", "NotebookEdit"):
        path = supplied.get("file_path") or supplied.get("notebook_path") or ""
        if governed(path):
            hits = [path]
    elif tool == "Bash":
        command = supplied.get("command", "")
        if WRITES.search(command):
            # Match the path as written, never its basename. A directory prefix
            # ends in "/", so its basename is "" -- and "" is a substring of
            # every string, which made the first draft of this hook fire on
            # every write-shaped command in the repository. A guard that fires
            # on everything trains the user to approve without reading, and an
            # approval nobody read is indistinguishable afterwards from a
            # considered one. That manufactures the exact artifact this hook
            # exists to prevent: a rule with consent attached and no author
            # behind it.
            hits = [
                name
                for name in (*GOVERNED, *GOVERNED_GLOBS, *GOVERNED_PREFIXES)
                if name.rstrip("/*") in command
            ]

    if not hits:
        return 0

    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "ask",
                    "permissionDecisionReason": (
                        f"This changes a file that states rules: {', '.join(sorted(set(hits)))}.\n"
                        "A rule here binds every future session. Before you say yes, ask what "
                        "the new or changed rule is, and who asked for it -- a commit or "
                        "something you said, not a nearby document that sounds similar.\n"
                        "This prompt is a tripwire, not a gate: it reads the command as text "
                        "and cannot parse shell, so a path built from a variable slips past it."
                    ),
                }
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
