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
would have watched d1db597 go past. The Bash arm judges the files a line writes
(tools/shell_scan.py; docs/increments/h4-guard-fixes.md §3), not the words in it,
so a governed name in a heredoc body, a quotation or a read passes. A line it
cannot read -- unparseable, a target that is wholly a variable, a writer whose
targets are unnamed -- is judged by the text rule it had before h4. It is a
tripwire, not a sandbox.

DECISION IS `ask`, NEVER `deny`. The user may well want the edit. What must not
happen is the edit arriving inside an unrelated commit with nobody having said
yes. Answering the prompt IS the authorship record. The exceptions: in
unattended mode (tools/harness_mode.py) the ask becomes a `deny` and is queued
for Ola, and a write of the harness state or a run of away.py is always denied
(docs/increments/h3-unattended-u1.md §3.6).
"""

import json
import re
import sys
from fnmatch import fnmatch
from pathlib import Path

# Appended, not put first: a tools/ file named after a stdlib module must not
# replace it here (docs/increments/h16-harness-fixes.md §2 G4).
sys.path.append(str(Path(__file__).resolve().parents[2] / "tools"))
try:
    import shell_scan
except ImportError:  # every line is then judged as text, as before h4
    shell_scan = None
try:
    import scratchpad
except ImportError:  # no exemption then: a scratchpad path is judged like any other
    scratchpad = None

#: Files whose content is normative. A change here changes what is allowed.
GOVERNED = (
    "CLAUDE.md",
    "docs/PRINCIPLES.md",
    ".claude/REQUIRED-READING.md",
    "docs/increments/README.md",
    # The self-protecting set: live before any review (unattended mode, §3.6).
    "tools/away.py",
    "tools/harness_mode.py",
    "tools/session_state.py",
    "tools/rule_sizes.py",
    "tools/shell_scan.py",
    "tools/brief.py",  # h9: guard_spawn.py imports it
    "tools/scratchpad.py",  # h16: this guard imports it
    "tools/count_loc.py",  # h16: it computes the arithmetic of CLAUDE.md §2
    ".claude/profile.toml",
)

#: A glob, not a literal: settings.json wires the hooks, settings.local.json is
#: the permission allow-list, and a settings.json.pending-... draft is guarded
#: too. A literal ".claude/settings.json" would miss the other two.
GOVERNED_GLOBS = (".claude/settings*.json*",)

#: Directory prefixes where every file is a gate: these turn prose into refusals.
GOVERNED_PREFIXES = ("tools/check_", ".claude/agents/", ".claude/hooks/", ".claude/skills/",
                     ".git/hooks/", ".claude/briefs/", ".git/remotes/", ".git/branches/")

#: Matched with endswith only: a basename rule would govern every file named config.
GOVERNED_SUFFIXES = (".git/config",)

#: The text rule, for a line the parser cannot read: shell constructs that write.
#: Anything else containing a governed path is a read -- `sed -n`, `grep`, `cat`,
#: `git show` -- and must pass silently, or the guard becomes noise and gets
#: disabled, which is how guards die.
WRITES = re.compile(
    r"(>>?\s|<<\s*['\"]?\w*EOF|\bsed\s+-i|\btee\b|\bcp\b|\bmv\b|\bdd\b|\btruncate\b"
    r"|\bgit\s+(checkout|restore|apply|revert)\b|\bpatch\b"
    r"|\bpython3?\s+-(?!m\s+(?:ruff|mypy|pytest)\b)\b)",  # python -m <gate> is a read
)


def governed(path: str, *, scratch_exempt: bool = True) -> bool:
    # h16 G3a, before every rule below: a scratchpad is temporary, and nothing
    # the harness reads lives there. By real path, so a symlink out is followed.
    if scratch_exempt and scratchpad is not None and scratchpad.under(path):
        return False
    # removeprefix, not lstrip: lstrip("./") strips CHARACTERS, so it eats the
    # leading dot of ".claude/..." and every dotfile path stops matching.
    norm = path.replace("\\", "/").removeprefix("./")
    parts = norm.split("/")
    tail = parts[-1]
    if any(norm.endswith(g) or tail == g for g in GOVERNED):
        return True
    if any(fnmatch(norm, f"*{pattern}") or fnmatch(tail, pattern.split("/")[-1])
           for pattern in GOVERNED_GLOBS):
        return True
    if any(norm.endswith(suffix) for suffix in GOVERNED_SUFFIXES):
        return True
    # h16 G4: `python3 tools/x.py` puts tools/ first on sys.path, so a tools/ file
    # named after a stdlib module (tools/json/__init__.py, tools/ast.py) replaces it.
    if any(parts[at - 1] == "tools" and part.split(".")[0] in sys.stdlib_module_names
           for at, part in enumerate(parts) if at):
        return True
    return any(prefix in norm for prefix in GOVERNED_PREFIXES)


#: Only Ola enters or leaves unattended mode; only hooks and away.py write its state.
HARNESS_PATH = re.compile(r"(^|/)\.git/harness(/|$)")
STATE_WRITES = re.compile(r"\b(rm|touch|ln|mkdir|unlink|install)\b")
#: A run of the away script in a line the parser cannot read: its name as a word.
AWAY_TEXT = re.compile(r"(^|[\s/'\"])away\.py\b|-m\s+(\S+\.)?away\b")
ALWAYS_DENIED = (
    "Only Ola enters or leaves unattended mode, and only hooks and away.py write the "
    "harness state. Nothing is queued: this act is not an agent's to wait for."
)
CANNOT_READ = "the guard cannot read this command's targets"


def runs_away(simple) -> bool:
    """argv[0], or the script or `-m` module of an interpreter, is away.py."""
    first = simple.argv[0] if simple.argv else ""
    return any(shell_scan.base(word or "") == "away.py" for word in (first, simple.script))


def judge_bash(command: str) -> tuple[str, list[str], str] | None:
    """("deny", [], "") or ("ask", governed targets, why) for a Bash line; None passes it."""
    simples = shell_scan.parse(command) if shell_scan else None
    if simples is not None:
        targets = [t for s in simples for t in (*s.writes, *shell_scan.candidates(s))]
        # An expansion counts as matching .git: $(git rev-parse --git-common-dir)/harness.
        # Backstop, as before h4: a line that writes anything and names the state
        # dir is denied, however the path reaches the write (export, read, cd).
        if (any(HARNESS_PATH.search(shell_scan.static(t, ".git")) for t in targets)
                or any(runs_away(s) for s in simples)
                or (".git/harness" in command and any(s.writes or s.unknown for s in simples))):
            return "deny", [], ""
        # A target is judged by its static tail: $D/CLAUDE.md is CLAUDE.md. The scratchpad
        # holds only for one plain command, which cannot link a path and write through it.
        one = simples[0] if len(simples) == 1 else None
        plain = one is not None and one.program is None and one.script is None and not one.unknown
        hits = list(dict.fromkeys(t for t in targets
                                  if governed(shell_scan.static(t), scratch_exempt=plain)))
        if hits:
            return "ask", hits, f"it writes {', '.join(hits)}"
        if not any(s.unknown or any(not shell_scan.static(w).strip("/") for w in s.writes)
                   for s in simples):
            return None
    elif AWAY_TEXT.search(command) or (
            ".git/harness" in command and (WRITES.search(command) or STATE_WRITES.search(command))):
        return "deny", [], ""
    hits = text_hits(command) if WRITES.search(command) else []
    return ("ask", hits, CANNOT_READ) if hits else None


def text_hits(command: str) -> list[str]:
    """The governed names in a line's text: the rule before h4, for a line it cannot read."""
    # Match the path as written, never its basename. A directory prefix ends in
    # "/", so its basename is "" -- and "" is a substring of every string, which
    # made the first draft of this hook fire on every write-shaped command in the
    # repository. A guard that fires on everything trains the user to approve
    # without reading, and an approval nobody read is indistinguishable afterwards
    # from a considered one. That manufactures the exact artifact this hook exists
    # to prevent: a rule with consent attached and no author behind it.
    hits = [name for name in (*GOVERNED, *GOVERNED_PREFIXES, *GOVERNED_SUFFIXES)
            if name.rstrip("/") in command]
    # A glob cannot be used as a substring: ".claude/settings*.json*" contains a
    # literal "*" that no real command does, so the settings files -- including
    # settings.local.json, the permission allow-list -- went unguarded on this arm
    # while the Edit arm caught them. Match the globs against the command's
    # whitespace-split tokens.
    hits += [pattern for pattern in GOVERNED_GLOBS
             if any(fnmatch(token.strip("\"'<>"), f"*{pattern}") for token in command.split())]
    return sorted(set(hits))


def emit(event: dict, verdict: str, reason: str, act: str, why: str) -> None:
    """Route the verdict through harness_mode; unable to import it, deny anyway."""
    try:
        import harness_mode
    except ImportError as error:
        failed = f"guard_governance cannot import harness_mode ({error}); refused: {why}."
        output = deny(failed if verdict == "ask" else reason)
    else:
        output = harness_mode.guard(event, "guard_governance", verdict, reason, act, why)
    if output is not None:
        print(json.dumps(output))


def deny(reason: str) -> dict:
    return {"hookSpecificOutput": {"hookEventName": "PreToolUse",
                                   "permissionDecision": "deny",
                                   "permissionDecisionReason": reason}}


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # A guard that crashes the session is worse than no guard.
    try:
        verdict(event)
    except Exception as error:  # after parsing, a crash must not turn an ask into a pass
        print(json.dumps(deny(f"guard_governance failed: {type(error).__name__}: {error}")))
    return 0


def verdict(event: dict) -> None:
    tool = event.get("tool_name", "")
    supplied = event.get("tool_input", {}) or {}
    hits: list[str] = []
    act = why = ""

    if tool in ("Edit", "Write", "NotebookEdit"):
        path = supplied.get("file_path") or supplied.get("notebook_path") or ""
        if HARNESS_PATH.search(path.replace("\\", "/")):
            emit(event, "deny", ALWAYS_DENIED, "", "")
            return
        if governed(path):
            hits, act, why = [path], f"{tool} {path}", f"it writes {path}"
    elif tool == "Bash":
        act = supplied.get("command", "")
        judged = judge_bash(act)
        if judged is not None and judged[0] == "deny":
            emit(event, "deny", ALWAYS_DENIED, "", "")
            return
        if judged is not None:
            _, hits, why = judged

    if not hits:
        return

    emit(event, "ask", (
        f"This changes a file that states rules: {', '.join(hits)}.\n"
        "A rule here binds every future session. Before you say yes, ask what "
        "the new or changed rule is, and who asked for it -- a commit or "
        "something you said, not a nearby document that sounds similar.\n"
        "This prompt is a tripwire, not a gate: it judges the files a shell line "
        "writes, and a line it cannot read by its text."
    ), act, why)


if __name__ == "__main__":
    sys.exit(main())
