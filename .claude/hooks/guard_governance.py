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
import os
import posixpath
import re
import sys
from fnmatch import fnmatch
from itertools import pairwise
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
    "tools/push_grant.py",  # h20: guard_push and record_review import it
    "tools/ci_changes.py",  # h20: push_grant imports it
    "tools/work_pipes.py",  # h31: session_state imports it
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
    r"(>>?&?[|!]?\s|<<\s*['\"]?\w*EOF|\bsed\s+-i|\btee\b|\bcp\b|\bmv\b|\bdd\b|\btruncate\b"
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
    # h16 G4: judged also with `.`, `..` and `//` resolved (lexically), as well as
    # as written: normpath drops a trailing `/`, which the prefixes below need.
    lexical = posixpath.normpath(norm)
    if lexical not in (norm, ".") and governed(lexical, scratch_exempt=scratch_exempt):
        return True
    # Names compare case-folded, as the macOS disk does; the module rule below stays
    # exact, since Python's import checks letter case there.
    parts, low = norm.split("/"), norm.lower()
    tail = low.split("/")[-1]
    if any(low.endswith(g.lower()) or tail == g.lower() for g in GOVERNED):
        return True
    if any(fnmatch(low, f"*{pattern}".lower()) or fnmatch(tail, pattern.split("/")[-1].lower())
           for pattern in GOVERNED_GLOBS):
        return True
    if any(low.endswith(suffix) for suffix in GOVERNED_SUFFIXES):
        return True
    # h16 G4: `python3 tools/x.py` puts tools/ first on sys.path, so a tools/ file
    # named after a stdlib module (tools/json/__init__.py, tools/ast.py) replaces it.
    if any(parts[at - 1] == "tools" and part.split(".")[0] in sys.stdlib_module_names
           for at, part in enumerate(parts) if at):
        return True
    return any(prefix.lower() in low for prefix in GOVERNED_PREFIXES)


GIT_DIR = re.compile(r"(^|/)\.git(/|$)", re.I)
#: A `NAME=value` the guard reads its words by: no quote, expansion or separator in it.
KNOWN = re.compile(r"[^\s'\"$`;&|()<>\\]*")
#: A link source that is a git folder by what it runs.
GIT_DIR_RUN = re.compile(r"rev-parse\b.*--(git-dir|git-common-dir|absolute-git-dir)\b")
#: On a line naming the record folder, the programs that pass: none writes a file or
#: runs a program, whatever its options (rg, git and python3 are judged further).
READERS = {"echo", "printf", "cat", "grep", "head", "tail", "wc", "ls", "stat", "diff", "cut"}
#: `${` other than `${NAME}`, `$((` and `$[` run code the reader does not see (zsh's `(e)`).
HIDDEN_RUN = re.compile(r"\$\{(?!\w+\})|\$\(\(|\$\[")
SHORT = re.compile(r"-[A-Za-z]+")
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


def in_folder(path: str, globs: bool = True) -> bool:
    """The path names the record folder: two adjacent parts stand for `.git` and `harness`.

    Case-folded, as the macOS disk is. A part stands for its name when it is that name,
    an expansion or, with `globs`, a pattern matching it; one of the two is written out,
    so `$X/*` and `*/*` never match.
    """
    text = shell_scan.static(path, "\0") if shell_scan else path
    parts = posixpath.normpath(text.replace("\\", "/")).lower().split("/")

    def stands(part: str, name: str) -> bool:
        pattern = re.sub(r"\{[^}]*\}", "*", part)
        return part == name or (part != "" and not part.strip("\0")) or (
            globs and fnmatch(name, pattern))  # fmt: skip

    return any((a == ".git" or b == "harness") and stands(a, ".git") and stands(b, "harness")
               for a, b in pairwise(parts))  # fmt: skip


def real(path: str, cwd: str) -> str | None:
    """The path's real path from `cwd`, a link followed; None for a path holding an expansion."""
    return None if shell_scan.substituted(path) else os.path.realpath(os.path.join(cwd, path))


def judge_bash(command: str, *, cwd: str | None = None) -> tuple[str, list[str], str] | None:
    """("deny", [], "") or ("ask", governed targets, why) for a Bash line; None passes it."""
    assigned: dict[str, str] = {}
    simples = shell_scan.parse(command, assigned=assigned) if shell_scan else None
    if simples is not None:
        known = {name: value for name, value in assigned.items() if KNOWN.fullmatch(value)}
        cwd, here, walked = cwd or os.getcwd(), "", []
        for simple in simples:  # each target joined to the cd path so far, known values in
            for target in (*simple.writes, *shell_scan.candidates(simple)):
                path = posixpath.join(here, value_of(target, known))
                walked.append((target, path, real(path, cwd)))
            if simple.written[:1] in (["cd"], ["pushd"]) and simple.written[1:]:
                here = posixpath.join(here, value_of(simple.written[1], known))
        # An expansion counts as matching .git: $(git rev-parse --git-common-dir)/harness.
        # Backstop, as before h4: a line that writes anything and names the state
        # dir is denied, however the path reaches the write (export, read, cd).
        if (any(in_folder(path) or (found is not None and in_folder(found))
                for _, path, found in walked)
                or any(runs_away(s) for s in simples) or reaches_folder(simples, known)
                or (names_folder(command, simples) and not only_reads(command, simples))):
            return "deny", [], ""
        # A target is judged by its static tail: $D/CLAUDE.md is CLAUDE.md, and by its real
        # path. The scratchpad holds only for one plain command, which cannot link a path
        # and write through it.
        one = simples[0] if len(simples) == 1 else None
        plain = one is not None and one.program is None and one.script is None and not one.unknown
        hits = list(dict.fromkeys(
            t for t, _, found in walked if governed(shell_scan.static(t), scratch_exempt=plain)
            or (found is not None and governed(found, scratch_exempt=plain))))  # fmt: skip
        if hits:
            return "ask", hits, f"it writes {', '.join(hits)}"
        if not any(s.unknown or any(not shell_scan.static(w).strip("/") for w in s.writes)
                   for s in simples):
            return None
    elif AWAY_TEXT.search(command) or names_folder(command, None):
        return "deny", [], ""
    unread = simples is not None and any(s.unknown for s in simples)
    hits = text_hits(command) if WRITES.search(command) or unread else []
    return ("ask", hits, CANNOT_READ) if hits else None


def value_of(word: str, known: dict[str, str]) -> str:
    """The word with each `$NAME` or `${NAME}` of a known assignment replaced by its value."""
    return shell_scan.PLAIN_PARAMETER.sub(
        lambda match: known.get(match.group(1) or match.group(2), match.group(0)), word)


def reaches_folder(simples, known: dict[str, str]) -> bool:
    """A command with a word in the folder that is not a plain reader, or an `ln` of a git
    folder outside a scratchpad: an earlier line's link would reach the folder."""
    moved = any(s.written[:1] in (["cd"], ["pushd"]) for s in simples)
    for simple in simples:
        words = [value_of(word, known) for word in simple.written]
        reads = allowed(simple, moved) and all(plain(t) and not in_folder(t) for t in simple.writes)
        if not reads and any(in_folder(word, globs=False) for word in words):
            return True
        if shell_scan.base(simple.argv[0] if simple.argv else "") != "ln":
            continue
        found = shell_scan.operands(simple.argv[1:], shell_scan.COPIERS["ln"])
        for source in (value_of(word, known) for word in found[:-1] or found):
            if scratchpad is not None and scratchpad.under(source):
                continue
            last = posixpath.normpath(shell_scan.static(source)).split("/")[-1].lower()
            if last == ".git" or GIT_DIR_RUN.search(source) or in_folder(source):
                return True
    return False


def names_folder(command: str, simples) -> bool:
    """The line names the record folder (any case, quotes removed) or changes into `.git`."""
    text = command.lower()
    if ".git/harness" in text or ".git/harness" in re.sub(r"['\"\\]", "", text):
        return True
    if simples is None:  # unread: any cd or pushd beside a `.git` path component
        return bool(re.search(r"\b(cd|pushd)\b", text)) and bool(re.search(
            r"(^|[\s/'\"=])\.git(/|\s|$)", text))  # fmt: skip
    return any(s.written[:1] in (["cd"], ["pushd"]) and any(map(GIT_DIR.search, s.written[1:]))
               for s in simples)  # fmt: skip


def plain(path: str) -> bool:
    """A path the shell writes as written: no expansion, glob, `!` or `|`, `~` only as `~/`."""
    return (bool(shell_scan.PLAIN_PATH.fullmatch(path)) and not set("!|") & set(path)
            and (not path.startswith("~") or path.startswith("~/")))  # fmt: skip


def allowed(simple, moved: bool) -> bool:
    """A command on the short list that cannot write a file or run a program."""
    words = simple.written
    if simple.unknown or simple.assigns or simple.runs_body or words != simple.argv:
        return False
    name, rest = (words[0], words[1:]) if words else ("", [])
    if not words or (len(words) == 1 and name in shell_scan.CLOSERS) or name in READERS:
        return True
    if name == "rg":
        return not any(w.startswith(("--pre", "--hostname-bin")) for w in rest)
    if name in ("cd", "pushd"):
        return len(rest) < 2 and all(plain(w) and not w.startswith("-") for w in rest)
    if name in ("python3", "python"):
        return rest[:1] == ["tools/brief.py"] and not moved
    if name != "git" or not rest:
        return False
    if rest[0] in ("log", "show", "diff", "status"):
        return not any(w.startswith("--output") for w in rest)
    shorts = [w for w in rest if SHORT.fullmatch(w)]
    message = any(w.startswith(("--message", "--file")) for w in rest) or any(
        set("mF") & set(w) for w in shorts)  # fmt: skip
    edits = any(w == "--edit" for w in rest) or any("e" in w for w in shorts)
    return rest[0] == "commit" and message and not edits


def only_reads(command: str, simples) -> bool:
    """A line naming the record folder: every command allowed, every write plain and outside
    it, no `cd` into it, and no expansion that runs code the reader did not parse."""
    if HIDDEN_RUN.search(command):
        return False
    moved = any(s.written[:1] in (["cd"], ["pushd"]) for s in simples)
    here = ""
    for simple in simples:
        if not allowed(simple, moved):
            return False
        for target in simple.writes:
            if not plain(target) or in_folder(posixpath.join(here, target)):
                return False
        if simple.written[:1] in (["cd"], ["pushd"]) and simple.written[1:]:
            here = posixpath.normpath(posixpath.join(here, simple.written[1]))
            if in_folder(here):
                return False
    return True


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
        found = os.path.realpath(os.path.join(event.get("cwd") or os.getcwd(), path))
        if in_folder(path) or in_folder(found):
            emit(event, "deny", ALWAYS_DENIED, "", "")
            return
        if governed(path) or governed(found):
            hits, act, why = [path], f"{tool} {path}", f"it writes {path}"
    elif tool == "Bash":
        act = supplied.get("command", "")
        judged = judge_bash(act, cwd=event.get("cwd"))
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
