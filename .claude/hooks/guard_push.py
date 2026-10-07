#!/usr/bin/env python3
"""PreToolUse: the working tree is the agent's, the remote is the user's.

PRINCIPLES.md E1 states this and the user stated it directly -- "ensure one agent
working at the time, and alert me before pushing". It held because the agent
remembered to ask. That is not a mechanism; it is a habit, and the same session
has forgotten comparable ones.

`gh pr create`, `gh pr merge` and `git push` all publish, and `gh pr update-branch`
writes to a pull request's branch (docs/increments/h10-merge-queue.md §4a). So
does a `git commit` carrying `--no-verify`, which is not a publish but is the
disabling of somebody else's guard, and belongs to the user for the same reason.
So do writes to refs, remotes and git config, and `gh api` / `curl` with a
writing method to the forge (docs/increments/h3-unattended-u1.md §3.5).

`ask`, not `deny`: the user says yes constantly. The point is that they say it.
The exception is unattended mode (tools/harness_mode.py): while Ola is away the
ask becomes a `deny`, and the act is queued for him.

The rules apply to each simple command's argv (tools/shell_scan.py;
docs/increments/h4-guard-fixes.md §3), never to text, so `git push` in a heredoc
body, a commit message or a PR body is not a push. A line the parser cannot
read is judged by the text rules below, as before h4.

docs/increments/h16-harness-fixes.md §2 adds: a fetch into a named ref and a
`git replace` write ask (G1); a git or gh command the guard does not know, such
as an alias, asks, since it cannot see what runs (G2); and gh's group and verb
are read past its `-R`/`--repo` option (G6), and glued or clustered `gh api`
and `curl` options are read as the flags they are (G6, round 7). A git or gh
command run by another program is judged too (G7, `runs`).
"""

import functools
import json
import re
import shlex
import subprocess
import sys
from itertools import pairwise
from pathlib import Path

# Appended, not put first: a tools/ file named after a stdlib module must not
# replace it here (h16 §2 G4).
sys.path.append(str(Path(__file__).resolve().parents[2] / "tools"))
try:
    import shell_scan
except ImportError:  # every line is then judged as text, as before h4
    shell_scan = None

PUBLISHES = (
    (re.compile(r"\bgit\b[^|;&]*\bpush\b"), "git push writes to the remote"),
    (re.compile(r"\bgh\b[^|;&]*\bpr\b[^|;&]*\b(create|new|merge|ready|edit|update-branch)\b"),
     "gh pr changes a pull request"),
    (re.compile(r"\bgh\b[^|;&]*\b(release|repo\b[^|;&]*\b(create|new|delete|edit))\b"),
     "gh publishes or alters the repo"),
    (re.compile(r"--no-verify\b"), "--no-verify disables git's own hooks"),
    (re.compile(r"\bgit\s+(rebase|reset\s+--hard|filter-branch)\b|\bgit\s+commit\b.*--amend"),
     "this rewrites history, which is destructive once anything is published"),
    (re.compile(r"\bgit\b[^|;&]*\bpush\b.*(--force|-f)\b"),
     "a force push can discard the user's commits"),
)
PUSH, PR, RELEASE, NO_VERIFY, HISTORY, FORCE = (why for _, why in PUBLISHES)
UNKNOWN = "a git or gh alias, or a command the guard does not know; it cannot see what it runs"

FORGE_HOST = "github.com"
SEGMENTS = re.compile(r"&&|\|\||;|\||\n")
#: Options whose next token is their argument, not a positional (§3.5).
TAKES_ARG = {"-f", "--file", "--blob", "--type", "--default", "--comment", "-m"}
#: git's own options, before the subcommand, whose value is the next word.
GIT_TAKES_ARG = {"-C", "-c", "--git-dir", "--work-tree", "--namespace", "--attr-source",
                 "--config-env"}  # fmt: skip
FETCH = "fetch writes a named ref"
REMOTE_WRITES = {"add", "set-url", "rename", "remove", "rm", "set-head", "set-branches"}
CONFIG_WRITE_FLAGS = {"--add", "--unset", "--unset-all", "--replace-all", "--rename-section",
                      "--remove-section", "--edit", "-e"}
CONFIG_WRITE_VERBS = {"set", "unset", "rename-section", "remove-section", "edit"}
CONFIG_READS = {"--list", "-l", "get", "list"}
REPLACE_LISTS = {"-l", "--list"}
GH_FIELDS = {"-f", "-F", "--field", "--raw-field", "--input"}
CURL_DATA = ("--data", "--form", "--json", "--upload-file")
#: Config keys that steer a fetch (G1): prefixes, and exact keys, compared lower-cased.
FETCH_KEYS = ("remote.", "url.", "include.", "includeif.")
FETCH_EXACT = ("core.sshcommand", "fetch.bundleuri")
#: Assignments stripped from the argv that steer a fetch, read from the text instead.
FETCH_TEXT = re.compile(r"GIT_CONFIG|GIT_SSH|\b(HOME|XDG_CONFIG_HOME)=")
#: Programs that take a command as one string (G7 c), and shells (G7 b).
STRING_RUNNERS, SHELLS = {"watch", "parallel", "flock"}, shell_scan.SHELLS if shell_scan else set()
#: gh's top-level commands (`gh help`, gh 2.101) less its alias `co`, which a user can redefine,
#: plus `help` itself, which only reads (Ola's ruling, 2026-10-05).
GH_COMMANDS = {
    "help", "auth", "browse", "codespace", "discussion", "gist", "issue", "org", "pr", "project",
    "release", "repo", "skill", "cache", "run", "workflow", "agent-task", "alias", "api",
    "attestation", "completion", "config", "copilot", "extension", "gpg-key", "label",
    "licenses", "preview", "ruleset", "search", "secret", "ssh-key", "status", "variable",
}  # fmt: skip


def tokens(segment: str) -> list[str]:
    try:
        return shlex.split(segment)
    except ValueError:
        return segment.split()


def positionals(args: list[str]) -> list[str]:
    found, skip = [], False
    for token in args:
        if skip:
            skip = False
        elif token in TAKES_ARG:
            skip = True
        elif not token.startswith("-"):
            found.append(token)
    return found


def git_call(words: list[str]) -> tuple[str, list[str]] | None:
    """(subcommand, its arguments) of a command that is `git`, skipping git's own options."""
    if not words or words[0].rsplit("/", 1)[-1] != "git":
        return None
    at = 1
    while at < len(words) and words[at].startswith("-"):
        at += 2 if words[at] in GIT_TAKES_ARG else 1
    return (words[at], words[at + 1:]) if at < len(words) else None


@functools.cache
def git_commands() -> frozenset[str]:
    """The current git commands, which no alias can hide: main less deprecated (G2)."""
    listed = [subprocess.run(["git", f"--list-cmds={kind}"], capture_output=True, text=True,
                             check=True, timeout=10).stdout.split()
              for kind in ("main", "deprecated")]  # fmt: skip
    return frozenset(listed[0]) - frozenset(listed[1])


def method(args: list[str], short: str, long: str) -> str | None:
    """The HTTP method given as `-X M`, `-XM`, `--long M` or `--long=M`, if any."""
    for at, token in enumerate(args):
        if token in (short, long) and at + 1 < len(args):
            return args[at + 1].upper()
        if token.startswith(short) and len(token) > len(short):
            return token[len(short):].upper()
        if token.startswith(long + "="):
            return token.split("=", 1)[1].upper()
    return None


def git_why(sub: str, args: list[str]) -> str | None:
    """The reason a git subcommand writes a ref, a remote or git config, or None."""
    pos = positionals(args)
    if sub not in git_commands():
        return UNKNOWN
    if sub == "update-ref":
        return "update-ref moves a ref directly"
    if sub == "remote" and pos and pos[0] in REMOTE_WRITES:
        return "this changes where the remote points"
    if sub == "config":
        reads = any(a.startswith("--get") or a in CONFIG_READS for a in args)
        if (CONFIG_WRITE_FLAGS & set(args) or (pos and pos[0] in CONFIG_WRITE_VERBS)
                or (len(pos) >= 2 and not reads)):
            return "this writes git configuration (hooks path, remote URLs)"
    if sub == "symbolic-ref" and (len(pos) >= 2 or {"-d", "--delete"} & set(args)):
        return "symbolic-ref rewrites a symbolic ref"
    # A refspec with a destination after the repository; every non-option word
    # counts, so an option's value can only push the refspec later, never hide it.
    if sub in ("fetch", "pull") and any(":" in a for a in [a for a in args if a[:1] != "-"][1:]):
        return FETCH
    if sub == "replace":
        flags = {a for a in args if a.startswith("-")}
        listing = all(f in REPLACE_LISTS or f.startswith("--format=") for f in flags)
        if not listing or (pos and not flags & REPLACE_LISTS):
            return "replace refs change what git reads for an object"
    return None


def fetch_named(words: list[str], sub: str, args: list[str], text: str) -> bool:
    """G1: a fetch steered into named refs by --refmap, --stdin, a config override or the env."""
    if not (sub in ("fetch", "pull") or (sub == "remote" and positionals(args)[:1] == ["update"])):
        return False
    keys = [a.split("=")[0] for a in args]  # git takes any unambiguous prefix of a long option
    if any((len(k) >= 5 and "--refmap".startswith(k)) or (len(k) >= 4 and "--stdin".startswith(k))
           for k in keys):  # fmt: skip
        return True
    head = words[1:len(words) - len(args) - 1]
    given = [v for f, v in pairwise(head) if f in ("-c", "--config-env")]
    given += [w.split("=", 1)[1] for w in head if w.startswith("--config-env=")]
    keys = [g.split("=", 1)[0].lower() for g in given]
    steered = any(k.startswith(FETCH_KEYS) or k in FETCH_EXACT for k in keys)
    return steered or FETCH_TEXT.search(text) is not None


def gh_words(words: list[str]) -> list[str] | None:
    """A gh command's words past its options, `-R`/`--repo`'s value included (G6); else None."""
    if not words or words[0].rsplit("/", 1)[-1] != "gh":
        return None
    found, skip = [], False
    for word in words[1:]:
        if not skip and not word.startswith("-"):
            found.append(word)
        skip = not skip and word in ("-R", "--repo")
    return found


def segment_why(words: list[str], text: str = "") -> str | None:
    """The reason a segment writes a ref, a remote, git config or the forge, or None."""
    call = git_call(words)
    if call is not None:
        why = git_why(*call)
        return FETCH if why is None and fetch_named(words, *call, text) else why
    gh = gh_words(words)
    if gh and gh[0] not in GH_COMMANDS:
        return UNKNOWN
    if gh and gh[0] == "api":
        # A run of `-i` (gh's one short flag without a value) glued before another flag is cut.
        rest = [re.sub(r"^-i+(?=[^-i])", "-", a) for a in words[words.index("api") + 1:]]
        given = method(rest, "-X", "--method")
        fields = any(a.split("=")[0] in GH_FIELDS or a[:2] in ("-f", "-F") for a in rest)
        if (given is not None and given != "GET") or (given is None and fields):
            return "gh api with a writing method changes the forge"
    if "curl" in words and any(FORGE_HOST in w for w in words):
        # `--expand-x` is `--x`; a short cluster has data letters before its `X`, its method after.
        longs = ["--" + w[9:] if w.startswith("--expand-") else w for w in words]
        clusters = [w for w in words if w[:1] == "-" and w[1:2] not in ("", "-")]
        given = method([f"-X{w.split('X', 1)[1]}" if w in clusters and "X" in w else w
                        for w in longs], "-X", "--request")  # fmt: skip
        data = any(w.startswith(CURL_DATA) for w in longs) or any(
            set(w.split("X")[0]) & set("dFT") for w in clusters)
        if (given is not None and given not in ("GET", "HEAD")) or data:
            return "curl with a writing method to the forge"
    return None


def runs(words: list[str], under=False, quoted=False) -> list[tuple[list[str], bool]]:
    """G7: (argv, bare) for the command and each command its later words run.

    under: this command or a runner before it is parallel; quoted: one is in
    STRING_RUNNERS. A runner or shell word hands the rest of the line to one
    call and the scan stops. Not linear: each shell word parses the rest of
    the line again, so time grows with shell words times line length. A line
    past Python's recursion limit (about 990 shell words) raises
    RecursionError, which main's except turns into a deny (fails closed).
    """
    found, name = [(words, False)], words[0].rsplit("/", 1)[-1] if words else ""
    under, quoted = under or name == "parallel", quoted or name in STRING_RUNNERS
    for at, word in enumerate(words[1:], 1):
        if (tail := word.rsplit("/", 1)[-1]) in ("git", "gh"):  # (a) a bare tail
            found.append((words[at:], not under))  # parallel builds it from inputs
        if tail in STRING_RUNNERS or tail in SHELLS:  # (d) a runner, (b) a shell
            program = shell_scan.parse(shlex.join(words[at:])) or [] if tail in SHELLS else []
            return found + runs(words[at:], under, quoted) + [
                pair for s in program[1:] for pair in runs(s.argv)]
        if quoted and len(word.split()) > 1 and shell_scan:  # (c) a command string
            found += [pair for s in shell_scan.parse(word) or [] for pair in runs(s.argv)]
    return found


def publishes(words: list[str]) -> list[str]:
    """The PUBLISHES reasons of one simple command, read from its argv."""
    found, call = [], git_call(words)
    if call is not None:
        sub, args = call
        found += [PUSH] if sub == "push" else []
        found += [NO_VERIFY] if "--no-verify" in args else []
        if (sub in ("rebase", "filter-branch") or (sub == "reset" and "--hard" in args)
                or (sub == "commit" and "--amend" in args)):
            found.append(HISTORY)
        if sub == "push" and any(a.startswith("--force") or a == "-f" for a in args):
            found.append(FORCE)
    gh = gh_words(words)
    if gh and len(words) > 2:
        group, verb = [*gh, ""][:2]
        pr_writes = ("create", "new", "merge", "ready", "edit", "update-branch")
        found += [PR] if group == "pr" and verb in pr_writes else []
        if group == "release" or (group == "repo" and verb in ("create", "new", "delete", "edit")):
            found.append(RELEASE)
    return found


def emit(event: dict, verdict: str, reason: str, act: str, why: str) -> None:
    """Route the verdict through harness_mode; unable to import it, deny what would ask."""
    output = None
    try:
        import harness_mode

        output = harness_mode.guard(event, "guard_push", verdict, reason, act, why)
    except ImportError as error:
        if verdict != "pass":
            output = deny(f"guard_push cannot import harness_mode ({error}); refused: {why}.")
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
        return 0
    try:
        if event.get("tool_name") != "Bash":
            return 0
        command = (event.get("tool_input", {}) or {}).get("command", "")

        simples = shell_scan.parse(command) if shell_scan else None
        if simples is None:  # the text rules
            reasons = [why for pattern, why in PUBLISHES if pattern.search(command)]
            commands = [tokens(segment) for segment in SEGMENTS.split(command)]
        else:
            reasons, commands = [], [s.argv for s in simples]
        for words, bare in [pair for argv in commands for pair in runs(argv)]:
            reasons += publishes(words)
            why = segment_why(words, command)
            if why is not None and not (bare and why == UNKNOWN):  # a bare tail never gives G2's
                reasons.append(why)
        reasons = list(dict.fromkeys(reasons))
        if not reasons:
            return 0
        why = "; ".join(reasons)
        reason = (
            "This reaches beyond the working tree: " + why + ".\nPRINCIPLES.md E1 -- the "
            "working tree is the agent's, the remote is the user's. Approval of an earlier "
            "push does not carry to this one."
        )
        emit(event, "ask", reason, command, why)
    except Exception as error:  # a crash must not turn an ask into a pass
        print(json.dumps(deny(f"guard_push failed: {type(error).__name__}: {error}")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
