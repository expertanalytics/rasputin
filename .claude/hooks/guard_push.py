#!/usr/bin/env python3
"""PreToolUse: the working tree is the agent's, the remote is the user's.

.claude/REQUIRED-READING.md (the approval boundary) states this and the user
stated it directly -- "ensure one agent working at the time, and alert me before
pushing". It held because the agent remembered to ask. That is not a mechanism;
it is a habit, and the same session has forgotten comparable ones.

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
command run by another program is judged too (G7, `runs`), by every rule
but one: a git or gh word in the argv of a program that tools/shell_scan.py
does not strip as a wrapper (`caffeinate git myalias`) skips G2's
unknown-command ask, except under `parallel`. `env`, `timeout`, `xargs` and
the other wrappers are stripped, so G2 applies to them. A quoted command
string (`sh -c '…'`) is judged as its own line by the same rules, so
`sh -c 'git myalias'` asks and `sh -c 'caffeinate git myalias'` passes.

An amend whose HEAD is held by no remote-tracking ref or tag passes: nobody
else can have that commit. It comes last on a line of simple commands joined
only by `&&`, after at most one `cd <absolute path>` and any `git add`
(docs/increments/h29-review-rules.md §3.3). Anything the guard cannot
establish about such a line asks, and names that form.
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
try:  # h20: a branch Ola granted for an away window
    import push_grant
except ImportError:  # no grant then passes anything
    push_grant = None

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
#: A program that may start a process, in any language, strings and comments included:
#: a call named in a string (`methodcaller('system')`) starts one as well.
SPAWN = re.compile(
    r"\b(subprocess|system|popen|exec\w*|spawn\w*|fork\w*|pty|posix|modules|import_module|"
    r"importlib|__import__|runpy|ctypes|getattr|methodcaller|eval|compile|child_process|open3|"
    r"qx|syscall)\b", re.IGNORECASE)
#: Perl's, Ruby's and a shell's other ways to start one; Python has none of them.
SPAWN_MARKS = re.compile(r"%x|`|['\"]\||\|['\"]")


def python(word: str, command: str) -> bool:
    """`word` runs Python: it names Python, and the line sets no PATH and defines no `word`."""
    name = re.escape(word)
    redefines = rf"\bPATH=|\bfunction\s+{name}\b|\b{name}\s*(\(\s*\)|=)"
    return bool(shell_scan.PYTHON.fullmatch(word)) and not re.search(redefines, command)

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
#: The flags a `gh api graphql` read may carry, each with its value as the next word.
GRAPHQL_FLAGS = {"-f", "--raw-field", "-F", "--field", "-X", "--method", "-q", "--jq"}
#: A line the shell expands nothing in: outside single quotes only these characters.
#: Each step takes one character or one whole quoted string, so it matches in linear time.
EXPANDS_NOTHING = re.compile(r"(?:[A-Za-z0-9_./= \t|-]|'[^']*')*")
#: The prefixes of `--amend` git accepts (`--a` is ambiguous to git and fails).
AMENDS = {"--am", "--ame", "--amen", "--amend"}
#: `&&` outside single quotes: an even number of quotes follows it.
AND = re.compile(r"&&(?=(?:[^']*'[^']*')*[^']*$)")
AMEND_FORM = (
    "an amend of a commit no remote holds passes when it comes last, after at most "
    "`cd <absolute path>` and `git add` commands joined by `&&`, or as "
    "`git -C <absolute path> commit --amend …` alone; it may be redone in that form"
)
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


def read_query(text: str) -> bool:
    """A GraphQL document that only reads: it opens with `{` or `query`, and the letters
    `mutation` and `subscription` appear nowhere, strings and comments included."""
    opens = re.match(r"\s*(\{|query(?!\w))", text) is not None
    return opens and "mutation" not in text and "subscription" not in text


def graphql_read(words: list[str], rest: list[str], raw: str) -> bool:
    """A `gh api graphql` call, every word read before expansion, whose one `query` reads.

    The whole line `raw` expands nothing (EXPANDS_NOTHING, and no word opens with zsh's
    `=cmd`), so the words are what gh receives. Only the flags of GRAPHQL_FLAGS, each
    followed by its value; method POST if given; fields exactly one `query` and at most
    one `operationName`, none read from a file.
    """
    if EXPANDS_NOTHING.fullmatch(raw) is None or any(w.startswith("=") for w in words):
        return False
    if any("$" in w or "`" in w for w in words) or rest[:1] != ["graphql"] or len(rest) % 2 == 0:
        return False
    fields: list[tuple[str, str]] = []
    for flag, value in zip(rest[1::2], rest[2::2], strict=True):
        name, equals, given = value.partition("=")
        if flag not in GRAPHQL_FLAGS or (flag in ("-X", "--method") and value != "POST"):
            return False
        if flag in ("-f", "--raw-field", "-F", "--field"):
            if not equals or (flag in ("-F", "--field") and given.startswith("@")):
                return False
            fields.append((name, given))
    names = sorted(name for name, _ in fields)
    if names not in (["query"], ["operationName", "query"]):
        return False
    return read_query(dict(fields)["query"])


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
        writes = (given is not None and given != "GET") or (given is None and fields)
        api = words[words.index("api") + 1 :]
        if writes and not graphql_read(words, api, text or shlex.join(words)):
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
                or (sub == "commit" and AMENDS & set(args))):
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


def git_out(where: str, *args: str) -> str | None:
    """The stripped output of `git -C where args`, or None on a non-zero exit."""
    done = subprocess.run(["git", "-C", where, *args], capture_output=True, text=True,
                          check=False, timeout=10)  # fmt: skip
    return done.stdout.strip() if done.returncode == 0 else None


def unpushed_amend(command: str, cwd: str | None) -> bool:
    """True when the line ends in a bare amend whose HEAD no remote-tracking ref or tag holds.

    The line is simple commands joined only by `&&`, each part's argv its own words
    and expanding nothing: at most one `cd <absolute path>` first, any `git add`, and
    the amend last with no git option before `commit` but one absolute `-C`, so HEAD
    is the one the amend replaces. Anything unknown, an exception included, is False.
    """
    try:
        simples = shell_scan.parse(command) if shell_scan else None
        parts = AND.split(command)
        if not simples or len(simples) != len(parts) or any(
                EXPANDS_NOTHING.fullmatch(p) is None or shlex.split(p) != s.argv
                for p, s in zip(parts, simples, strict=True)):  # fmt: skip
            return False
        *before, argv = [s.argv for s in simples]
        if before and before[0][:1] == ["cd"]:
            cd = before.pop(0)
            if len(cd) != 2 or not Path(cd[1]).is_absolute():
                return False
            cwd = cd[1]
        if any(b[:2] != ["git", "add"] for b in before):
            return False
        head, rest = (argv[1:3], argv[3:]) if argv[1:2] == ["-C"] else ([], argv[1:])
        where = head[1] if head else cwd
        if (argv[:1] != ["git"]
                or rest[:1] != ["commit"] or not AMENDS & set(rest) or not isinstance(where, str)
                or not where or (head and not Path(where).is_absolute())):  # fmt: skip
            return False
        sha = git_out(where, "rev-parse", "--verify", "-q", "HEAD^{commit}")
        held = sha and git_out(where, "for-each-ref", "--contains", sha, "--count=1",
                               "--format=%(refname)", "refs/remotes", "refs/tags")  # fmt: skip
        return bool(sha) and held == ""
    except Exception:
        return False


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
        else:  # a runner fed text the scanner did not read, and a program that may spawn
            texts = [command] if any(s.unread for s in simples) else []
            texts += [s.program for s in simples if s.program and (SPAWN.search(s.program) or (
                not python(s.argv[0], command) and SPAWN_MARKS.search(s.program)))]
            reasons = [why for text in texts for pattern, why in PUBLISHES if pattern.search(text)]
            commands = [s.argv for s in simples]
        for words, bare in [pair for argv in commands for pair in runs(argv)]:
            reasons += publishes(words)
            why = segment_why(words, command)
            if why is not None and not (bare and why == UNKNOWN):  # a bare tail never gives G2's
                reasons.append(why)
        reasons = list(dict.fromkeys(reasons))
        if not reasons or (reasons == [HISTORY] and unpushed_amend(command, event.get("cwd"))):
            return 0
        why = "; ".join(reasons)
        if reasons == [HISTORY] and any(AMENDS & set(words) for words in commands):
            why += "; " + AMEND_FORM
        grant = push_grant.judge(event, command, reasons) if push_grant else None
        if grant is not None and grant.passed:
            return 0
        if grant is not None:
            why += f"; {grant.lane} branch {grant.branch}, but {grant.why}"
        reason = (
            "This reaches beyond the working tree: " + why + ".\n.claude/REQUIRED-READING.md, "
            "the approval boundary -- the working tree is the agent's, the remote is the "
            "user's. Approval of an earlier push does not carry to this one."
        )
        emit(event, "ask", reason, command, why)
    except Exception as error:  # a crash must not turn an ask into a pass
        print(json.dumps(deny(f"guard_push failed: {type(error).__name__}: {error}")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
