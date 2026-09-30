#!/usr/bin/env python3
"""PreToolUse: the working tree is the agent's, the remote is the user's.

PRINCIPLES.md E1 states this and the user stated it directly -- "ensure one agent
working at the time, and alert me before pushing". It held because the agent
remembered to ask. That is not a mechanism; it is a habit, and the same session
has forgotten comparable ones.

`gh pr create`, `gh pr merge` and `git push` all publish. So does a `git commit`
carrying `--no-verify`, which is not a publish but is the disabling of somebody
else's guard, and belongs to the user for the same reason. So do writes to refs,
remotes and git config, and `gh api` / `curl` with a writing method to the forge
(docs/increments/h3-unattended-u1.md §3.5).

`ask`, not `deny`: the user says yes constantly. The point is that they say it.
The exception is unattended mode (tools/harness_mode.py): while Ola is away the
ask becomes a `deny`, and the act is queued for him.
"""

import json
import re
import shlex
import sys
from pathlib import Path

PUBLISHES = (
    (re.compile(r"\bgit\b[^|;&]*\bpush\b"), "git push writes to the remote"),
    (re.compile(r"\bgh\s+pr\s+(create|merge|ready|edit)\b"), "gh pr changes a pull request"),
    (re.compile(r"\bgh\s+(release|repo\s+(create|delete|edit))\b"),
     "gh publishes or alters the repo"),
    (re.compile(r"--no-verify\b"), "--no-verify disables git's own hooks"),
    (re.compile(r"\bgit\s+(rebase|reset\s+--hard|filter-branch)\b|\bgit\s+commit\b.*--amend"),
     "this rewrites history, which is destructive once anything is published"),
    (re.compile(r"\bgit\b[^|;&]*\bpush\b.*(--force|-f)\b"),
     "a force push can discard the user's commits"),
)

FORGE_HOST = "github.com"
SEGMENTS = re.compile(r"&&|\|\||;|\||\n")
#: Options whose next token is their argument, not a positional (§3.5).
TAKES_ARG = {"-f", "--file", "--blob", "--type", "--default", "--comment", "-m"}
REMOTE_WRITES = {"add", "set-url", "rename", "remove", "rm", "set-head", "set-branches"}
CONFIG_WRITE_FLAGS = {"--add", "--unset", "--unset-all", "--replace-all", "--rename-section",
                      "--remove-section", "--edit", "-e"}
CONFIG_WRITE_VERBS = {"set", "unset", "rename-section", "remove-section", "edit"}
CONFIG_READS = {"--list", "-l", "get", "list"}
GH_FIELDS = {"-f", "-F", "--field", "--raw-field", "--input"}
CURL_DATA = {"-d", "-F", "--form", "--json", "-T", "--upload-file"}


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
    """(subcommand, its arguments) of a `git` in the segment, skipping -C/-c options."""
    if "git" not in words:
        return None
    at = words.index("git") + 1
    while at < len(words) and words[at].startswith("-"):
        at += 2 if words[at] in ("-C", "-c") else 1
    return (words[at], words[at + 1:]) if at < len(words) else None


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


def segment_why(words: list[str]) -> str | None:
    """The reason a segment writes a ref, a remote, git config or the forge, or None."""
    call = git_call(words)
    if call is not None:
        sub, args = call
        pos = positionals(args)
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
    if words[:2] == ["gh", "api"]:
        given = method(words[2:], "-X", "--method")
        fields = any(a in GH_FIELDS or a.split("=")[0] in GH_FIELDS for a in words[2:])
        if (given is not None and given != "GET") or (given is None and fields):
            return "gh api with a writing method changes the forge"
    if "curl" in words and any(FORGE_HOST in w for w in words):
        given = method(words, "-X", "--request")
        data = any(w in CURL_DATA or w.startswith("--data") for w in words)
        if (given is not None and given not in ("GET", "HEAD")) or data:
            return "curl with a writing method to the forge"
    return None


def emit(event: dict, verdict: str, reason: str, act: str, why: str) -> None:
    """Route the verdict through harness_mode; unable to import it, deny what would ask."""
    output = None
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
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

        reasons = [why for pattern, why in PUBLISHES if pattern.search(command)]
        for segment in SEGMENTS.split(command):
            why = segment_why(tokens(segment))
            if why is not None and why not in reasons:
                reasons.append(why)
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
