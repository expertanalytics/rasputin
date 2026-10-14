"""Away-mode push grants: the granted forms and the reviewed rule (harness increment h20).

Spec: docs/increments/h20-away-push-grants.md §3.2 (the three forms and the
expansion test), §3.3 (`reviewed`) and §3.5 (`judge`, called by
.claude/hooks/guard_push.py). While Ola is away, a branch he named with
`away.py --grant` may be pushed by sha, opened as a pull request and enqueued
once per head, provided `@reviewer`'s last recorded verdict approved its head.
Standard library only.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

sys.path.append(str(Path(__file__).resolve().parent))  # appended: the stdlib wins (h16 G4)
import ci_changes
import harness_mode as hm
import shell_scan

#: guard_push's reason strings for a push and a pull request write.
PUSH, PR = "git push writes to the remote", "gh pr changes a pull request"
SHA = re.compile(r"[0-9a-f]{40}")
BRANCH = re.compile(r"[A-Za-z0-9._/-]+")
WORD = re.compile(r"[A-Za-z0-9._/:@-]+")
LONG = re.compile(r"(--[a-z-]+)=(.*)", re.S)
REVIEW = re.compile(
    r"^Review: (APPROVED|APPROVED WITH WORDINGS|CHANGES REQUESTED) ([0-9a-f]{40})$", re.M
)
PASSING = ("APPROVED", "APPROVED WITH WORDINGS")
NOT_GRANTABLE = ("master", "main", "HEAD")
NO_REVIEW = "no review recorded for this branch"
#: gh pr create's options: each spelling to its name; `draft` takes no value.
CREATE = {"--head": "head", "-H": "head", "--base": "base", "-B": "base",
          "--title": "title", "-t": "title", "--body-file": "body", "-F": "body",
          "--draft": "draft", "-d": "draft"}  # fmt: skip
FLAGS = ("draft", "merge")  # options that take no value


@dataclass(frozen=True)
class Act:
    kind: str  # "push", "pull request" or "enqueue"
    branch: str
    sha: str | None = None
    directory: str | None = None


@dataclass(frozen=True)
class Grant:
    passed: bool
    branch: str
    act: str
    head: str | None
    why: str


def branch_name(name: str) -> bool:
    return BRANCH.fullmatch(name) is not None and name not in NOT_GRANTABLE


def options(words: list[str], names: dict[str, str]) -> tuple[dict[str, str], list[str]] | None:
    """Each named option once (`--opt value`, `--opt=value`, or a flag), and the positionals."""
    found: dict[str, str] = {}
    rest: list[str] = []
    at = 0
    while at < len(words):
        word, at = words[at], at + 1
        glued = LONG.fullmatch(word)
        spelling, value = (glued[1], glued[2]) if glued else (word, None)
        name = names.get(spelling)
        if name is None:
            if word.startswith("-"):
                return None
            rest.append(word)
            continue
        if name in found or (name in FLAGS and value is not None):
            return None
        if name not in FLAGS and value is None:
            if at == len(words):
                return None
            value, at = words[at], at + 1
        found[name] = value or ""
    return found, rest


def unexpandable(command: str, argv: list[str]) -> bool:
    """No word the shell could expand; the title alone stands in single quotes (§3.2)."""
    try:
        lexer = shlex.shlex(command, posix=False)
        lexer.whitespace_split = True
        raw = list(lexer)
    except ValueError:
        return False
    if len(raw) != len(argv):
        return False
    after_title = False
    for text, word in zip(raw, argv, strict=True):
        glued = LONG.fullmatch(word)
        prefix, value = (glued[1] + "=", glued[2]) if glued else ("", word)
        if after_title or prefix == "--title=":
            if "'" in value or text != f"{prefix}'{value}'":
                return False
        elif text != word or not WORD.fullmatch(value):
            return False
        after_title = word in ("--title", "-t")
    return True


def granted_act(command: str) -> Act | None:
    """The act of §3.2's three forms that `command` is, alone and unexpanded; else None. Pure."""
    simples = shell_scan.parse(command)
    try:
        if simples is None or len(simples) != 1 or shlex.split(command) != simples[0].argv:
            return None
    except ValueError:
        return None
    argv = simples[0].argv
    act = None
    if argv[:1] == ["git"]:
        directory = argv[2] if argv[1:2] == ["-C"] and len(argv) > 2 else None
        rest = argv[3:] if directory is not None else argv[1:]
        sha, _, ref = rest[2].partition(":") if len(rest) == 3 else ("", "", "")
        branch = ref.removeprefix("refs/heads/")
        absolute = directory is None or directory.startswith("/")
        if rest[:2] == ["push", "origin"] and SHA.fullmatch(sha) and ref != branch and absolute:
            act = Act("push", branch, sha, directory) if branch_name(branch) else None
    elif argv[:3] == ["gh", "pr", "create"]:
        parsed = options(argv[3:], CREATE)
        if parsed is not None and not parsed[1]:
            found = parsed[0]
            if branch_name(found.get("head", "")) and found.get("base", "master") == "master":
                act = Act("pull request", found["head"])
    elif argv[:3] == ["gh", "pr", "merge"]:
        names = {"--match-head-commit": "sha", "--merge": "merge"}
        parsed = options(argv[3:], names)
        if parsed is not None and len(parsed[1]) == 1:
            branch, sha = parsed[1][0], parsed[0].get("sha", "")
            if branch_name(branch) and not branch.isdigit() and SHA.fullmatch(sha):
                act = Act("enqueue", branch, sha)
    return act if act is not None and unexpandable(command, argv) else None


def grants(flag: dict[str, Any]) -> tuple[str, ...]:
    """The flag's granted branches; anything but a list of strings grants nothing."""
    value = flag.get("grants")
    if isinstance(value, list) and all(isinstance(item, str) for item in value):
        return tuple(value)
    return ()


def parse_review(text: str) -> tuple[str, str] | None:
    """(verdict, sha) of the one `Review:` line of a handback, or None for none or several."""
    found = REVIEW.findall(text)
    return (found[0][0], found[0][1]) if len(found) == 1 else None


def git(repo: Path | str, *args: str, codes: tuple[int, ...] = (0,)) -> str:
    done = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True,
                          check=False)  # fmt: skip
    if done.returncode not in codes:
        raise subprocess.CalledProcessError(done.returncode, args, done.stdout, done.stderr)
    return done.stdout


def changes(repo: Path | str, commit: str) -> str | None:
    """What a commit after the review changes beyond prose, or None (§3.3 item 4)."""
    parents = git(repo, "rev-list", "--parents", "-n", "1", commit).split()[1:]
    if len(parents) == 1:
        base = parents[0]
    elif len(parents) == 2 and subprocess.run(
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", parents[1],
         "refs/remotes/origin/master"], capture_output=True, check=False,
    ).returncode == 0:  # fmt: skip
        base = git(repo, "merge-tree", "--write-tree", *parents, codes=(0, 1)).splitlines()[0]
    else:
        return "is neither a prose commit nor a merge of origin/master"
    paths = git(repo, "diff", "--name-only", "--no-renames", "-z", base, commit).split("\0")
    code = [path for path in paths if path and not ci_changes.is_prose(path)]
    return f"changes {code[0]}" if code else None


def review_of(
    repo: Path | str, head: str, records: list[Any]
) -> tuple[str | None, dict[str, Any] | None]:
    """(None, the record relied on) when `head` is reviewed; else (the reason, None)."""
    try:
        chain = git(repo, "rev-list", "--first-parent", f"refs/remotes/origin/master..{head}")
        commits = chain.split()
        record = next((r for r in reversed(records)
                       if isinstance(r, dict) and r.get("sha") in commits), None)  # fmt: skip
        if record is None:
            return NO_REVIEW, None
        if record.get("verdict") not in PASSING:
            return f"last review round: {record.get('verdict')}", None
        reviewed_at = record["sha"]
        for commit in reversed(commits[: commits.index(reviewed_at)]):
            what = changes(repo, commit)
            if what is not None:
                return f"commit {commit[:7]} {what} after the review of {reviewed_at[:7]}", None
    except subprocess.CalledProcessError as error:
        return f"git failed: {' '.join(error.cmd)}: {str(error.stderr).strip()}", None
    return None, record


def reviewed(repo: Path | str, head: str, records: list[Any]) -> str | None:
    """None when `head` is reviewed (§3.3), else the reason in plain words."""
    return review_of(repo, head, records)[0]


def jsonl(path: Path) -> list[Any]:
    """The parsed lines of a JSON-lines file, oldest first; unreadable lines skipped."""
    if not path.exists():
        return []
    found = []
    for line in path.read_text(errors="replace").splitlines():
        try:
            found.append(json.loads(line))
        except ValueError:
            continue
    return found


def judge(event: dict[str, Any], command: str, reasons: list[str]) -> Grant | None:
    """A granted branch's act, passed (and logged) or refused with why; None if not one (§3.5)."""
    state = hm.state_dir(hm.ROOT)
    if reasons not in ([PUSH], [PR]) or state is None:
        return None
    if hm.read_mode(state, datetime.now(UTC)).state != "on":
        return None
    flag = json.loads((state / "unattended.json").read_text())
    act = granted_act(command)
    if act is None or act.branch not in grants(flag):
        return None

    def refused(why: str, head: str | None = None) -> Grant:
        return Grant(False, act.branch, act.kind, head, why)

    where = act.directory or event.get("cwd") or os.getcwd()
    common = git(where, "rev-parse", "--path-format=absolute", "--git-common-dir").strip()
    if Path(common).resolve() != state.parent.resolve():
        return refused("the command runs in another repository")
    try:
        head = git(where, "rev-parse", "--verify", f"refs/heads/{act.branch}").strip()
    except subprocess.CalledProcessError:
        return refused("the branch does not exist here")
    if act.sha is not None and act.sha != head:
        return refused(f"the command names {act.sha[:7]}, not the head {head[:7]}", head)
    reason, record = review_of(where, head, jsonl(state / "reviews.jsonl"))
    if reason is not None or record is None:
        return refused(f"head {head[:7]} is not reviewed: {reason}", head)
    since = flag.get("since")
    if act.kind == "enqueue" and any(
        isinstance(u, dict) and (u.get("window"), u.get("branch"), u.get("act"), u.get("head"))
        == (since, act.branch, "enqueue", head) for u in jsonl(state / "grants-used.jsonl")
    ):  # fmt: skip
        return refused(f"head {head[:7]} was enqueued already in this window", head)
    line = {"at": datetime.now(UTC).isoformat(timespec="seconds"), "window": since,
            "branch": act.branch, "act": act.kind, "head": head,
            "review": [record.get("sha"), record.get("at")]}  # fmt: skip
    error = hm.append_queue(state, line, "grants-used.jsonl")
    if error is not None:
        return refused(f"the grants log cannot be written ({error})", head)
    return Grant(True, act.branch, act.kind, head, "")
