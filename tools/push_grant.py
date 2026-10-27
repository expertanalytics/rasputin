"""Away-mode push grants: the granted forms and the reviewed rule (harness increment h20).

Spec: docs/increments/h20-away-push-grants.md §3.2 (the three forms and the
expansion test), §3.3 (`reviewed`) and §3.5 (`judge`, called by
.claude/hooks/guard_push.py). While Ola is away, a branch he named with
`away.py --grant` may be pushed by sha, opened as a pull request and enqueued
once per head, provided `@reviewer`'s last recorded verdict approved its head.
A branch holding commits of `private/master` goes only to the private
repository, where it may also be merged once per head when its checks are
green; any other branch goes only to origin. Standard library only.
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
PRIVATE_REMOTE, PRIVATE_REPO = "private", "expertanalytics/rasputin-private"
PRIVATE_MASTER = "refs/remotes/private/master"
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
          "--draft": "draft", "-d": "draft", "-R": "repo", "--repo": "repo"}  # fmt: skip
MERGE = {"--match-head-commit": "sha", "--merge": "merge", "-R": "repo", "--repo": "repo"}
FLAGS = ("draft", "merge")  # options that take no value
VIEW = ("state", "baseRefName", "headRefOid", "statusCheckRollup")
UNREAD = "the pull request's checks could not be read"
UNSURE = "so it cannot be told whether this branch holds private commits; run git fetch private"
#: `@perf`'s dated evidence folders; a file a commit adds there counts like prose.
BENCH = re.compile(r"docs/benchmarks/\d{4}-\d{2}-\d{2}/.+", re.S)
NOT_BARE = ("the grant takes only one bare command, every word written out "
            "(no variable, no &&, ;, pipe, redirect or heredoc); run it alone; "
            "it may be redone in that form")  # fmt: skip


@dataclass(frozen=True)
class Act:
    kind: str  # "push", "pull request", "enqueue" (origin) or "merge" (private)
    branch: str
    sha: str | None = None
    directory: str | None = None
    remote: str = "origin"


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


def repo_remote(found: dict[str, str]) -> str:
    """The remote a gh `-R` option names: origin when absent, "" for any other repository."""
    return {None: "origin", PRIVATE_REPO: PRIVATE_REMOTE}.get(found.get("repo"), "")


def granted_act(command: str) -> Act | None:
    """The act of the granted forms that `command` is, alone and unexpanded; else None. Pure."""
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
        remote = rest[1] if rest[:1] == ["push"] and len(rest) == 3 else ""
        pinned = SHA.fullmatch(sha) and ref != branch and absolute
        if remote in ("origin", PRIVATE_REMOTE) and pinned:
            act = Act("push", branch, sha, directory, remote) if branch_name(branch) else None
    elif argv[:3] == ["gh", "pr", "create"]:
        parsed = options(argv[3:], CREATE)
        if parsed is not None and not parsed[1]:
            found, remote = parsed[0], repo_remote(parsed[0])
            base = found.get("base", "master")
            if remote and base == "master" and branch_name(found.get("head", "")):
                act = Act("pull request", found["head"], remote=remote)
    elif argv[:3] == ["gh", "pr", "merge"]:
        parsed = options(argv[3:], MERGE)
        if parsed is not None and len(parsed[1]) == 1:
            branch, sha, remote = parsed[1][0], parsed[0].get("sha", ""), repo_remote(parsed[0])
            # Without a merge queue, gh would ask for a method; the private form names it.
            method = remote == "origin" or "merge" in parsed[0]
            named = branch_name(branch) and not branch.isdigit() and SHA.fullmatch(sha)
            if remote and method and named:
                kind = "enqueue" if remote == "origin" else "merge"
                act = Act(kind, branch, sha, remote=remote)
    return act if act is not None and unexpandable(command, argv) else None


def loose_act(command: str, granted: tuple[str, ...]) -> tuple[str, str] | None:
    """(branch, kind) of a push, pull request or merge in `command` naming a granted branch.

    Loose on purpose: it only chooses the words of a refusal, never passes anything.
    """
    for simple in shell_scan.parse(command) or []:
        argv = simple.argv
        verb = argv[2] if argv[:2] == ["gh", "pr"] and len(argv) > 2 else ""
        kind = {"create": "pull request", "merge": "enqueue"}.get(verb, "")
        kind = "push" if argv[:1] == ["git"] and "push" in argv else kind
        for word in argv[1:] if kind else []:
            ref = word.rpartition(":")[2].removeprefix("refs/heads/")
            named = [b for b in (word, word.removeprefix("--head="), ref) if b in granted]
            if named:
                return named[0], kind
    return None


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


def git(
    repo: Path | str, *args: str, codes: tuple[int, ...] = (0,), timeout: float | None = None
) -> str:
    done = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True,
                          check=False, timeout=timeout)  # fmt: skip
    if done.returncode not in codes:
        raise subprocess.CalledProcessError(done.returncode, args, done.stdout, done.stderr)
    return done.stdout


def home(repo: Path | str, head: str) -> tuple[str, str]:
    """(remote, "") for the repository `head` belongs to, or ("", the reason it cannot be told).

    Fails closed: while the private remote is set up, a head is public only when
    `private/master` was read and the head holds none of its private-only commits.
    A config read error is never taken for "not set up".
    """
    try:
        argv = ["git", "-C", str(repo), "config", "--get", "remote.private.url"]
        config = subprocess.run(argv, capture_output=True, text=True, check=False, timeout=10)
    except subprocess.TimeoutExpired:
        return "", ("git config did not answer within 10 s, so it cannot be told whether the "
                    "private remote is set up")  # fmt: skip
    if config.returncode not in (0, 1):
        code = config.returncode
        return "", (f"git config could not be read (it exited {code}), so it cannot be told "
                    "whether the private remote is set up; fix the file git names")  # fmt: skip
    url = config.returncode == 0
    try:
        if not git(repo, "rev-parse", "--verify", "-q", PRIVATE_MASTER, codes=(0, 1)).strip():
            missing = f"the private remote is set up but {PRIVATE_MASTER} is missing, {UNSURE}"
            return ("", missing) if url else ("origin", "")
        count = ("rev-list", "--count", PRIVATE_MASTER, "^refs/remotes/origin/master")
        n, m = (int(git(repo, *count, *extra, timeout=10)) for extra in ((), (f"^{head}",)))
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, ValueError):
        return "", "the private and public histories could not be compared"
    if n == 0:
        none = f"{PRIVATE_MASTER} holds no commit that origin/master lacks, {UNSURE}"
        return ("", none) if url else ("origin", "")
    return (PRIVATE_REMOTE if m < n else "origin"), ""


def changes(repo: Path | str, commit: str, remote: str = "origin") -> str | None:
    """What a commit after the review changes beyond prose and added dated evidence, or None."""
    parents = git(repo, "rev-list", "--parents", "-n", "1", commit).split()[1:]
    if len(parents) == 1:
        base = parents[0]
    elif len(parents) == 2 and subprocess.run(
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", parents[1],
         f"refs/remotes/{remote}/master"], capture_output=True, check=False,
    ).returncode == 0:  # fmt: skip
        base = git(repo, "merge-tree", "--write-tree", *parents, codes=(0, 1)).splitlines()[0]
    else:
        return f"is neither a prose commit nor a merge of {remote}/master"
    paths = git(repo, "diff", "--name-only", "--no-renames", "-z", base, commit).split("\0")
    added = git(repo, "diff", "--name-only", "--no-renames", "--diff-filter=A", "-z", base, commit)
    evidence = {path for path in added.split("\0") if BENCH.fullmatch(path)}
    code = [p for p in paths if p and not ci_changes.is_prose(p) and p not in evidence]
    return f"changes {code[0]}" if code else None


def review_of(
    repo: Path | str, head: str, records: list[Any], remote: str = "origin"
) -> tuple[str | None, dict[str, Any] | None]:
    """(None, the record relied on) when `head` is reviewed; else (the reason, None)."""
    try:
        chain = git(repo, "rev-list", "--first-parent", f"refs/remotes/{remote}/master..{head}")
        commits = chain.split()
        record = next((r for r in reversed(records)
                       if isinstance(r, dict) and r.get("sha") in commits), None)  # fmt: skip
        if record is None:
            return NO_REVIEW, None
        if record.get("verdict") not in PASSING:
            return f"last review round: {record.get('verdict')}", None
        reviewed_at = record["sha"]
        for commit in reversed(commits[: commits.index(reviewed_at)]):
            what = changes(repo, commit, remote)
            if what is not None:
                return f"commit {commit[:7]} {what} after the review of {reviewed_at[:7]}", None
    except subprocess.CalledProcessError as error:
        return f"git failed: {' '.join(error.cmd)}: {str(error.stderr).strip()}", None
    return None, record


def reviewed(repo: Path | str, head: str, records: list[Any], remote: str = "origin") -> str | None:
    """None when `head` is reviewed (§3.3), else the reason in plain words."""
    return review_of(repo, head, records, remote)[0]


def checks(view: Any, head: str) -> str | None:
    """Why a private merge's pull request (`gh pr view --json` of VIEW) is not green, or None."""
    if not isinstance(view, dict) or any(key not in view for key in VIEW):
        return UNREAD
    rollup = view["statusCheckRollup"]
    if view["state"] != "OPEN":
        return f"the pull request is {view['state']}, not open"
    if view["baseRefName"] != "master":
        return f"the pull request is into {view['baseRefName']}, not master"
    if view["headRefOid"] != head:
        return f"the pull request's head is {str(view['headRefOid'])[:7]}, not {head[:7]}"
    if not isinstance(rollup, list) or not rollup:
        return "no checks reported" if rollup == [] else UNREAD
    for item in rollup:
        if not isinstance(item, dict):
            return UNREAD
        if item.get("__typename") == "StatusContext":
            green = item.get("state") == "SUCCESS"
        else:
            green = item.get("status") == "COMPLETED" and item.get("conclusion") in (
                "SUCCESS", "SKIPPED", "NEUTRAL")  # fmt: skip
        if not green:
            name = item.get("name") or item.get("context")
            seen = item.get("state") or item.get("conclusion") or item.get("status")
            return f"check {name} is {seen}" if name and seen else UNREAD
    if not any("SUCCESS" in (i.get("state"), i.get("conclusion")) for i in rollup):
        return "no check succeeded: every check was skipped or neutral"
    return None


def pr_view(where: str, branch: str) -> Any:
    """The private pull request's state and checks as gh reports them; None when unread."""
    argv = ["gh", "pr", "view", branch, "-R", PRIVATE_REPO, "--json", ",".join(VIEW)]
    try:
        done = subprocess.run(argv, cwd=where, capture_output=True, text=True, check=False,
                              timeout=20)  # fmt: skip
        return json.loads(done.stdout) if done.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired, ValueError):
        return None


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
    if act is None:
        loose = loose_act(command, grants(flag))
        return Grant(False, *loose, None, NOT_BARE) if loose else None
    if act.branch not in grants(flag):
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
    remote, why = home(where, head)
    if why:
        return refused(why, head)
    if remote != act.remote:
        holds, only = (("holds commits of", "the private repository") if act.remote == "origin"
                       else ("holds no commit of", "origin"))  # fmt: skip
        return refused(f"head {head[:7]} {holds} private/master, so it goes only to {only}", head)
    reason, record = review_of(where, head, jsonl(state / "reviews.jsonl"), act.remote)
    if reason is not None or record is None:
        return refused(
            f"head {head[:7]} is not reviewed: {reason}; it may be redone once "
            f"@reviewer's verdict approving {head[:7]} is recorded",
            head,
        )
    since = flag.get("since")
    if act.kind in ("enqueue", "merge") and any(
        isinstance(u, dict) and (u.get("window"), u.get("branch"), u.get("act"), u.get("head"))
        == (since, act.branch, act.kind, head) for u in jsonl(state / "grants-used.jsonl")
    ):  # fmt: skip
        done = "enqueued" if act.kind == "enqueue" else "merged"
        return refused(f"head {head[:7]} was {done} already in this window", head)
    if act.kind == "merge" and (why := checks(pr_view(where, act.branch), head) or ""):
        return refused(why, head)
    line = {"at": datetime.now(UTC).isoformat(timespec="seconds"), "window": since,
            "branch": act.branch, "act": act.kind, "head": head, "remote": act.remote,
            "review": [record.get("sha"), record.get("at")]}  # fmt: skip
    error = hm.append_queue(state, line, "grants-used.jsonl")
    if error is not None:
        return refused(f"the grants log cannot be written ({error})", head)
    return Grant(True, act.branch, act.kind, head, "")
