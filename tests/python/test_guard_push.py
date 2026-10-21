"""`.claude/hooks/guard_push.py` in and out of unattended mode (h3 U1, T4, T5, T8-T10; h10 §4a).

The spec is `docs/increments/h3-unattended-u1.md` §3.5 (the new patterns),
§3.9 (failure direction) and §4. Each command is run through a copy of the hook
in a temporary repository whose flag the test sets, so the hook's state dir is
`<tmp>/.git/harness`. By day an asked act asks; at night it is denied and
queued; a pass is silent in both.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import NamedTuple

import pytest

from harness_fixtures import (
    GUARD_GOVERNANCE,
    GUARD_PUSH,
    REAL,
    SUBAGENT,
    Granted,
    Tool,
    add_worktree,
    bash_event,
    clean_env,
    commit_file,
    file_event,
    flag_granting,
    flag_path,
    git,
    grants_used_path,
    iso,
    jsonl,
    load_tool,
    make_repo,
    point_scratchpad,
    pretool_decision,
    queue_lines,
    run_script,
    set_mode,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

REMOTE = "this changes where the remote points"
CONFIG = "this writes git configuration (hooks path, remote URLs)"
GH_API = "gh api with a writing method changes the forge"
CURL = "curl with a writing method to the forge"
GH_PR = "gh pr changes a pull request"
PUSH = "git push writes to the remote"

#: T4's ask rows, each with the `why` §3.5 gives it (None: an existing pattern,
#: whose reason is unchanged).
ASKED: dict[str, str | None] = {
    "git update-ref refs/remotes/origin/master HEAD": "update-ref moves a ref directly",
    "git -C /x update-ref -d refs/heads/y": "update-ref moves a ref directly",
    "git remote add up u": REMOTE,
    "git remote set-url origin u": REMOTE,
    "git remote rename a b": REMOTE,
    "git remote remove a": REMOTE,
    "git config core.hooksPath x": CONFIG,
    "git config --unset remote.origin.url": CONFIG,
    "git config --global --add a.b c": CONFIG,
    "git config set a.b c": CONFIG,
    "git symbolic-ref HEAD refs/heads/x": "symbolic-ref rewrites a symbolic ref",
    "git symbolic-ref -d HEAD": "symbolic-ref rewrites a symbolic ref",
    "gh api -X PUT repos/o/r/pulls/1/merge": GH_API,
    "gh api --method=POST repos/o/r/pulls": GH_API,
    "gh api -XDELETE repos/o/r/git/refs/heads/x": GH_API,
    "gh api repos/o/r/pulls -f title=t": GH_API,
    "gh api graphql -f query=q": GH_API,
    "curl -X POST https://api.github.com/repos/o/r/pulls": CURL,
    "curl -d x https://api.github.com/x": CURL,
    "echo ok && git update-ref a b": "update-ref moves a ref directly",
    "git push": None,
    "gh pr merge 1": None,
    # h10 §4a: update-branch writes to the PR's branch on the remote; --auto
    # enqueues into the merge queue, the same act as a plain merge.
    "gh pr update-branch 1": GH_PR,
    "gh pr update-branch --rebase 1": GH_PR,
    "gh pr merge --auto 1": GH_PR,
    # h27: `git commit --amend` moved to the pushed-HEAD rows below (an amend
    # in this fixture, which has no remote, passes).
    "git rebase main": None,
}

#: T4's pass rows: read forms of the same commands, and requests to other hosts.
PASSED = (
    "git remote -v",
    "git remote get-url origin",
    "git config user.email",
    "git config --get remote.origin.url",
    "git config --list",
    "git symbolic-ref --short HEAD",
    "gh api repos/o/r/pulls",
    "gh api -X GET repos/o/r/pulls -f state=open",
    "curl https://api.github.com/x",
    "curl -X POST https://example.org/x",
    "git status",
    # h10 §4a: the update-branch ask must not widen to the read-only gh pr verbs.
    "gh pr checks 1",
    "gh pr view 1",
)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


# ---------------------------------------------------------------- T4


@pytest.mark.parametrize("mode", ["off", "expired"])
@pytest.mark.parametrize("command", ASKED)
def test_an_asked_command_asks_while_attended(repo: Path, mode: str, command: str) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert reason.startswith("This reaches beyond the working tree: ")
    why = ASKED[command]
    if why is not None:
        assert why in reason
    assert queue_lines(repo) == []


@pytest.mark.parametrize("command", ASKED)
def test_an_asked_command_is_denied_and_queued_while_unattended(repo: Path, command: str) -> None:
    set_mode(repo, "on")
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "deny"
    assert reason.startswith("Refused: unattended mode is on until ")
    assert f"; {command} waits for Ola (" in reason
    why = ASKED[command]
    if why is not None:
        assert why in reason
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_push", command)
    assert (line["agent_type"], line["agent_id"]) == (None, None)


@pytest.mark.parametrize("mode", ["off", "on", "expired"])
@pytest.mark.parametrize("command", PASSED)
def test_a_passed_command_is_silent_in_every_mode(repo: Path, mode: str, command: str) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
def test_a_non_bash_tool_is_silent(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    event = file_event(repo, "Write", str(repo / "notes.txt"))
    assert pretool_decision(run_script(repo, GUARD_PUSH, event)) is None


# ---------------------------------------------------------------- T5


def test_a_push_with_a_broken_flag_is_denied_with_the_broken_reason(repo: Path) -> None:
    set_mode(repo, "broken")
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, "git push")))
    assert found is not None
    kind, reason = found
    assert kind == "deny"
    assert "cannot be read" in reason
    assert "away.py --back" in reason
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_push", "git push")


# ---------------------------------------------------------------- T8


def test_a_subagent_push_asks_while_attended(repo: Path) -> None:
    event = bash_event(repo, "git push", **SUBAGENT)
    found = pretool_decision(run_script(repo, GUARD_PUSH, event))
    assert found is not None and found[0] == "ask"


def test_a_subagent_push_while_unattended_is_queued_with_its_identity(repo: Path) -> None:
    set_mode(repo, "on")
    event = bash_event(repo, "git push", **SUBAGENT)
    found = pretool_decision(run_script(repo, GUARD_PUSH, event))
    assert found is not None and found[0] == "deny"
    [line] = queue_lines(repo)
    assert (line["agent_id"], line["agent_type"]) == ("a1b2c3", "tester")


# ---------------------------------------------------------------- T9


@pytest.mark.parametrize("hook", [GUARD_PUSH, GUARD_GOVERNANCE])
@pytest.mark.parametrize("mode", ["off", "on"])
def test_a_guard_that_crashes_on_the_event_denies(repo: Path, hook: str, mode: str) -> None:
    set_mode(repo, mode)
    event = {
        "hook_event_name": "PreToolUse",
        "tool_name": "Bash",
        "tool_input": ["git push > CLAUDE.md"],
        "cwd": str(repo),
    }
    found = pretool_decision(run_script(repo, hook, event))
    assert found is not None, "a crash turned into a pass"
    kind, reason = found
    assert kind == "deny"
    assert reason.strip()


@pytest.mark.parametrize("hook", [GUARD_PUSH, GUARD_GOVERNANCE])
@pytest.mark.parametrize("mode", ["off", "on"])
def test_stdin_that_is_not_json_prints_nothing(repo: Path, hook: str, mode: str) -> None:
    set_mode(repo, mode)
    result = run_script(repo, hook, "git push {")
    assert result.returncode == 0
    assert result.stdout == ""
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- T10


def _without_harness_mode(repo: Path) -> None:
    (repo / "tools" / "harness_mode.py").unlink(missing_ok=True)


def test_guard_push_without_harness_mode_denies_what_it_would_ask(repo: Path) -> None:
    _without_harness_mode(repo)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, "git push")))
    assert found is not None
    kind, reason = found
    assert kind == "deny"
    assert "harness_mode" in reason


def test_guard_push_without_harness_mode_still_passes_a_read(repo: Path) -> None:
    _without_harness_mode(repo)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, "git status"))) is None


def test_guard_governance_without_harness_mode_denies_what_it_would_ask(repo: Path) -> None:
    _without_harness_mode(repo)
    event = file_event(repo, "Edit", str(repo / "CLAUDE.md"))
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found is not None
    kind, reason = found
    assert kind == "deny"
    assert "harness_mode" in reason


# ---------------------------------------------------------------- h16 G1
#
# docs/increments/h16-harness-fixes.md §2 G1: a fetch or pull whose refspec
# names a destination writes a ref, and `git replace` changes what git reads
# for an object, unless it only lists.

FETCH_WRITE = "fetch writes a named ref"
REPLACE = "replace refs change what git reads for an object"

G1_ASKED: dict[str, str] = {
    "git fetch . abc:refs/remotes/origin/master": FETCH_WRITE,
    "git fetch origin +master:refs/heads/x": FETCH_WRITE,
    "git pull . a:b": FETCH_WRITE,
    "git pull origin a:b": FETCH_WRITE,
    # The pinned false positive: `1` is taken for the repository, so the URL is
    # the second positional and its `:` reads as a refspec. `--depth=1` is not.
    "git fetch --depth 1 git@github.com:a/b.git": FETCH_WRITE,
    "git replace HEAD HEAD~1": REPLACE,
    "git replace -d x": REPLACE,
    "git replace -f a b": REPLACE,
    "git replace --graft a b": REPLACE,
    "git replace --edit a": REPLACE,
    "git replace --convert-graft-file": REPLACE,
    # The amendment after review round 6 (Ola's option C): a fetch, a pull or a
    # `git remote update` writes a named ref with no `src:dst` word when given
    # `--refmap` or `--stdin` (any prefix git takes, down to `--ref` and `--st`),
    # a `remote.*` or `url.*` override among git's own options (`-c`,
    # `--config-env`, either case), or `GIT_CONFIG*` in the line's text.
    "git fetch --refmap=+a:refs/heads/x o master": FETCH_WRITE,
    "git fetch --refmap '+a:refs/heads/x' o master": FETCH_WRITE,
    "git fetch --refm=+a:b o master": FETCH_WRITE,
    "git pull --ref=+a:b o master": FETCH_WRITE,
    "git fetch --stdin o": FETCH_WRITE,
    "git fetch --st o": FETCH_WRITE,
    "git -c remote.s.fetch=+a:b fetch s": FETCH_WRITE,
    "git -c REMOTE.s.FETCH=+a:b fetch s": FETCH_WRITE,
    "git -c remote.s.url=/x fetch s": FETCH_WRITE,
    "git -c url./x.insteadOf=https://github.com/ fetch origin": FETCH_WRITE,
    "git --config-env=remote.s.fetch=RF fetch s": FETCH_WRITE,
    "git --config-env remote.s.fetch=RF fetch s": FETCH_WRITE,
    "git -c remote.s.fetch=a:b pull s": FETCH_WRITE,
    "git -c remote.s.fetch=a:b remote update s": FETCH_WRITE,
    "GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=remote.s.fetch GIT_CONFIG_VALUE_0=a:b git fetch s": (
        FETCH_WRITE
    ),
}

G1_PASSED = (
    "git fetch origin",
    "git fetch -q origin master",
    "git fetch --all",
    "git fetch --depth=1 git@github.com:a/b.git",
    "git replace",
    "git replace -l",
    "git replace --list 'a*'",
    "git replace --format=short",
    "git replace -l --format=long 'a*'",
    # The amendment's passes: a different long option, a `-c` with another key,
    # a `remote update` with no override, and an override on a command that
    # fetches nothing.
    "git fetch --refetch origin",
    "git -c protocol.version=2 fetch origin",
    "git remote update",
    "git -c remote.s.fetch=a:b status",
)


# ---------------------------------------------------------------- h16 G2
#
# §2 G2: an alias hides what runs, so a git subcommand that is not a current
# git command (`git --list-cmds=main` less `--list-cmds=deprecated`), and a gh
# command outside gh's fixed top-level set, asks.

UNKNOWN = "cannot see what it runs"

G2_ASKED = (
    "git -c alias.p=push p origin",
    "git p origin",
    "git whatchanged",  # deprecated, so an alias may take its name
    "git pack-redundant",  # the other deprecated name in git 2.55
    "gh pm 12",
    "gh co 12",  # gh's own alias, which a user can redefine
)

G2_PASSED = (
    "git status",
    "git log -1",
    "git worktree list",
    "git -C /x log -1",
    "gh pr view 12",
    "gh api repos/x",
    "gh run list",
    "gh auth status",
    # Ola ruled 2026-10-05 ("yes, add help"): gh's own help command reads only.
    "gh help",
    "gh help pr",
)


@pytest.mark.parametrize("mode", ["off", "expired"])
@pytest.mark.parametrize("command", [*G1_ASKED, *G2_ASKED])
def test_a_named_ref_write_or_an_unknown_command_asks_while_attended(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert reason.startswith("This reaches beyond the working tree: ")
    assert G1_ASKED.get(command, UNKNOWN) in reason
    assert queue_lines(repo) == []


@pytest.mark.parametrize("command", [*G1_ASKED, *G2_ASKED])
def test_a_named_ref_write_or_an_unknown_command_is_queued_while_unattended(
    repo: Path, command: str
) -> None:
    set_mode(repo, "on")
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "deny"
    assert G1_ASKED.get(command, UNKNOWN) in reason
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_push", command)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", [*G1_PASSED, *G2_PASSED])
def test_a_listing_a_plain_fetch_or_a_known_command_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


def test_git_lists_its_main_commands() -> None:
    """G2 rests on `--list-cmds`; a git without it must fail here, not ask on everything."""
    listed = subprocess.run(
        ["git", "--list-cmds=main"], capture_output=True, text=True, env=clean_env(), check=True
    ).stdout.split()
    assert "push" in listed
    assert "status" in listed


# ---------------------------------------------------------------- h16 G3
#
# §2 G3: G3b, which passed local ref, remote and config writes made with
# `git -C <dir>` in a scratchpad repository, was dropped by Ola's option C
# after review round 6. A git write in a scratch repository now asks like any
# other; the rows that G3b passed are kept here, asked, with review round 6's
# spellings that slipped past it (`--glo`, a glued `-f<path>`).


@pytest.fixture
def scratch(repo: Path, tmp_path: Path) -> dict[str, Path]:
    """The copy's scratchpad, a repository in it, and a linked worktree of an outside one."""
    pad = point_scratchpad(repo, tmp_path / "faketmp")
    inside = pad / "r"
    git(make_plain_repo(inside), "commit", "-q", "--allow-empty", "-m", "two")
    outside = make_plain_repo(tmp_path / "outside")
    linked = add_worktree(outside, pad / "linked", "wt-linked")
    return {"pad": pad, "inside": inside, "linked": linked, "outside": outside}


def make_plain_repo(root: Path) -> Path:
    """A repository with one empty commit and its own identity, as `make_repo` sets.

    The identity is in the repository's config, not given per command, so a
    later commit in it (or in a worktree of it) works on a host with no global
    identity, such as the CI runner, where git exits 128.
    """
    root.mkdir(parents=True)
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    git(root, "commit", "-q", "--allow-empty", "-m", "root")
    return root


#: Git writes in or around the scratch repository, `{inside}` its absolute
#: path, each with the reason it asks.
G3_ASKED: dict[str, str] = {
    # G3b's former passes.
    "git -C {inside} config user.name x": CONFIG,
    "git -C {inside} config core.hooksPath x": CONFIG,
    "git -C {inside} remote add o /x": REMOTE,
    "git -C {inside} remote set-url o /y": REMOTE,
    "git -C {inside} update-ref refs/heads/y HEAD": "update-ref moves a ref directly",
    "git -C {inside} symbolic-ref HEAD refs/heads/y": "symbolic-ref rewrites a symbolic ref",
    "git -C {inside} fetch . HEAD:refs/remotes/origin/master": FETCH_WRITE,
    "git -C {inside} replace HEAD HEAD~1": REPLACE,
    "git config --file {inside}/.git/config user.name x": CONFIG,
    "cd {pad} && git -C {inside} config user.name x": CONFIG,
    # Review round 6, finding 1: an abbreviated `--global`, and a glued `-f`
    # naming a config outside the scratchpad.
    "git -C {inside} config --glo user.name x": CONFIG,
    "git -C {inside} config -f{outside}/.git/config user.name x": CONFIG,
    # Asked before option C too: something sends git outside the scratchpad.
    "git -C {linked} config user.name x": CONFIG,
    "git -C {linked} remote add o /x": REMOTE,
    "git -C {linked} update-ref refs/heads/y HEAD": "update-ref moves a ref directly",
    "git -C {outside} config user.name x": CONFIG,
    "GIT_DIR=/x git -C {inside} config user.name x": CONFIG,
    "GIT_CONFIG_GLOBAL=/x git -C {inside} config user.name x": CONFIG,
    "export GIT_DIR={outside}/.git; git -C {inside} config user.name x": CONFIG,
    "git -C {inside} --git-dir={outside}/.git config user.name x": CONFIG,
    "git -C {inside} --work-tree={outside} config user.name x": CONFIG,
    "git -C {inside} config --global user.name x": CONFIG,
    "git -C {inside} config --system user.name x": CONFIG,
    "git config user.name x": CONFIG,
    "git -C r config user.name x": CONFIG,  # relative: the guard does not track cd
    "git -C {inside} -C {outside} config user.name x": CONFIG,
    "git config --file {outside}/.git/config user.name x": CONFIG,
    "git -C {inside} push {pad}/bare master": PUSH,
    "git -C {inside} push file://{pad}/bare master": PUSH,
}


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("template", G3_ASKED)
def test_a_git_write_in_a_scratch_repository_asks(
    repo: Path, scratch: dict[str, Path], mode: str, template: str
) -> None:
    set_mode(repo, mode)
    command = template.format(**scratch)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert G3_ASKED[template] in reason
    assert len(queue_lines(repo)) == (1 if mode == "on" else 0)


# ---------------------------------------------------------------- h16 G6
#
# §2 G6 (review round 6, finding 4): gh's group and verb are read after its
# options, so `-R`/`--repo` (with its value, separate or glued) anywhere before
# them hides neither; `gh pr new` and `gh repo new` are the aliases of create.

RELEASE = "gh publishes or alters the repo"

G6_ASKED: dict[str, str] = {
    "gh -R o/r pr merge 12": GH_PR,
    "gh --repo o/r pr merge 12": GH_PR,
    "gh --repo=o/r pr merge 12": GH_PR,
    "gh -Ro/r pr merge 12": GH_PR,
    "gh pr -R o/r merge 12": GH_PR,
    "gh pr new": GH_PR,
    "gh -R o/r pr new": GH_PR,
    'gh -R o/r pr merge 12 "': GH_PR,  # unreadable: the text rule
    "gh repo new x": RELEASE,
    "gh -R o/r release create v1": RELEASE,
    "gh -R o/r pm 12": UNKNOWN,  # G2, judged on the word after gh's options
}

G6_PASSED = (
    "gh -R o/r pr view 12",
    "gh pr -R o/r view 12",
    "gh --version",
)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G6_ASKED)
def test_a_gh_write_behind_its_repo_option_or_alias_asks(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert G6_ASKED[command] in reason
    if mode == "on":
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_push", command)
    else:
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G6_PASSED)
def test_a_gh_read_behind_its_repo_option_is_silent(repo: Path, mode: str, command: str) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h16 G1, round 7
#
# §2 G1, amendment after review round 7: an `include.`/`includeIf.` key (a
# whole file of config), `core.sshCommand` (the transport for an ssh
# `origin`) and `fetch.bundleURI` (unpacked into `refs/bundles/*`) also steer
# a fetch into named refs, as do `GIT_SSH*`, `HOME` and `XDG_CONFIG_HOME` set
# on the line, which the guard reads from the line's text.

G1_ROUND7_ASKED = (
    "git -c include.path=/x/f fetch s",
    "git -c INCLUDE.PATH=/x/f fetch s",
    "git -c includeIf.onbranch:master.path=/x/f fetch s",
    "git -c includeif.gitdir:/x/.path=/x/f pull s",
    "git --config-env=include.path=F fetch s",
    "git --config-env include.path=F remote update s",
    "git -c core.sshCommand=/x/p fetch origin",
    "git -c fetch.bundleURI=file:///x/b fetch origin",
    "HOME=/x git fetch origin",
    "XDG_CONFIG_HOME=/x git fetch origin",
    "env HOME=/x git fetch origin",
    "GIT_SSH_COMMAND=/x/p git fetch origin",
    "GIT_SSH=/x/p git pull origin",
)

#: Keys the design passes (`fetch.prune` deletes stale tracking refs as
#: `--prune` does; `branch.<name>.remote` chooses what a plain pull takes, as
#: the command line does), a variable that only ends in HOME, and the asked
#: keys on a command that fetches nothing.
G1_ROUND7_PASSED = (
    "git -c fetch.prune=true fetch origin",
    "git -c branch.master.remote=origin pull",
    "JAVA_HOME=/x git fetch origin",
    "git -c include.path=/x/f status",
    "git -c core.sshCommand=/x/p log -1",
)


# ---------------------------------------------------------------- h16 G6, round 7
#
# §2 G6, amendment after review round 7 (Ola's ruling): glued and clustered
# options to `gh api` (a run of `-i` before a value-taking flag, `-fquery=…`)
# and to `curl` (a short cluster holding `d`, `F` or `T`, or `X` with the
# method; `--form-string`, `--expand-*`) are read as the flags they are.

G6_ROUND7_ASKED: dict[str, str] = {
    "gh api graphql -fquery=x": GH_API,
    "gh api graphql -Fquery=@f": GH_API,
    "gh api -ifquery=x graphql": GH_API,
    "gh api -iFquery=x graphql": GH_API,
    "gh api -iXPUT repos/o/r/pulls/1/merge": GH_API,
    "gh api -iiXPUT repos/o/r/pulls/1/merge": GH_API,
    "gh -R o/r api -fquery=x graphql": GH_API,
    "curl -d@f https://api.github.com/graphql": CURL,
    "curl -sd x https://api.github.com/x": CURL,
    "curl -Tf https://api.github.com/x": CURL,
    "curl -sTf https://api.github.com/x": CURL,
    "curl -sFa=b https://api.github.com/x": CURL,
    "curl --form-string a=b https://api.github.com/x": CURL,
    "curl --expand-data x https://api.github.com/x": CURL,
    "curl -sXPUT https://api.github.com/repos/o/r/pulls/1/merge": CURL,
    "curl -sX PUT https://api.github.com/repos/o/r/pulls/1/merge": CURL,
    "curl --expand-request PUT https://api.github.com/x": CURL,
    # `-X=PUT` reads as the method `=PUT`, not GET, so it asks, as it should.
    "gh api -X=PUT repos/o/r/pulls/1/merge": GH_API,
    # Pinned false positives: a data letter inside a glued value (`d` in
    # `data.json`, `T` in `Type`) asks; so does `-X=GET` with a field.
    "curl -o/tmp/data.json https://github.com/x": CURL,
    "curl -HContent-Type:x https://github.com/x": CURL,
    "gh api -X=GET search/issues -f q=x": GH_API,
}

G6_ROUND7_PASSED = (
    "curl -fsSL https://github.com/x",
    "curl -sI https://github.com/x",
    "curl -s -o out https://api.github.com/x",
    # The `T` of `GET` is the method's value, not curl's `-T`: a cluster's
    # data letters are read before its `X` only.
    "curl -sXGET https://api.github.com/x",
    "curl -sX GET https://api.github.com/x",
    "gh api repos/x",
    "gh api -i repos/x",
    "gh api --paginate repos/o/r/pulls",
    "gh api -q .name repos/x",
    "gh api -iXGET repos/x",
    "gh api -X GET search/issues -f q=x",
)


# ---------------------------------------------------------------- h16 G7
#
# §2 G7 (Ola's ruling on §7 question 4): a git or gh command run by another
# program is judged too: (a) the words from any later bare `git`/`gh` word,
# without G2's unknown-command reason; (b) a later `sh`/`bash`/`zsh` word's
# command line, parsed and judged in full; (c) for `watch`, `parallel` and
# `flock` only, each later word holding whitespace, parsed and judged in full.

UPDATE_REF = "update-ref moves a ref directly"

G7_ASKED: dict[str, str] = {
    # (a) a later bare git word
    r"find . -maxdepth 0 -exec git push origin HEAD \;": PUSH,
    "caffeinate -i git push": PUSH,
    "stdbuf -o0 git push": PUSH,
    "watch -n1 git push": PUSH,
    "flock /tmp/l git push": PUSH,
    "parallel git push ::: a": PUSH,
    "arch -arm64 git push": PUSH,
    "caffeinate -i /usr/bin/git push": PUSH,
    "nohup caffeinate stdbuf -o0 git push": PUSH,
    # (c) a quoted command line given to watch, parallel or flock
    "watch -n1 'git push'": PUSH,
    "parallel 'git push' ::: a": PUSH,
    "flock /tmp/l -c 'git push'": PUSH,
    "watch 'caffeinate git push'": PUSH,
    # (b) a shell's command line
    "caffeinate sh -c 'git push'": PUSH,
    "caffeinate sh -c 'cd x && git push'": PUSH,
    r"find . -maxdepth 0 -exec sh -c 'git push' \;": PUSH,
    # every other reason applies through (a)
    "caffeinate -i gh pr merge 12": GH_PR,
    "caffeinate -i gh -R o/r pr merge 12": GH_PR,
    "caffeinate -i gh api -X PUT repos/o/r/pulls/1/merge": GH_API,
    "caffeinate -i git -c remote.s.url=/x fetch s": FETCH_WRITE,
    "caffeinate -i git update-ref refs/heads/x HEAD": UPDATE_REF,
    # (b) judges in full, so an alias given to a shell keeps G2's reason
    "caffeinate sh -c 'git p origin'": UNKNOWN,
    # Both paths: a line the parser cannot read is split into words, and those
    # words go through the same three rules.
    'caffeinate -i git update-ref refs/heads/x HEAD "': UPDATE_REF,
    # Pinned false positives: the words `git push` as arguments ask.
    "echo git push": PUSH,
    "grep git push file": PUSH,
    "man git push": PUSH,
    "parallel 'echo git push' ::: a": PUSH,
}

G7_PASSED = (
    # A later git or gh word that is only an argument: a bare tail never
    # gives G2's reason, and quoted words are read only for the three runners.
    "grep git file",
    "grep -rn git .",
    "grep -c gh tools/guard.py",
    "rg -n gh tools/",
    "echo gh",
    "echo git",
    "which git gh",
    "brew upgrade git gh",
    "git grep -n git -- tools",
    "git log --author git",
    'git commit -m "git push is guarded"',
    'grep -rn "git push" .claude/',
    "rg 'gh pr merge' tools/",
    "grep -rn sh .",
    "grep bash -c x",
    # Reads through a runner.
    "caffeinate -i git status",
    "caffeinate -i gh pr view 12",
    "watch -n5 'gh pr checks 185'",
    "watch -n5 gh pr checks 185",
    # Pinned passes, the residual §6 keeps: an alias through a bare tail. An
    # interpreter program that runs git asks since h25 (test_h25_* below).
    "caffeinate -i git p origin",
    "caffeinate -i gh pm 12",
)


ROUND7_ASKED: dict[str, str] = {
    **dict.fromkeys(G1_ROUND7_ASKED, FETCH_WRITE),
    **G6_ROUND7_ASKED,
    **G7_ASKED,
}


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", ROUND7_ASKED)
def test_an_override_a_glued_forge_write_or_a_run_command_asks(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert ROUND7_ASKED[command] in reason
    if mode == "on":
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_push", command)
    else:
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", [*G1_ROUND7_PASSED, *G6_ROUND7_PASSED, *G7_PASSED])
def test_a_harmless_key_a_read_or_a_word_that_only_names_git_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h16 G7, round 8
#
# §2 G7, amended after review round 8 (Ola's ruling of 2026-10-06): rule (b)
# takes every shell in one set, `shell_scan.SHELLS`, which `guard_push.py`
# reads rather than copies; and under `parallel` a later git or gh word is not
# a bare tail, so it keeps G2's unknown-command reason, since `parallel` builds
# the command it runs from its template and inputs.

shell_scan = Tool("shell_scan")

G7_ROUND8_ASKED: dict[str, str] = {
    # (b) any shell, as the first word or through a runner
    "dash -c 'git push'": PUSH,
    "ksh -c 'git push'": PUSH,
    "csh -c 'git push'": PUSH,
    "tcsh -c 'git push'": PUSH,
    "fish -c 'git push'": PUSH,
    "/bin/ksh -c 'cd x && git push'": PUSH,
    "caffeinate dash -c 'git push'": PUSH,
    r"find . -exec tcsh -c 'git push' \;": PUSH,
    "watch \"dash -c 'git push'\"": PUSH,
    # (a) under parallel: the tail as written is not what runs
    "parallel git ::: push": UNKNOWN,
    "parallel ::: git ::: push": UNKNOWN,
    "parallel git {} ::: push": UNKNOWN,
    "parallel -j2 git ::: push": UNKNOWN,
    # Pinned false positives: harmless lines under parallel ask all the same.
    "parallel git ::: status": UNKNOWN,
    "parallel grep git ::: a": UNKNOWN,
}

#: Controls. `parallel git push ::: a` still asks with PUSH: it is a row of G7_ASKED.
G7_ROUND8_PASSED = (
    "grep -rn dash .",
    "grep ksh file",
    "echo tcsh",
    "which dash ksh",
    "ls /bin/*sh",
    "man csh",
    "dash -c 'git status'",
    "parallel 'git status' ::: a",
)

#: The shells the design names: those `/etc/shells` lists on Ola's Mac, and the
#: common shells taking `-c` that are not installed there.
ON_OLAS_MAC = {"bash", "csh", "dash", "ksh", "sh", "tcsh", "zsh"}
NAMED_SHELLS = ON_OLAS_MAC | {"rbash", "fish", "mksh", "ash", "yash"}


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G7_ROUND8_ASKED)
def test_a_command_run_by_any_shell_or_under_parallel_asks(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert G7_ROUND8_ASKED[command] in reason
    if mode == "on":
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_push", command)
    else:
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G7_ROUND8_PASSED)
def test_a_word_that_only_names_a_shell_or_a_read_under_one_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


def test_the_shell_set_holds_every_shell_the_design_names() -> None:
    assert NAMED_SHELLS - set(shell_scan.SHELLS) == set()


def test_every_shell_takes_its_program_with_c() -> None:
    flags = {name: shell_scan.PROGRAM_FLAG.get(name) for name in shell_scan.SHELLS}
    assert {name: flag for name, flag in flags.items() if flag != "c"} == {}


def test_the_shell_set_is_exactly_the_shells_the_design_names() -> None:
    """A fixed list, not this host's `/etc/shells`: what a machine has installed
    (PowerShell on GitHub's Ubuntu runner) must not decide the result. A shell
    added to or dropped from the set is a design change, made red first here.
    PowerShell is a named known gap, not in the set."""
    assert set(shell_scan.SHELLS) == NAMED_SHELLS


def test_guard_push_holds_no_shell_list_of_its_own(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "path", [*sys.path])  # the hook appends tools/ on import
    hook = load_guard_push()
    assert hook.shell_scan is not None, "guard_push could not import tools/shell_scan.py"
    if hasattr(hook, "SHELLS"):
        assert hook.SHELLS is hook.shell_scan.SHELLS


def load_guard_push() -> ModuleType:
    """Import `.claude/hooks/guard_push.py` from the real checkout, for its module state."""
    path = REAL / ".claude" / "hooks" / "guard_push.py"
    spec = importlib.util.spec_from_file_location("h16_guard_push", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------- h16 G7, round 9
#
# §2 G7, amended after review round 9 (Ola's option A, and his default yes on
# §7 question 5): `runs` carries two flags down the chain of programs, `under`
# (this command, or a runner before it, is `parallel`) and `quoted` (it, or a
# runner before it, is `watch`, `parallel` or `flock`). So a runner in front of
# one of the three no longer hides it, and each word is read once per chain,
# so a line of many shell words is judged in linear time.

#: The heads, each with the reason it asks with when it is the first word.
G7_ROUND9_HEADS: dict[str, str] = {
    "parallel git ::: push": UNKNOWN,
    "parallel ::: git ::: push": UNKNOWN,
    "parallel gh ::: pr ::: merge": UNKNOWN,
    "parallel 'git push' ::: a": PUSH,
    "watch 'git push'": PUSH,
    "flock /tmp/l -c 'git push'": PUSH,
}

#: Runners that are not wrappers `shell_scan` strips: red on `a0001d00`.
#: `{}` is where the head goes.
G7_ROUND9_RUNNERS = (
    "caffeinate {}",
    "caffeinate -i {}",
    "stdbuf -o0 {}",
    "arch -arm64 {}",
    r"find . -maxdepth 0 -exec {} \;",
    "env FOO=1 caffeinate {}",
)

#: The arguments a wrapper needs before the command it runs; a wrapper not
#: named here takes the command at once.
WRAPPER_PREFIX = {"timeout": "timeout 5", "script": "script -q /dev/null"}


def wrapper_runners() -> list[str]:
    """One runner per name in `shell_scan.WRAPPERS`, read from the module, not copied.

    Read at collection, so a missing or broken module gives no rows here rather
    than a collection error; `test_the_matrix_holds_every_wrapper` then fails
    in its body, naming the cause.
    """
    if not (REAL / "tools" / "shell_scan.py").exists():
        return []
    try:
        names = sorted(load_tool("shell_scan").WRAPPERS)
    except Exception:  # judged in the test body, not at collection
        return []
    return [WRAPPER_PREFIX.get(name, name) + " {}" for name in names]


#: The runner matrix: each runner before each head asks with the head's reason.
#: The wrapper rows ask on `a0001d00` too, and guard the stripping.
G7_ROUND9_MATRIX: dict[str, str] = {
    runner.format(head): why
    for runner in (*G7_ROUND9_RUNNERS, *wrapper_runners())
    for head, why in G7_ROUND9_HEADS.items()
}

G7_ROUND9_ASKED: dict[str, str] = {
    **G7_ROUND9_MATRIX,
    # A runner between runners: the flags pass down the whole chain.
    "watch parallel git ::: push": UNKNOWN,
    "caffeinate watch caffeinate parallel git ::: push": UNKNOWN,
    # Asks on `a0001d00`; keeps (c) alive through a shell under the new flags.
    "parallel sh -c {} ::: 'git push'": PUSH,
    # Pinned false positives: a runner's name as an argument starts a chain.
    "grep parallel git ::: x": UNKNOWN,
    "echo parallel gh ::: pr": UNKNOWN,
    "caffeinate parallel git ::: status": UNKNOWN,
    "grep watch 'git push' file": PUSH,
}

#: Controls: reads through a runner before one of the three, and quoted words
#: that only name a runner. Round 7's and round 8's controls stay as they are.
G7_ROUND9_PASSED = (
    "caffeinate parallel 'git status' ::: a",
    "stdbuf -o0 flock /tmp/l -c 'git status'",
    "caffeinate watch -n5 'gh pr checks 185'",
    "caffeinate watch -n5 git status",
    "watch sh -c 'git commit -m \"git push\"'",
    "git log --grep 'watch git push'",
)


def test_the_matrix_holds_every_wrapper() -> None:
    names = set(shell_scan.WRAPPERS)
    assert names, "shell_scan.WRAPPERS is empty: the wrapper rows would measure nothing"
    runners = {runner.split()[0] for runner in wrapper_runners()}
    assert runners == names
    assert len(G7_ROUND9_MATRIX) == len(G7_ROUND9_HEADS) * (len(G7_ROUND9_RUNNERS) + len(names))


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G7_ROUND9_ASKED)
def test_a_runner_in_front_of_parallel_watch_or_flock_asks_with_the_heads_reason(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert G7_ROUND9_ASKED[command] in reason
    if mode == "on":
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_push", command)
    else:
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G7_ROUND9_PASSED)
def test_a_read_through_a_runner_before_parallel_watch_or_flock_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# The slow line (§7 question 5). Rule (b) as on `a0001d00` parses the rest of
# the line at every shell word, so the cost doubles per word: through the hook,
# 23 `sh` words before `-c 'git push'` took 6.4 s in one run and 50.6 s in
# another, and 40 did not finish in 120 s. Claude Code lets the command through
# when a hook times out (600 s), so a slow guard is an open guard.
#
# The bound, measured from outside the hook process: 5 s per line, for lines
# of up to 400 runner or shell words (the largest checked here). The design's
# prototype of the linear `runs` answered each row in under 0.2 s through the
# hook, Python's start included, so 5 s leaves room for a loaded CI runner and
# still fails the doubling cost by orders of magnitude at 40 words.
SLOW_LINE_BOUND_S = 5.0

#: {line: the reason it asks with}; each has 30 or 40 runner or shell words.
SLOW_LINES_ASKED: dict[str, str] = {
    "sh " * 40 + "-c 'git push'": PUSH,
    "watch " * 40 + "'git push'": PUSH,
    "caffeinate sh -c " * 30 + "'git push'": PUSH,
}

#: 400 shell words before a read: judged without a deny from Python's
#: recursion limit (1000 frames), which the design puts at one frame per word.
SLOW_LINE_PASSED = "bash " * 400 + "-c 'git status'"


def timed_hook(repo: Path, command: str) -> subprocess.CompletedProcess[str]:
    """Run the copied push guard on `command`, failing if it takes longer than the bound."""
    script = repo / GUARD_PUSH
    try:
        return subprocess.run(
            [sys.executable, str(script)],
            input=json.dumps(bash_event(repo, command)),
            capture_output=True,
            text=True,
            cwd=repo,
            env=clean_env(),
            timeout=SLOW_LINE_BOUND_S,
            check=False,
        )
    except subprocess.TimeoutExpired:
        pytest.fail(
            f"guard_push took longer than {SLOW_LINE_BOUND_S} s on a line of "
            f"{len(command.split())} words: past the hook's time limit the command runs unasked"
        )


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize(
    "command", SLOW_LINES_ASKED, ids=["40-sh", "40-watch", "30-caffeinate-sh-c"]
)
def test_a_line_of_many_runner_or_shell_words_is_judged_within_the_bound(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(timed_hook(repo, command))
    assert found is not None, "the line passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert SLOW_LINES_ASKED[command] in reason


@pytest.mark.parametrize("mode", ["off", "on"])
def test_four_hundred_shell_words_before_a_read_pass_within_the_bound(
    repo: Path, mode: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(timed_hook(repo, SLOW_LINE_PASSED)) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h20 §4.3
#
# docs/increments/h20-away-push-grants.md §3.2 and §3.5: while away, a branch
# Ola granted may be pushed (by sha), opened as a pull request and enqueued
# once per reviewed head. Each case runs on `Granted`'s fixture repository.
# Pinned here beyond the design's words: each line of `grants-used.jsonl`
# carries the keys §3.5 lists, and test 15's refusal carries the grant clause
# of §3.5 and the word "enqueued".

GRANT_CLAUSE = "granted branch feat, but"


@pytest.fixture
def granted(tmp_path: Path) -> Granted:
    return Granted(tmp_path / "repo")


def granted_acts(head: str) -> list[str]:
    """The three forms of §3.2 for `feat` at `head`, in the order a night uses them."""
    return [
        f"git push origin {head}:refs/heads/feat",
        "gh pr create --head feat --title 't' --body-file f",
        f"gh pr merge feat --match-head-commit {head} --merge",
    ]


def guard(repo: Path, command: str, env: dict[str, str] | None = None) -> tuple[str, str] | None:
    if env is None:
        return pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    result = subprocess.run(
        [sys.executable, str(repo / GUARD_PUSH)], input=json.dumps(bash_event(repo, command)),
        capture_output=True, text=True, cwd=repo, env=env, timeout=60, check=False,
    )  # fmt: skip
    return pretool_decision(result)


def denied(repo: Path, command: str, env: dict[str, str] | None = None) -> str:
    """The `why` of the queue line a refused command left; fails if it was not denied."""
    before = len(queue_lines(repo))
    found = guard(repo, command, env)
    assert found is not None and found[0] == "deny", f"{command!r}: {found}"
    lines = queue_lines(repo)
    assert len(lines) == before + 1, "the refusal was not queued"
    why: str = lines[-1]["why"]
    return why


def test_10_a_granted_reviewed_branch_is_pushed_opened_and_enqueued(granted: Granted) -> None:
    head = granted.head
    for command in granted_acts(head):
        assert guard(granted.repo, command) is None, command
    used = jsonl(grants_used_path(granted.repo))
    assert len(used) == 3
    assert len({line["act"] for line in used}) == 3
    for line in used:
        assert set(line) == {"at", "window", "branch", "act", "head", "review", "remote"}
        assert (line["branch"], line["head"], line["remote"]) == ("feat", head, "origin")
        assert line["window"] == iso(granted.since)
        assert line["review"][0] == granted.x
    assert queue_lines(granted.repo) == []


def test_11_with_the_mode_off_a_granted_act_asks_as_today(granted: Granted) -> None:
    flag_path(granted.repo).unlink()
    for command in granted_acts(granted.head):
        found = guard(granted.repo, command)
        assert found is not None and found[0] == "ask", command
    assert not grants_used_path(granted.repo).exists()


def test_12_with_the_flag_broken_a_granted_act_is_denied(granted: Granted) -> None:
    flag_path(granted.repo).write_text("{not json")
    for command in granted_acts(granted.head):
        denied(granted.repo, command)
    assert not grants_used_path(granted.repo).exists()


def test_13_a_branch_not_granted_is_denied_without_the_grant_clause(granted: Granted) -> None:
    flag_granting(granted.repo, ["other"])
    for command in granted_acts(granted.head):
        assert "granted branch" not in denied(granted.repo, command)
    assert not grants_used_path(granted.repo).exists()


@pytest.mark.parametrize(
    "wrap",
    ["env A=1 {}", "{} 2>&1", "{} && echo done", "sh -c '{}'", "timeout 5 {}"],
    ids=["env", "redirection", "and-list", "sh-c", "timeout"],
)
def test_14_a_granted_act_inside_any_wrapper_is_denied(granted: Granted, wrap: str) -> None:
    denied(granted.repo, wrap.format(f"git push origin {granted.head}:refs/heads/feat"))
    assert not grants_used_path(granted.repo).exists()


def test_15_one_enqueue_per_reviewed_head(granted: Granted) -> None:
    merge = f"gh pr merge feat --match-head-commit {granted.head}"
    assert guard(granted.repo, merge) is None
    why = denied(granted.repo, merge)
    assert GRANT_CLAUSE in why and "enqueued" in why, why
    granted.on_master("src/y.py", "y = 1\n")
    second = granted.merge("refs/remotes/origin/master")
    again = f"gh pr merge feat --match-head-commit {second}"
    assert guard(granted.repo, again) is None
    denied(granted.repo, again)
    assert [line["head"] for line in jsonl(grants_used_path(granted.repo))] == [
        granted.rec, second
    ]  # fmt: skip


@pytest.mark.parametrize("act", ["push", "enqueue"])
def test_16_a_sha_other_than_the_head_is_denied(granted: Granted, act: str) -> None:
    command = {
        "push": f"git push origin {granted.x}:refs/heads/feat",
        "enqueue": f"gh pr merge feat --match-head-commit {granted.x} --merge",
    }[act]
    denied(granted.repo, command)
    assert not grants_used_path(granted.repo).exists()


def test_17_a_push_run_in_another_repository_is_denied(granted: Granted, tmp_path: Path) -> None:
    other = tmp_path / "other"
    git(tmp_path, "clone", "-q", str(granted.repo), str(other))
    git(other, "update-ref", "refs/heads/feat", granted.head)
    git(other, "update-ref", "refs/remotes/origin/master", granted.base)
    denied(granted.repo, f"git -C {other} push origin {granted.head}:refs/heads/feat")
    assert not grants_used_path(granted.repo).exists()


@pytest.mark.parametrize("directory", ["~/other", "=other"])
def test_17_a_directory_the_shell_would_expand_is_denied(
    granted: Granted, tmp_path: Path, directory: str
) -> None:
    """`~/other` (through HOME) and `=other` (through PATH) reach the guard's own repository."""
    home, bin_dir = tmp_path / "home", tmp_path / "bin"
    home.mkdir()
    bin_dir.mkdir()
    (home / "other").symlink_to(granted.repo)
    (bin_dir / "other").symlink_to(granted.repo)
    env = {**clean_env(), "HOME": str(home), "PATH": f"{bin_dir}{os.pathsep}{clean_env()['PATH']}"}
    command = f"git -C {directory} push origin {granted.head}:refs/heads/feat"
    denied(granted.repo, command, env)
    assert not grants_used_path(granted.repo).exists()


def test_18_a_grants_log_that_cannot_be_written_denies(granted: Granted) -> None:
    grants_used_path(granted.repo).mkdir(parents=True)
    denied(granted.repo, f"git push origin {granted.head}:refs/heads/feat")


# ---------------------------------------------------------------- h24 item 2
#
# docs/increments/h24-harness-fixes.md §3 item 2: a granted act written in a
# line that is not one bare command is refused, as before, but the refusal
# says so (NOT_BARE), so the session can rerun the command alone. The words
# only choose the refusal's text; the rerun is judged by the grant as today.

NOT_BARE = (
    "the grant takes only one bare command, every word written out "
    "(no variable, no &&, ;, pipe, redirect or heredoc); run it alone"
)

#: The 00:22:54 line of 2026-10-18 (body cut to one line), for `feat`.
PR_BODY_THEN_CREATE = (
    "S=/abs/scratchpad && cat > $S/pr68.md <<'EOF'\n"
    "One line of body.\n"
    "EOF\n"
    "gh pr create --head feat --title 'T' --body-file $S/pr68.md"
)
BARE_CREATE = "gh pr create --head feat --title 'T' --body-file /abs/pr.md"


def test_h24_2_a_granted_act_in_a_line_of_several_commands_says_not_bare(
    granted: Granted,
) -> None:
    why = denied(granted.repo, PR_BODY_THEN_CREATE)
    assert f"{GRANT_CLAUSE} {NOT_BARE}" in why, why
    assert not grants_used_path(granted.repo).exists()


def test_h24_2_the_same_line_for_a_branch_not_granted_is_the_plain_refusal(
    granted: Granted,
) -> None:
    flag_granting(granted.repo, ["other"])
    why = denied(granted.repo, PR_BODY_THEN_CREATE)
    assert GH_PR in why, why
    assert NOT_BARE not in why and "granted branch" not in why, why


def test_h24_2_a_push_naming_its_sha_by_a_variable_says_not_bare(granted: Granted) -> None:
    why = denied(granted.repo, "git push origin $SHA:refs/heads/feat")
    assert f"{GRANT_CLAUSE} {NOT_BARE}" in why, why
    assert not grants_used_path(granted.repo).exists()


def test_h24_2_the_bare_command_rerun_after_a_not_bare_refusal_passes(granted: Granted) -> None:
    """The rerun the refusal asks for: one bare command, judged by the grant as before."""
    denied(granted.repo, PR_BODY_THEN_CREATE)
    assert guard(granted.repo, BARE_CREATE) is None
    [line] = jsonl(grants_used_path(granted.repo))
    assert (line["branch"], line["act"], line["head"]) == ("feat", "pull request", granted.head)


# ---------------------------------------------------------------- h24 item 5a
#
# docs/increments/h24-harness-fixes.md §3 item 5a: a `gh api graphql` segment
# whose one `query` field reads (starts with `{` or `query`, and holds the
# letters `mutation` and `subscription` nowhere) passes; anything the guard
# cannot read word for word keeps today's rule, where a field means a write.
# `segment_why` runs on the words as the shell would split them, before
# expansion, so `$X` is the word the guard sees.

#: The 23:03:56 read of 2026-10-17.
PR_EDITS = (
    'query={ repository(owner:"expertanalytics", name:"rasputin") '
    "{ pullRequest(number:285) { userContentEdits(first:5) { nodes { editedAt } } } } }"
)
GRAPHQL = ["gh", "api", "graphql"]

GRAPHQL_READS: dict[str, list[str]] = {
    "the-23-03-56-read": [*GRAPHQL, "-f", PR_EDITS],
    "shorthand": [*GRAPHQL, "-f", "query={ viewer { login } }"],
    "named-with-operation-name": [*GRAPHQL, "-f", "query=query A { viewer { login } }",
                                  "-f", "operationName=A"],
    "method-post": [*GRAPHQL, "-X", "POST", "-f", "query={ viewer { login } }"],
    "raw-field-and-jq": [*GRAPHQL, "--raw-field", "query={ viewer { login } }",
                         "--jq", ".data.viewer.login"],
}  # fmt: skip

GRAPHQL_ASKED: dict[str, list[str]] = {
    # The design review's writes: a `#` inside a string is no comment, an
    # expansion hides what gh receives, and a second plain word is unread.
    "hash-in-a-string-then-a-mutation": [
        *GRAPHQL, "-f", 'query=query A { a(x:"#") } mutation B { addComment(input:{}) { c } }',
        "-f", "operationName=B"],
    "an-expansion-inside-the-query": [*GRAPHQL, "-f", "query=query A {a} $M",
                                      "-f", "operationName=B"],
    "an-expansion-as-a-word": [*GRAPHQL, "$X", "-f", "query={ a }"],
    "a-second-plain-word": [*GRAPHQL, "repos/a/b/issues", "-f", "query={ a }"],
    # Each a write, or a query the guard cannot read.
    "variables": [*GRAPHQL, "-f", "query=query($n:Int!){ a }", "-F", "n=285"],
    "mutation": [*GRAPHQL, "-f", "query=mutation { addComment(input:{}) { c } }"],
    "query-then-mutation": [*GRAPHQL, "-f", "query=query { a } mutation { b }"],
    "subscription": [*GRAPHQL, "-f", "query=subscription { a }"],
    "query-from-a-file": [*GRAPHQL, "-F", "query=@q.graphql"],
    "input-file": [*GRAPHQL, "--input", "q.json"],
    "method-put": [*GRAPHQL, "-X", "PUT", "-f", "query={ a }"],
    "two-queries": [*GRAPHQL, "-f", "query={ a }", "-f", "query={ b }"],
    "a-variable-beside-the-query": [*GRAPHQL, "-f", "query={ a }", "-f", "n=1"],
    "method-override-header": [*GRAPHQL, "-H", "X-HTTP-Method-Override: PUT",
                               "-f", "query={ a }"],
    "paginate": [*GRAPHQL, "--paginate", "-f", "query={ a }"],
    "a-flag-with-its-value-glued": [*GRAPHQL, "--raw-field=query={ a }"],
    "a-backtick": [*GRAPHQL, "-f", "query={ a(x:`id`) }"],
    "a-leading-comment": [*GRAPHQL, "-f", "query=# c\n{ a }"],
    "rest-with-a-field": ["gh", "api", "repos/a/b/issues", "-f", "title=x"],
}  # fmt: skip


@pytest.fixture
def guard_push_hook(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    monkeypatch.setattr(sys, "path", [*sys.path])  # the hook appends tools/ on import
    return load_guard_push()


@pytest.mark.parametrize("words", GRAPHQL_READS.values(), ids=GRAPHQL_READS.keys())
def test_h24_5a_a_graphql_read_passes(guard_push_hook: ModuleType, words: list[str]) -> None:
    assert guard_push_hook.segment_why(words) is None


@pytest.mark.parametrize("words", GRAPHQL_ASKED.values(), ids=GRAPHQL_ASKED.keys())
def test_h24_5a_a_graphql_write_or_unreadable_query_asks(
    guard_push_hook: ModuleType, words: list[str]
) -> None:
    assert guard_push_hook.segment_why(words) == GH_API


# ------------------------------------------------- h24 item 5a, code review 1
#
# The guard reads the words before the shell expands them. It refuses `$` and
# a backtick, but brace expansion and a glob also change the words gh gets:
# each line below passes as a read today, yet the shell hands gh a mutation,
# a `-X PUT`, or the name of a file holding a mutation. The fix proposed in
# review: a read only when every flag value is plain ([A-Za-z0-9_.=-]) or
# single-quoted in the command text, so these run through the whole hook.

#: A mutation the brace expansion turns into a second `-f query=` field.
BRACED_MUTATION = (
    '\'query=mutation { addComment(input:{subjectId:"X",body:"y"}) { clientMutationId } }\''
)
EXPANDED_WRITES: dict[str, str] = {
    "braces-give-a-second-query": (
        f"gh api graphql -f query='{{ a }}' -q {{.,-f,{BRACED_MUTATION}}}"
    ),
    "braces-give-a-put": "gh api graphql -f query='{ a }' -q {x,-X,PUT}",
    "a-glob-names-a-mutation-file": "gh api graphql -f query={a}*",
}
#: The file the glob of "a-glob-names-a-mutation-file" matches in the working folder.
GLOB_TARGET = "query={a}mutation{b}"

GRAPHQL_READ_LINES: dict[str, str] = {
    "single-quoted-query": "gh api graphql -f query='{ viewer { login } }'",
    "whole-field-single-quoted": "gh api graphql -f 'query={ viewer { login } }'",
    "plain-values": (
        "gh api graphql -X POST -f 'query=query A { viewer { login } }' "
        "-f operationName=A --jq .data.viewer.login"
    ),
}


@pytest.mark.parametrize("command", EXPANDED_WRITES.values(), ids=EXPANDED_WRITES.keys())
def test_h24_5a_a_graphql_line_the_shell_expands_asks(repo: Path, command: str) -> None:
    set_mode(repo, "off")
    (repo / GLOB_TARGET).write_text("")
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert GH_API in reason, reason


@pytest.mark.parametrize("command", GRAPHQL_READ_LINES.values(), ids=GRAPHQL_READ_LINES.keys())
def test_h24_5a_a_plain_or_single_quoted_graphql_read_passes(repo: Path, command: str) -> None:
    set_mode(repo, "off")
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None


# ------------------------------------------------- h24 item 5a, design review 3
#
# `segment_why` gets the whole line as `text`, as `main` passes it, and the
# words the scanner gives for the gh segment. zsh expands a word starting
# with `=` to a path, also when written `''=ls`. A line holding a redirect
# expands nothing gh sees but still fails the whole-line pattern, so it asks;
# it must ask fast, since a hook that times out lets the call through.


def gh_segment(hook: ModuleType, line: str) -> list[str]:
    """The words of the line's one gh simple command, as `main` hands them on."""
    simples = hook.shell_scan.parse(line)
    assert simples is not None, line
    [argv] = [s.argv for s in simples if s.argv[:1] == ["gh"]]
    words: list[str] = argv
    return words


def test_h24_5a_a_word_starting_with_equals_after_quote_removal_asks(
    guard_push_hook: ModuleType,
) -> None:
    line = "gh api graphql -f query='{ a }' -q ''=ls"
    words = gh_segment(guard_push_hook, line)
    assert words[-1] == "=ls", words
    assert guard_push_hook.segment_why(words, line) == GH_API


#: Design review round 3's failing line: 72 characters, ending in a redirect.
#: On a backtracking pattern a 36-character failing line took 5.8 s.
REDIRECT_LINE = "gh api graphql -f query='{ a }' --jq .data.repository.id > /tmp/out.json"


def test_h24_5a_a_long_failing_line_asks_in_under_a_second(guard_push_hook: ModuleType) -> None:
    assert len(REDIRECT_LINE) == 72
    words = gh_segment(guard_push_hook, REDIRECT_LINE)
    start = time.perf_counter()
    why = guard_push_hook.segment_why(words, REDIRECT_LINE)
    elapsed = time.perf_counter() - start
    assert why == GH_API
    assert elapsed < 1.0, f"segment_why took {elapsed:.3f} s on a 72-character line"


# ---------------------------------------------------------------- h24 item 15
#
# A here-string fed to a shell is the shell's program, parsed as commands, so a
# publishing act inside it is judged like one on the line. zsh joins every
# stdin redirect of one command, so a second here-string counts too.

#: {command: the reason it must name}.
HERE_STRING_ASKED: dict[str, str] = {
    "bash <<< 'git push origin x'": PUSH,
    'sh <<< "gh pr merge 1"': GH_PR,
    "env bash <<< 'git push origin x'": PUSH,
    "timeout 5 bash <<< 'git push origin x'": PUSH,
    "zsh <<< 'echo a' <<< 'git push origin x'": PUSH,
    # A shell parses every stdin body whatever its -c program: the -c program
    # here is itself a shell reading stdin.
    "bash -c bash <<'EOF'\ngit push origin x\nEOF": PUSH,
    "sh -c 'cat | sh' <<'EOF'\ngit push origin x\nEOF": PUSH,
    "bash -c bash <<< 'git push origin x'": PUSH,
    # A heredoc on a shell that has any argument feeds what the script runs,
    # so the whole line is judged by its text
    # (docs/increments/h25-stdin-routes.md §3.2, line feeds; §3.5).
    "bash script.sh <<'EOF'\ngit push origin x\nEOF": PUSH,
}

HERE_STRING_PASSED = (
    "cat <<< 'git push origin x'",
    "bash <<< 'git status'",
    "python3 -c 'print(1)' <<< 'git push origin x'",
    # An interpreter's inline program keeps a heredoc as data.
    "python3 -c 'print(1)' <<'EOF'\ngit push origin x\nEOF",
)


@pytest.mark.parametrize("command", HERE_STRING_ASKED)
def test_h24_15_a_publishing_act_in_a_shell_here_string_asks(repo: Path, command: str) -> None:
    set_mode(repo, "off")
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert HERE_STRING_ASKED[command] in reason, reason
    assert queue_lines(repo) == []


@pytest.mark.parametrize("command", HERE_STRING_ASKED)
def test_h24_15_a_publishing_act_in_a_shell_here_string_is_queued_while_unattended(
    repo: Path, command: str
) -> None:
    set_mode(repo, "on")
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "deny"
    assert HERE_STRING_ASKED[command] in reason, reason
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_push", command)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", HERE_STRING_PASSED)
def test_h24_15_a_here_string_that_is_data_or_reads_passes(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h25
#
# docs/increments/h25-stdin-routes.md §3.5: an unread runner (fed text the
# scanner did not read) puts the whole line under the PUBLISHES text rules,
# and a program the scanner reads that may start a process (SPAWN) is matched
# against PUBLISHES, strings and comments included.

H25_ASKED = (
    "python3 -c \"import os; os.system('git push origin x')\"",
    "python3 -c \"import subprocess; subprocess.run(['git','push','origin','x'])\"",
    "python3 -c \"import subprocess; subprocess.run(['git', 'push'])\"",
    "perl -e 'system(\"git push origin x\")'",
    "bash -c python3 <<< \"import os; os.system('git push origin x')\"",
    "echo 'git push origin x' | sh",
    "echo 'git push origin x' | bash -s",
    "echo 'git push origin x' | { true; sh; }",
    "sh <(echo 'git push origin x')",
    "python3 -c \"import operator,sys; "
    "operator.methodcaller('system','git push origin x')(sys.modules['posix'])\"",
    "ruby -e '%x(git push origin x)'",
    "perl -e 'open(F, \"|git push origin x\")'",
    "read -r l <<< 'git push origin x'; eval \"$l\"",
)  # fmt: skip

H25_PASSED = (
    "echo 'git push origin x' | cat",
    "echo 'git push origin x' | python3 x.py",
    "python3 -c \"print('git push origin x')\"",
    "git apply x.patch; echo 'git push origin x'",
)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", H25_ASKED)
def test_h25_a_fed_runner_or_a_spawning_program_naming_a_push_asks(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert PUSH in reason, reason
    if mode == "on":
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_push", command)
    else:
        assert reason.startswith("This reaches beyond the working tree: ")
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", H25_PASSED)
def test_h25_a_push_named_as_data_or_printed_is_silent(repo: Path, mode: str, command: str) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# ------------------------------------------------- h25 code review round 1
#
# docs/increments/h25-stdin-routes.md §3.4: a script operand "not under
# /dev/" is judged after `.` and `..` are resolved, so a path that reaches
# /dev/stdin through them runs stdin as `python3 /dev/stdin` does: the
# interpreter is fed, and the whole line goes under the PUBLISHES rules.

#: Script operands that resolve to stdin; the first two are the controls.
H25_STDIN_SCRIPTS = (
    "-",
    "/dev/stdin",
    "/./dev/stdin",
    "../../../../../../../../dev/stdin",
    "/usr/../dev/stdin",
    "/usr/lib/../../dev/fd/0",
)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("script", H25_STDIN_SCRIPTS)
def test_h25_a_script_path_resolving_to_stdin_feeds_the_interpreter(
    repo: Path, mode: str, script: str
) -> None:
    set_mode(repo, mode)
    command = f"echo 'git push origin x' | python3 {script}"
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert PUSH in reason, reason


# §3.5: "A program that names `git push` and no way to start a process (a
# `print`, an edit of a design file) publishes nothing." Python has no
# backtick, `%x` or pipe-in-a-string route to a process, so those three
# SPAWN alternatives are Perl's, Ruby's and a shell's; a Python program
# editing Markdown that names `git push` is caught by them alone today.

H25_PYTHON_MARKDOWN_EDITS: dict[str, str] = {
    "c-table-row-quote-next-to-pipe": (
        "python3 -c 'p = \"docs/x.md\"; s = open(p).read(); "
        "open(p, \"w\").write(s + \"| git push origin x | publishes |\\n\")'"),
    "heredoc-table-row-quote-next-to-pipe": (
        "python3 - <<'EOF'\np = 'docs/x.md'\ns = open(p).read()\n"
        "open(p, 'w').write(s + \"| git push origin x | publishes the branch |\\n\")\nEOF"),
    "heredoc-percent-x-format": (
        "python3 - <<'EOF'\np = 'docs/x.md'\n"
        "open(p, 'a').write('After git push origin x, the head is %x.\\n' % 3054)\nEOF"),
    "heredoc-backtick-and-table": (
        "python3 <<'EOF'\nimport pathlib\np = pathlib.Path('docs/increments/x.md')\n"
        "p.write_text(p.read_text() + '| `git push origin x` | asks |\\n')\nEOF"),
}  # fmt: skip

#: The same Markdown edits in languages where the backtick, `%x` or a quote
#: next to `|` can start a process: they still ask.
H25_OTHER_MARKDOWN_EDITS: dict[str, str] = {
    "ruby-backtick": "ruby -e 'File.write(\"docs/x.md\", \"run `git push origin x`\")'",
    "ruby-percent-x": "ruby -e 'File.write(\"docs/x.md\", %x(git push origin x))'",
    "perl-quote-next-to-pipe": (
        "perl -e 'open(F, \">>docs/x.md\"); print F \"| git push origin x |\\n\"'"),
    "perl-heredoc-table-row": (
        "perl <<'EOF'\nopen(F, '>>docs/x.md');\n"
        "print F \"| git push origin x | publishes |\\n\";\nEOF"),
}  # fmt: skip


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize(
    "command", H25_PYTHON_MARKDOWN_EDITS.values(), ids=H25_PYTHON_MARKDOWN_EDITS.keys()
)
def test_h25_a_python_markdown_edit_naming_a_push_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize(
    "command", H25_OTHER_MARKDOWN_EDITS.values(), ids=H25_OTHER_MARKDOWN_EDITS.keys()
)
def test_h25_the_same_edit_in_perl_or_ruby_still_asks(repo: Path, mode: str, command: str) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert PUSH in reason, reason


# ------------------------------------------------- h25 code review round 2
#
# §3.4, as in test_guard_governance.py: on macOS `/DEV/stdin` and `/Dev/fd/0`
# are stdin (the root is case-insensitive), and a relative script on a line
# that runs cd or pushd is resolved against a directory the guard cannot
# know. Each was run for real on 2026-10-20 (zsh, bash); each feeds the
# interpreter, so the whole line goes under the PUBLISHES rules.

H25_STDIN_ROUTES_ROUND2: dict[str, str] = {
    "case-dev-stdin": "echo 'git push origin x' | python3 /DEV/stdin",
    "case-dev-fd-0": "echo 'git push origin x' | python3 /Dev/fd/0",
    "cd-and": "cd /dev && echo 'git push origin x' | python3 stdin",
    "cd-semicolon-fd": "cd /dev/fd; echo 'git push origin x' | python3 0",
    "cd-variable": "D=/dev; cd \"$D\" && echo 'git push origin x' | python3 stdin",
    "pushd-and": "pushd /dev && echo 'git push origin x' | python3 stdin",
    "pushd-semicolon-fd": "pushd /dev/fd; echo 'git push origin x' | python3 0",
}  # fmt: skip

# §3.5 exempts Python from SPAWN's three marks (backtick, `%x`, a quote next
# to `|`) because Python has no such route. The exemption must follow what
# runs, not the word written: when the line redefines `python3` as a
# function, an alias or through PATH, the body is Ruby's or Perl's and the
# marks start a process. Each was run for real on 2026-10-20 with a harmless
# command in place of the push: the function forms and the PATH prefix (a
# directory whose `python3` is ruby) under zsh and bash; the alias with a
# newline after it under bash with expand_aliases (an alias on the same line
# as its use is not expanded, and is not pinned).

H25_REDEFINED_PYTHON: dict[str, str] = {
    "function-ruby-backtick": "python3() { ruby; }; python3 <<'EOF'\nputs `git push origin x`\nEOF",
    "function-ruby-percent-x": (
        "python3() { ruby; }; python3 <<'EOF'\nputs %x(git push origin x)\nEOF"),
    "function-perl-quote-next-to-pipe": (
        "python3() { perl; }; python3 <<'EOF'\nopen(F, \"|git push origin x\"); close F;\nEOF"),
    "alias-ruby-backtick": "alias python3=ruby\npython3 <<'EOF'\nputs `git push origin x`\nEOF",
    "function-keyword-perl-backtick": (
        "function python3 { perl; }; python3 <<'EOF'\nprint `git push origin x`;\nEOF"),
    "path-prefix-backtick": "PATH=/tmp:$PATH python3 <<'EOF'\nputs `git push origin x`\nEOF",
}  # fmt: skip


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize(
    "command",
    [*H25_STDIN_ROUTES_ROUND2.values(), *H25_REDEFINED_PYTHON.values()],
    ids=[*H25_STDIN_ROUTES_ROUND2.keys(), *H25_REDEFINED_PYTHON.keys()],
)
def test_h25_round2_stdin_route_or_redefined_python_naming_a_push_asks(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert PUSH in reason, reason


# ---------------------------------------------------------------- h27
#
# docs/increments/h27-unpushed-amend.md §3: an amend whose HEAD no
# remote-tracking ref and no tag holds, on a bare one-command line (h29 below
# widens the line), passes
# silently in every mode; anything else asks by day and is refused and queued
# while away. Each fixture pushes into a bare `origin` with git itself, so the
# tracking ref is the one a real push leaves.

HISTORY = "this rewrites history, which is destructive once anything is published"
NO_VERIFY = "--no-verify disables git's own hooks"
#: h29 §3.3's message, without its leading "; ", word for word.
AMEND_FORM = (
    "an amend of a commit no remote holds passes when it comes last, after at most "
    "`cd <absolute path>` and `git add` commands joined by `&&`, or as "
    "`git -C <absolute path> commit --amend …` alone; it may be redone in that form"
)


class Line(NamedTuple):
    """A command and the event's cwd; None means the event has no `cwd` key."""

    command: str
    cwd: Path | None


#: A case: given the pushed repository and tmp_path, it sets up and returns the line.
Case = Callable[[Path, Path], Line]


@pytest.fixture
def pushed(tmp_path: Path) -> Path:
    """`make_repo` with `master` pushed by `git push` into a bare `origin`."""
    repo = make_repo(tmp_path / "repo")
    git(tmp_path, "init", "-q", "--bare", "origin.git")
    git(repo, "remote", "add", "origin", str(tmp_path / "origin.git"))
    git(repo, "push", "-q", "origin", "master")
    # §1: the push itself leaves the tracking ref; checked, not assumed.
    assert git(repo, "rev-parse", "refs/remotes/origin/master") == git(repo, "rev-parse", "HEAD")
    # The exemption needs a line the shell expands nothing in, paths included.
    assert all(c.isalnum() or c in "/._-" for c in str(tmp_path)), tmp_path
    return repo


def local_commit(repo: Path) -> str:
    return commit_file(repo, "notes.txt", "local work\n", "local work, never pushed")


def plain_amend(repo: Path, tmp_path: Path) -> Line:
    local_commit(repo)
    return Line("git commit --amend --no-edit", repo)


def amend_by_path(repo: Path, tmp_path: Path) -> Line:
    """The refused form of 2026-10-19; cwd is not a repository, so `-C` must decide."""
    local_commit(repo)
    return Line(f"git -C {repo} commit --amend --no-edit -a", tmp_path)


def amend_in_worktree(repo: Path, tmp_path: Path) -> Line:
    """cwd is the main checkout, whose HEAD is pushed; the worktree's is not."""
    tree = add_worktree(repo, tmp_path / "tree", "side")
    local_commit(tree)
    return Line(f"git -C {tree} commit --amend --no-edit", repo)


def abbreviated(option: str) -> Case:
    def case(repo: Path, tmp_path: Path) -> Line:
        local_commit(repo)
        return Line(f"git commit {option} --no-edit", repo)

    return case


def not_bare(template: str, cwd: str = "repo") -> Case:
    """An unpushed HEAD on a line other than one bare `git commit --amend`."""

    def case(repo: Path, tmp_path: Path) -> Line:
        local_commit(repo)
        return Line(template.format(repo=repo), repo if cwd == "repo" else tmp_path)

    return case


H27_PASSED = {
    "1-no-C": plain_amend,
    "2-absolute-C-with-a": amend_by_path,
    "3-worktree": amend_in_worktree,
    "4-amen": abbreviated("--amen"),
    "4-am": abbreviated("--am"),
    # Recomputed for h29 §3.3: asked under h27, passes once `cd` may come first.
    "9-cd-first": not_bare("cd {repo} && git commit --amend --no-edit"),
}


def pushed_head(repo: Path, tmp_path: Path) -> Line:
    return Line("git commit --amend --no-edit", repo)


def pushed_head_by_path(repo: Path, tmp_path: Path) -> Line:
    return Line(f"git -C {repo} commit --amend --no-edit -a", tmp_path)


def tagged(repo: Path, tmp_path: Path) -> Line:
    local_commit(repo)
    git(repo, "tag", "t")
    return Line("git commit --amend --no-edit", repo)


def pushed_elsewhere(repo: Path, tmp_path: Path) -> Line:
    local_commit(repo)
    git(repo, "push", "-q", "origin", "HEAD:refs/heads/x")
    assert git(repo, "rev-parse", "refs/remotes/origin/x") == git(repo, "rev-parse", "HEAD")
    return Line("git commit --amend --no-edit", repo)


def pushed_from_worktree(repo: Path, tmp_path: Path) -> Line:
    """§2: worktrees share refs, so a push from the worktree is seen from anywhere."""
    tree = add_worktree(repo, tmp_path / "tree", "side")
    local_commit(tree)
    git(tree, "push", "-q", "origin", "side")
    return Line(f"git -C {tree} commit --amend --no-edit", repo)


def no_head(repo: Path, tmp_path: Path) -> Line:
    git(tmp_path, "init", "-q", "fresh")
    return Line("git commit --amend --no-edit", tmp_path / "fresh")


def no_head_by_path(repo: Path, tmp_path: Path) -> Line:
    git(tmp_path, "init", "-q", "fresh")
    return Line(f"git -C {tmp_path / 'fresh'} commit --amend --no-edit", repo)


def no_cwd(repo: Path, tmp_path: Path) -> Line:
    local_commit(repo)
    return Line("git commit --amend --no-edit", None)


def pushed_head_abbreviated(option: str) -> Case:
    def case(repo: Path, tmp_path: Path) -> Line:
        return Line(f"git commit {option} --no-edit", repo)

    return case


H27_ASKED = {
    "5-pushed-head": pushed_head,
    "5-pushed-head-absolute-C": pushed_head_by_path,
    "6-tagged": tagged,
    "7-pushed-as-another-branch": pushed_elsewhere,
    "7-pushed-from-a-worktree": pushed_from_worktree,
    "8-no-head": no_head,
    "8-no-head-absolute-C": no_head_by_path,
    "9-GIT_DIR-assignment": not_bare("GIT_DIR={repo}/.git git commit --amend --no-edit"),
    "9-config-override": not_bare("git -c a.b=c commit --amend --no-edit"),
    # `repo` relative to cwd names the unpushed repository, so only the rule refuses it.
    "9-relative-C": not_bare("git -C repo commit --amend --no-edit", cwd="parent"),
    "9-double-quotes": not_bare('git commit --amend -m "x"'),
    "9-no-verify": not_bare("git commit --amend --no-edit --no-verify"),
    "10-no-cwd": no_cwd,
    "11-amen-pushed": pushed_head_abbreviated("--amen"),
    "11-am-pushed": pushed_head_abbreviated("--am"),
}


def assert_amend_form(command: str, reason: str) -> None:
    """h29 §3.3: the form is named exactly when the history rewrite is the only reason."""
    if "--no-verify" in command:
        assert AMEND_FORM not in reason, reason
    else:
        assert AMEND_FORM in reason, reason


def judge(repo: Path, line: Line, mode: str) -> tuple[str, str] | None:
    set_mode(repo, mode)
    event = bash_event(repo, line.command)
    if line.cwd is None:
        del event["cwd"]
    else:
        event["cwd"] = str(line.cwd)
    return pretool_decision(run_script(repo, GUARD_PUSH, event))


@pytest.mark.parametrize("mode", ["off", "on", "expired"])
@pytest.mark.parametrize("case", H27_PASSED.values(), ids=H27_PASSED.keys())
def test_h27_an_amend_of_an_unpushed_head_passes_in_every_mode(
    pushed: Path, tmp_path: Path, case: Case, mode: str
) -> None:
    line = case(pushed, tmp_path)
    assert judge(pushed, line, mode) is None, line.command
    assert queue_lines(pushed) == []


@pytest.mark.parametrize("mode", ["off", "expired"])
@pytest.mark.parametrize("case", H27_ASKED.values(), ids=H27_ASKED.keys())
def test_h27_an_amend_not_shown_unpushed_asks_while_attended(
    pushed: Path, tmp_path: Path, case: Case, mode: str
) -> None:
    line = case(pushed, tmp_path)
    found = judge(pushed, line, mode)
    assert found is not None, f"{line.command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert reason.startswith("This reaches beyond the working tree: ")
    assert HISTORY in reason
    if "--no-verify" in line.command:
        assert NO_VERIFY in reason
    assert_amend_form(line.command, reason)
    assert queue_lines(pushed) == []


@pytest.mark.parametrize("case", H27_ASKED.values(), ids=H27_ASKED.keys())
def test_h27_an_amend_not_shown_unpushed_is_refused_and_queued_while_away(
    pushed: Path, tmp_path: Path, case: Case
) -> None:
    line = case(pushed, tmp_path)
    found = judge(pushed, line, "on")
    assert found is not None, f"{line.command!r} passed silently"
    kind, reason = found
    assert kind == "deny"
    assert reason.startswith("Refused: unattended mode is on until ")
    assert HISTORY in reason
    assert_amend_form(line.command, reason)
    [queued] = queue_lines(pushed)
    assert (queued["hook"], queued["act"]) == ("guard_push", line.command)
    assert HISTORY in queued["why"]


# ---------------------------------------------------------------- h29
#
# docs/increments/h29-review-rules.md §3.3 and §4 rows 6 to 17: besides h27's
# bare form, an amend passes on a line of simple commands joined only by `&&`:
# at most one `cd <absolute path>` first, any `git add`, and the amend last.
# The repository judged is the amend's `-C` path, else the `cd` path, else cwd.
# Every other shape asks by day and is refused and queued while away, and its
# reason names the passing form (AMEND_FORM). Same fixture as h27: one pushed
# commit and, where a case calls `local_commit`, an unpushed HEAD.


def pushed_head_after_cd(repo: Path, tmp_path: Path) -> Line:
    return Line(f"cd {repo} && git commit --amend --no-edit", repo)


H29_PASSED = {
    # cwd is the parent, which is no repository: the `cd` path is the one judged.
    "6-cd-then-amend": not_bare("cd {repo} && git commit --amend --no-edit", "parent"),
    "7-cd-adds-then-amend": not_bare(
        "cd {repo} && git add f && git add g && git commit -q --amend --no-edit -a", "parent"
    ),
    "8-add-then-amend-by-C": not_bare(
        "git add f && git -C {repo} commit --amend --no-edit", "repo"
    ),
    "9-and-inside-single-quotes": not_bare("cd {repo} && git commit --amend -m 'a && b'", "parent"),
}

H29_ASKED = {
    "10-pushed-head-after-cd": pushed_head_after_cd,
    # `repo` relative to the parent names the unpushed repository: only the rule refuses it.
    "11-relative-cd": not_bare("cd repo && git commit --amend --no-edit", "parent"),
    "12-two-cd": not_bare("cd {repo} && cd {repo} && git commit --amend", "repo"),
    "13-semicolon": not_bare("cd {repo}; git commit --amend", "repo"),
    "13-or-list": not_bare("cd {repo} || git commit --amend", "repo"),
    "13-pipe-before": not_bare("git add f | cat && git commit --amend", "repo"),
    "14-other-command-before": not_bare(
        "cd {repo} && git checkout -q HEAD~0 && git commit --amend", "repo"
    ),
    "15-not-last": not_bare(
        "cd {repo} && git commit --amend --no-edit && git log --oneline -1", "repo"
    ),
    "16-rewrite-after": not_bare("cd {repo} && git commit --amend && git rebase master", "repo"),
}


@pytest.mark.parametrize("mode", ["off", "on", "expired"])
@pytest.mark.parametrize("case", H29_PASSED.values(), ids=H29_PASSED.keys())
def test_h29_an_amend_last_after_cd_and_git_add_passes_in_every_mode(
    pushed: Path, tmp_path: Path, case: Case, mode: str
) -> None:
    line = case(pushed, tmp_path)
    assert judge(pushed, line, mode) is None, line.command
    assert queue_lines(pushed) == []


@pytest.mark.parametrize("mode", ["off", "expired"])
@pytest.mark.parametrize("case", H29_ASKED.values(), ids=H29_ASKED.keys())
def test_h29_any_other_amend_line_asks_naming_the_passing_form(
    pushed: Path, tmp_path: Path, case: Case, mode: str
) -> None:
    line = case(pushed, tmp_path)
    found = judge(pushed, line, mode)
    assert found is not None, f"{line.command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert HISTORY in reason and AMEND_FORM in reason, reason
    assert queue_lines(pushed) == []


@pytest.mark.parametrize("case", H29_ASKED.values(), ids=H29_ASKED.keys())
def test_h29_any_other_amend_line_is_refused_and_queued_while_away(
    pushed: Path, tmp_path: Path, case: Case
) -> None:
    line = case(pushed, tmp_path)
    found = judge(pushed, line, "on")
    assert found is not None, f"{line.command!r} passed silently"
    kind, reason = found
    assert kind == "deny"
    assert reason.startswith("Refused: unattended mode is on until ")
    assert HISTORY in reason and AMEND_FORM in reason, reason
    [queued] = queue_lines(pushed)
    assert (queued["hook"], queued["act"]) == ("guard_push", line.command)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("unpushed", [False, True], ids=["pushed-head", "unpushed-head"])
def test_h29_17_a_rewrite_that_is_no_amend_does_not_name_the_amend_form(
    pushed: Path, mode: str, unpushed: bool
) -> None:
    if unpushed:
        local_commit(pushed)
    found = judge(pushed, Line("git rebase master", pushed), mode)
    assert found is not None, "git rebase master passed silently"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert HISTORY in reason
    assert AMEND_FORM not in reason, reason
