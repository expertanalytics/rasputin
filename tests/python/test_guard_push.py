"""`.claude/hooks/guard_push.py` in and out of unattended mode (h3 U1, T4, T5, T8-T10; h10 §4a).

The spec is `docs/increments/h3-unattended-u1.md` §3.5 (the new patterns),
§3.9 (failure direction) and §4. Each command is run through a copy of the hook
in a temporary repository whose flag the test sets, so the hook's state dir is
`<tmp>/.git/harness`. By day an asked act asks; at night it is denied and
queued; a pass is silent in both.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from harness_fixtures import (
    GUARD_GOVERNANCE,
    GUARD_PUSH,
    SUBAGENT,
    add_worktree,
    bash_event,
    clean_env,
    file_event,
    git,
    make_repo,
    point_scratchpad,
    pretool_decision,
    queue_lines,
    run_script,
    set_mode,
)

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
    "git commit --amend": None,
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


# ---------------------------------------------------------------- h16 G3b
#
# §2 G3b: local ref, remote and config writes made with `git -C <dir>` pass
# when both of `<dir>`'s git dirs lie under a session scratchpad, and nothing
# in the line redirects git elsewhere. A push is never exempt (G3c dropped).


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
    root.mkdir(parents=True)
    git(root, "init", "-q", "-b", "master")
    git(root, "-c", "user.name=T", "-c", "user.email=t@example.invalid",
        "commit", "-q", "--allow-empty", "-m", "root")  # fmt: skip
    return root


#: Local writes in the scratch repository, `{inside}` its absolute path.
G3B_PASSED = (
    "git -C {inside} config user.name x",
    "git -C {inside} config core.hooksPath x",
    "git -C {inside} remote add o /x",
    "git -C {inside} remote set-url o /y",
    "git -C {inside} update-ref refs/heads/y HEAD",
    "git -C {inside} symbolic-ref HEAD refs/heads/y",
    "git -C {inside} fetch . HEAD:refs/remotes/origin/master",
    "git -C {inside} replace HEAD HEAD~1",
    "git config --file {inside}/.git/config user.name x",
    "cd {pad} && git -C {inside} config user.name x",
)

#: The same writes where something sends git outside the scratchpad, with the reason.
G3B_ASKED: dict[str, str] = {
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
@pytest.mark.parametrize("template", G3B_PASSED)
def test_a_local_write_in_a_scratch_repository_is_silent(
    repo: Path, scratch: dict[str, Path], mode: str, template: str
) -> None:
    set_mode(repo, mode)
    command = template.format(**scratch)
    result = run_script(repo, GUARD_PUSH, bash_event(repo, command))
    assert pretool_decision(result) is None, f"{command!r} was not passed"
    assert queue_lines(repo) == []


@pytest.mark.parametrize("template", G3B_ASKED)
def test_a_scratch_write_that_reaches_outside_still_asks(
    repo: Path, scratch: dict[str, Path], template: str
) -> None:
    command = template.format(**scratch)
    found = pretool_decision(run_script(repo, GUARD_PUSH, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert G3B_ASKED[template] in reason
