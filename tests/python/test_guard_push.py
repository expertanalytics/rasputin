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
    root.mkdir(parents=True)
    git(root, "init", "-q", "-b", "master")
    git(root, "-c", "user.name=T", "-c", "user.email=t@example.invalid",
        "commit", "-q", "--allow-empty", "-m", "root")  # fmt: skip
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
    "find . -maxdepth 0 -exec git push origin HEAD \;": PUSH,
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
    "find . -maxdepth 0 -exec sh -c 'git push' \;": PUSH,
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
    # Pinned passes, the residual §6 keeps: an alias through a bare tail, and
    # an interpreter program that runs git.
    "caffeinate -i git p origin",
    "caffeinate -i gh pm 12",
    "python3 -c \"import subprocess; subprocess.run(['git', 'push'])\"",
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
