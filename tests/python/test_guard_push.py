"""`.claude/hooks/guard_push.py` in and out of unattended mode (h3 U1, T4, T5, T8-T10; h10 §4a).

The spec is `docs/increments/h3-unattended-u1.md` §3.5 (the new patterns),
§3.9 (failure direction) and §4. Each command is run through a copy of the hook
in a temporary repository whose flag the test sets, so the hook's state dir is
`<tmp>/.git/harness`. By day an asked act asks; at night it is denied and
queued; a pass is silent in both.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from harness_fixtures import (
    GUARD_GOVERNANCE,
    GUARD_PUSH,
    SUBAGENT,
    bash_event,
    file_event,
    make_repo,
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
