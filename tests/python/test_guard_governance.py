"""`.claude/hooks/guard_governance.py` in and out of unattended mode (h3 U1, T6, T7).

The spec is `docs/increments/h3-unattended-u1.md` §3.6 and §4. Two changes are
pinned: the self-protecting set joins the governed files (asked by day, denied
and queued at night), and the harness state and `away.py` are denied to every
agent in both modes, with nothing queued.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from harness_fixtures import (
    GUARD_GOVERNANCE,
    bash_event,
    file_event,
    make_repo,
    pretool_decision,
    queue_lines,
    run_script,
    set_mode,
)

ALWAYS_DENIED = (
    "Only Ola enters or leaves unattended mode, and only hooks and away.py write the "
    "harness state. Nothing is queued: this act is not an agent's to wait for."
)

#: T6: CLAUDE.md, and each new member of the self-protecting set, relative to the repo.
GOVERNED_NOW = (
    "CLAUDE.md",
    "tools/away.py",
    "tools/harness_mode.py",
    "tools/session_state.py",
    ".claude/profile.toml",
    ".claude/skills/x/SKILL.md",
    ".git/hooks/pre-commit",
    ".git/config",
)

#: Files named `config` that are not a git config: GOVERNED_SUFFIXES is matched
#: with endswith, never on the basename.
NOT_GOVERNED = ("src/app/config", "docs/config")

#: T7: Bash commands that write the harness state or run away.py.
DENIED_COMMANDS = (
    "echo {} > .git/harness/unattended.json",
    "rm .git/harness/unattended.json",
    "python3 tools/away.py 8h",
    "tools/away.py --back",
    "script -q /dev/null python3 tools/away.py 8h",
)

#: T7: reads of the same paths, which must stay silent.
SILENT_COMMANDS = (
    "cat .git/harness/queue.jsonl",
    "cat tools/away.py",
    "git diff tools/away.py",
    "pytest tests/python/test_away.py",
)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


# ---------------------------------------------------------------- T6


@pytest.mark.parametrize("tool", ["Edit", "Write"])
@pytest.mark.parametrize("relative", GOVERNED_NOW)
def test_a_governed_file_asks_while_attended(repo: Path, tool: str, relative: str) -> None:
    path = str(repo / relative)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, file_event(repo, tool, path)))
    assert found is not None, f"{relative} is not governed"
    assert found[0] == "ask"
    assert queue_lines(repo) == []


@pytest.mark.parametrize("tool", ["Edit", "Write"])
@pytest.mark.parametrize("relative", GOVERNED_NOW)
def test_a_governed_file_is_denied_and_queued_while_unattended(
    repo: Path, tool: str, relative: str
) -> None:
    set_mode(repo, "on")
    path = str(repo / relative)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, file_event(repo, tool, path)))
    assert found is not None, f"{relative} is not governed"
    kind, reason = found
    assert kind == "deny"
    assert reason.startswith("Refused: unattended mode is on until ")
    assert "a rule file changes" in reason
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"], line["why"]) == (
        "guard_governance",
        f"{tool} {path}",
        "a rule file changes",
    )


def test_a_governed_bash_write_is_queued_with_the_command_as_act(repo: Path) -> None:
    set_mode(repo, "on")
    command = "echo x >> CLAUDE.md"
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found is not None and found[0] == "deny"
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_governance", command)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("relative", NOT_GOVERNED)
def test_a_file_merely_named_config_is_not_governed(repo: Path, mode: str, relative: str) -> None:
    set_mode(repo, mode)
    event = file_event(repo, "Edit", str(repo / relative))
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- T7


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("tool", ["Write", "Edit"])
def test_a_file_tool_write_to_the_harness_state_is_always_denied(
    repo: Path, mode: str, tool: str
) -> None:
    set_mode(repo, mode)
    event = file_event(repo, tool, str(repo / ".git" / "harness" / "unattended.json"))
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", DENIED_COMMANDS)
def test_a_bash_write_of_the_state_or_a_run_of_away_is_always_denied(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
def test_a_subagent_running_away_is_denied_and_not_queued(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    event = bash_event(repo, "python3 tools/away.py 8h", agent_id="a1", agent_type="developer")
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", SILENT_COMMANDS)
def test_a_read_of_the_state_or_of_away_is_silent(repo: Path, mode: str, command: str) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []
