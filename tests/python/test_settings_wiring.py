"""`.claude/settings.json` wires the harness hooks, and each one runs as wired.

The spec is `docs/increments/h3-unattended-u1.md` §3.11 (the three new
entries for `guard_unattended.py`, the existing wiring kept) and §5's fallback
paragraph (the settings PR also restores the `guard_unattended.py` clause to
`.claude/REQUIRED-READING.md`, *The harness*).

Every test reads the real files of the checkout. The by-path runs do not: each
hook is run from a copy in a temporary repository (`harness_fixtures.make_repo`),
with `$CLAUDE_PROJECT_DIR` pointing at that copy, so a hook's state dir is
`<tmp>/.git/harness` and the real repository's harness state is never read or
written. The copy keeps the real file's mode (`shutil.copy2`), and the real
file's mode is asserted separately, both on disk and in git's index.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import REAL, bash_event, clean_env, git, is_work_tree_top, make_repo

SETTINGS = REAL / ".claude" / "settings.json"
REQUIRED_READING = REAL / ".claude" / "REQUIRED-READING.md"

PROJECT = "$CLAUDE_PROJECT_DIR"
UNATTENDED = f"{PROJECT}/.claude/hooks/guard_unattended.py"
GOVERNANCE = f"{PROJECT}/.claude/hooks/guard_governance.py"
PUSH = f"{PROJECT}/.claude/hooks/guard_push.py"
GATES = f"{PROJECT}/.claude/hooks/gates_after_commit.py"
SPAWN = f"{PROJECT}/.claude/hooks/guard_spawn.py"
SESSION_STATE = 'python3 "$CLAUDE_PROJECT_DIR/tools/session_state.py"'

CONFIG_SOURCES = "user_settings|project_settings|local_settings|skills"

#: (event, matcher, command); a matcher of None means the entry has none.
NEW_WIRING = [
    pytest.param("PreToolUse", "AskUserQuestion", UNATTENDED, id="pretooluse-askuserquestion"),
    pytest.param("PermissionRequest", None, UNATTENDED, id="permissionrequest-no-matcher"),
    pytest.param("ConfigChange", CONFIG_SOURCES, UNATTENDED, id="configchange-sources"),
]

#: h9 §6 (docs/increments/h9-spawn-briefs.md): one entry, after AskUserQuestion's.
H9_WIRING = [
    pytest.param("PreToolUse", "Agent|SendMessage", SPAWN, id="guard-spawn"),
]

EXISTING_WIRING = [
    pytest.param("SessionStart", None, SESSION_STATE, id="sessionstart-session-state"),
    pytest.param("PreToolUse", "Edit|Write|NotebookEdit|Bash", GOVERNANCE, id="guard-governance"),
    pytest.param("PreToolUse", "Bash", PUSH, id="guard-push"),
    pytest.param("PostToolUse", "Bash", GATES, id="gates-after-commit"),
]

#: Exit codes a shell returns when it cannot execute the command at all.
NOT_EXECUTED = {126: "found but not executable", 127: "not found"}


def settings() -> dict[str, Any]:
    loaded = json.loads(SETTINGS.read_text())
    assert isinstance(loaded, dict)
    return loaded


def wired_hooks() -> list[tuple[str, str | None, str]]:
    """Every command hook in settings.json as (event, matcher or None, command)."""
    found: list[tuple[str, str | None, str]] = []
    for event, entries in settings().get("hooks", {}).items():
        for entry in entries:
            matcher = entry.get("matcher") or None  # absent and "" both match every event
            for hook in entry.get("hooks", []):
                if hook.get("type") == "command":
                    found.append((event, matcher, hook["command"]))
    return found


@pytest.mark.parametrize(("event", "matcher", "command"), NEW_WIRING)
def test_guard_unattended_is_wired(event: str, matcher: str | None, command: str) -> None:
    """§3.11 items 1-3: each entry, with its matcher exactly as specified."""
    wired = wired_hooks()
    assert (event, matcher, command) in wired, (
        f"{event} (matcher {matcher!r}) does not run {command}; "
        f"{event} runs {[(m, c) for e, m, c in wired if e == event]}"
    )


@pytest.mark.parametrize(("event", "matcher", "command"), H9_WIRING)
def test_guard_spawn_is_wired(event: str, matcher: str | None, command: str) -> None:
    """h9 test 24: the entry, with its matcher exactly as §6 gives it."""
    wired = wired_hooks()
    assert (event, matcher, command) in wired, (
        f"{event} (matcher {matcher!r}) does not run {command}; "
        f"{event} runs {[(m, c) for e, m, c in wired if e == event]}"
    )


def test_guard_spawn_is_the_only_new_pretooluse_entry() -> None:
    """h9 test 24: the existing PreToolUse entries are unchanged, and kept in order."""
    entries = settings()["hooks"]["PreToolUse"]
    matchers = [entry.get("matcher") for entry in entries]
    assert matchers.count("Agent|SendMessage") == 1
    spawn = matchers.index("Agent|SendMessage")
    assert matchers.index("AskUserQuestion") == spawn - 1, matchers


@pytest.mark.parametrize(("event", "matcher", "command"), EXISTING_WIRING)
def test_existing_wiring_is_kept(event: str, matcher: str | None, command: str) -> None:
    """§3.11: `guard_push.py`, `guard_governance.py` and the rest keep their wiring."""
    assert (event, matcher, command) in wired_hooks()


# -- by-path invocation -------------------------------------------------------


def bare_path_hooks() -> list[Any]:
    """The command hooks whose command is a path alone, with no interpreter prefix.

    Claude Code hands such a command to a shell, which executes the file itself;
    a file without its executable bit then fails with 126 before any Python runs.
    The spec's entries are included whether or not settings.json has them yet,
    so `guard_unattended.py` is run by path before it is wired.
    """
    specified = [tuple(p.values) for p in (*NEW_WIRING, *H9_WIRING, *EXISTING_WIRING)]
    params = []
    for event, matcher, command in dict.fromkeys([*wired_hooks(), *specified]):
        assert isinstance(event, str) and isinstance(command, str)
        assert matcher is None or isinstance(matcher, str)
        if command.startswith(f"{PROJECT}/") and " " not in command:
            name = Path(command).name
            params.append(pytest.param(event, matcher, command, id=f"{event}-{name}"))
    return params


def harmless_event(repo: Path, event: str, matcher: str | None) -> dict[str, Any]:
    """An event of the hook's kind that no hook should act on, by day."""
    tools = set((matcher or "").split("|"))
    if event in {"PreToolUse", "PostToolUse"} and "Bash" in tools:
        return {**bash_event(repo, "ls"), "hook_event_name": event}
    if event == "PreToolUse" and tools == {"AskUserQuestion"}:
        question = {"question": "Which one?", "header": "Q", "options": [{"label": "A"}]}
        return {
            "hook_event_name": event,
            "tool_name": "AskUserQuestion",
            "tool_input": {"questions": [question]},
            "cwd": str(repo),
        }
    if event == "PreToolUse" and tools == {"Agent", "SendMessage"}:
        return {
            "hook_event_name": event,
            "tool_name": "SendMessage",
            "tool_input": {"to": "a1b2c3d4", "message": "Carry on."},
            "cwd": str(repo),
        }
    if event == "PermissionRequest":
        return {
            "hook_event_name": event,
            "tool_name": "Bash",
            "tool_input": {"command": "ls", "description": "List files"},
            "cwd": str(repo),
        }
    if event == "ConfigChange":
        return {
            "hook_event_name": event,
            "source": "project_settings",
            "file_path": str(repo / ".claude" / "settings.json"),
            "cwd": str(repo),
        }
    pytest.fail(f"no harmless event is defined for {event} (matcher {matcher!r}); add one")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


def test_every_hook_script_is_covered_by_path() -> None:
    """The five hooks under `.claude/hooks/` are all wired as bare paths.

    Without this, a hook rewired with an interpreter prefix would drop out of
    the by-path test below unnoticed.
    """
    scripts = {Path(command).name for _, _, command in wired_hooks() if " " not in command}
    assert scripts == {
        "guard_governance.py",
        "guard_push.py",
        "gates_after_commit.py",
        "guard_unattended.py",
        "guard_spawn.py",
    }


@pytest.mark.parametrize(("event", "matcher", "command"), bare_path_hooks())
def test_hook_is_executable_in_the_checkout(event: str, matcher: str | None, command: str) -> None:
    """The real file has its executable bit, on disk and in git's index (100755)."""
    relative = command.removeprefix(f"{PROJECT}/")
    path = REAL / relative
    assert path.is_file(), f"{relative} is wired on {event} but missing"
    assert os.access(path, os.X_OK), f"{relative} is not executable (mode {path.stat().st_mode:o})"
    if not is_work_tree_top(REAL):  # h16 P9: a `git archive` copy has no index to read
        pytest.skip(f"{REAL} is not the top of a git work tree, so git's index mode is not checked")
    staged = git(REAL, "ls-files", "-s", "--", relative).split()
    assert staged, f"{relative} is not tracked"
    assert staged[0] == "100755", f"{relative} is tracked at mode {staged[0]}, not 100755"


@pytest.mark.parametrize(("event", "matcher", "command"), bare_path_hooks())
def test_hook_runs_by_path(repo: Path, event: str, matcher: str | None, command: str) -> None:
    """The command, run by a shell as Claude Code runs it, exits 0 on a harmless event.

    `$CLAUDE_PROJECT_DIR` is the temporary copy, whose flag is absent (mode
    off), so the run cannot touch the real repository's harness state.
    """
    relative = command.removeprefix(f"{PROJECT}/")
    copy = repo / relative
    if not copy.exists():  # make_repo copies the U1 scripts; the rest come here
        copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REAL / relative, copy)
    result = subprocess.run(
        ["/bin/sh", "-c", command],
        input=json.dumps(harmless_event(repo, event, matcher)),
        capture_output=True,
        text=True,
        cwd=repo,
        env={**clean_env(), "CLAUDE_PROJECT_DIR": str(repo)},
        timeout=60,
        check=False,
    )
    code = result.returncode
    assert code not in NOT_EXECUTED, f"{relative}: exit {code}, {NOT_EXECUTED.get(code)}"
    assert code == 0, f"{relative}: exit {code}; stderr: {result.stderr}"
    assert not (repo / ".git" / "harness" / "queue.jsonl").exists(), "an attended run queued"


def test_guard_spawn_by_path_denies_a_blockless_persona_spawn(repo: Path) -> None:
    """h9 test 24: the wired command, run by a shell from the fixture's copy."""
    relative = SPAWN.removeprefix(f"{PROJECT}/")
    assert (repo / relative).exists(), f"{relative} is missing from the copy"
    event = {
        "hook_event_name": "PreToolUse",
        "tool_name": "Agent",
        "tool_input": {
            "description": "red",
            "prompt": "Write the tests.",
            "subagent_type": "tester",
        },
        "cwd": str(repo),
    }
    result = subprocess.run(
        ["/bin/sh", "-c", SPAWN],
        input=json.dumps(event),
        capture_output=True,
        text=True,
        cwd=repo,
        env={**clean_env(), "CLAUDE_PROJECT_DIR": str(repo)},
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    output = json.loads(result.stdout)["hookSpecificOutput"]
    assert output["permissionDecision"] == "deny"
    assert output["permissionDecisionReason"].startswith("brief: no block")


# -- the rule text ------------------------------------------------------------


def harness_hook_list() -> str:
    """The paragraph of *The harness* that lists the hooks active in settings.json."""
    text = REQUIRED_READING.read_text()
    section = re.search(r"^## The harness\n(.*?)(?=^## )", text, re.MULTILINE | re.DOTALL)
    assert section is not None, "REQUIRED-READING.md has no '## The harness' section"
    paragraphs = [p for p in section.group(1).split("\n\n") if p.strip()]
    listing = [p for p in paragraphs if p.lstrip().startswith("Active in `.claude/settings.json`")]
    assert len(listing) == 1, "The harness section has no single 'Active in' hook list"
    return listing[0]


def test_required_reading_lists_guard_unattended() -> None:
    """§5: the clause naming `guard_unattended.py` returns with the wiring."""
    assert "guard_unattended.py" in harness_hook_list()


def test_required_reading_lists_the_existing_hooks() -> None:
    listing = harness_hook_list()
    for name in ("guard_push.py", "guard_governance.py", "gates_after_commit.py"):
        assert name in listing
