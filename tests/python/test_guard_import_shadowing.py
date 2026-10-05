"""Guard import shadowing: a `tools/` file named after a stdlib module (h16 G4).

`docs/increments/h16-harness-fixes.md` §2 G4. The hooks put `tools/` on
`sys.path` to import their helpers; put first, a `tools/dataclasses.py`
replaces the standard library's for every later import, and a guard that
crashes at import exits 1 with no decision, which Claude Code treats as a
non-blocking error. The fix appends `tools/` instead, so the probe below must
leave each hook's decision unchanged; and a stdlib-named `tools/` file is
governed (tested in `test_guard_governance.py`) and absent from the repository.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

from brief_fixtures import make_brief_repo, run_hook
from harness_fixtures import (
    GUARD_GOVERNANCE,
    GUARD_PUSH,
    GUARD_UNATTENDED,
    REAL,
    bash_event,
    file_event,
    pretool_decision,
    queue_lines,
    run_script,
)

#: The planted module: importing it at all is the failure.
PROBE = 'raise RuntimeError("planted tools/dataclasses.py was imported")\n'

#: The seven `sys.path.insert(0, <tools>)` lines §2 G4 part 1 turns into appends.
APPENDERS = (
    ".claude/hooks/guard_push.py",
    ".claude/hooks/guard_governance.py",
    ".claude/hooks/guard_spawn.py",
    ".claude/hooks/guard_unattended.py",
    "tools/session_state.py",
    "tools/brief.py",
    "tools/away.py",
)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """The fixture repository (it holds the hooks and their tools), with the probe planted."""
    root = make_brief_repo(tmp_path.resolve() / "repo")
    (root / "tools" / "dataclasses.py").write_text(PROBE)
    return root


def test_guard_push_still_asks_before_a_push(repo: Path) -> None:
    result = run_script(repo, GUARD_PUSH, bash_event(repo, "git push"))
    assert "planted" not in result.stderr
    found = pretool_decision(result)
    assert found is not None and found[0] == "ask"
    assert "git push writes to the remote" in found[1]


def test_guard_governance_still_asks_before_a_rule_file_write(repo: Path) -> None:
    event = file_event(repo, "Write", str(repo / "CLAUDE.md"))
    result = run_script(repo, GUARD_GOVERNANCE, event)
    assert "planted" not in result.stderr
    found = pretool_decision(result)
    assert found is not None and found[0] == "ask"
    assert found[1].startswith(f"This changes a file that states rules: {repo / 'CLAUDE.md'}.")


def test_guard_unattended_stays_silent_while_attended(repo: Path) -> None:
    """Attended, a question passes; a crash in its import would refuse it instead."""
    event = {
        "hook_event_name": "PreToolUse",
        "tool_name": "AskUserQuestion",
        "tool_input": {"questions": [{"question": "Which?", "header": "Q", "options": []}]},
        "cwd": str(repo),
    }
    result = run_script(repo, GUARD_UNATTENDED, event)
    assert pretool_decision(result) is None, result.stdout
    assert queue_lines(repo) == []


def test_guard_spawn_still_runs_the_brief_check(repo: Path) -> None:
    """With brief.py importable, a tester spawn without a block is refused, not waved through."""
    event = {
        "hook_event_name": "PreToolUse",
        "tool_name": "Agent",
        "tool_input": {"description": "a step", "prompt": "Write the tests.",
                       "subagent_type": "tester"},
        "cwd": str(repo),
    }  # fmt: skip
    result = run_hook(repo, event, repo)
    assert "brief check unavailable" not in result.stdout, result.stdout
    found = pretool_decision(result)
    assert found is not None, result.stdout
    kind, reason = found
    assert (kind, reason.split(";")[0]) == ("deny", "brief: no block for @tester")


@pytest.mark.parametrize("relative", APPENDERS)
def test_tools_is_appended_to_sys_path_not_put_first(relative: str) -> None:
    source = (REAL / relative).read_text()
    assert not re.search(r"sys\.path\.insert\(\s*0\s*,", source), (
        f"{relative} puts a directory first on sys.path; append tools/ instead (h16 G4)"
    )
    assert "sys.path.append(" in source


@pytest.mark.parametrize("directory", ["tools", ".claude/hooks"])
def test_no_file_in_the_import_directories_is_named_after_a_stdlib_module(
    directory: str,
) -> None:
    shadowing = sorted(
        entry.name
        for entry in (REAL / directory).iterdir()
        if entry.name.split(".")[0] in sys.stdlib_module_names
    )
    assert shadowing == [], f"{directory}/ shadows the standard library: {shadowing}"
