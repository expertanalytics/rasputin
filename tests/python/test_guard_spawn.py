"""`.claude/hooks/guard_spawn.py`: spawns and resumes are checked (h9, tests 10 to 22).

Spec: `docs/increments/h9-spawn-briefs.md` §3.5 to §3.8. The hook is run by
path from the fixture's copy (`brief_fixtures.make_brief_repo`), with
`$CLAUDE_PROJECT_DIR` set to the fixture, as `test_settings_wiring.py` runs the
other hooks. Blocks are made by `brief_fixtures.make_block`, which writes §3.2
from the design's text, so this suite needs `tools/brief.py` only where the
hook itself imports it; one test feeds the hook `brief.py`'s own output.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

import pytest

from brief_fixtures import (
    BRIEF,
    END,
    HEADER,
    human,
    make_block,
    make_brief_repo,
    run_brief,
    run_hook,
    write_transcript,
)
from harness_fixtures import add_worktree, flag_on, git, queue_path, state_dir

NO_BLOCK = "brief: no block"


@pytest.fixture
def project(tmp_path: Path) -> Path:
    root = make_brief_repo(tmp_path.resolve() / "repo")
    (root / "docs").mkdir(exist_ok=True)
    return root


@pytest.fixture
def tree(tmp_path: Path, project: Path) -> Path:
    return add_worktree(project, tmp_path.resolve() / "A", "wt-a")


def agent(cwd: Path | str | None, prompt: str, subagent_type: str | None = "tester",
          **tool: Any) -> dict[str, Any]:  # fmt: skip
    """An `Agent` PreToolUse event with the fields §3.5 reads."""
    tool_input: dict[str, Any] = {"description": "a step", "prompt": prompt, **tool}
    if subagent_type is not None:
        tool_input["subagent_type"] = subagent_type
    event: dict[str, Any] = {
        "hook_event_name": "PreToolUse",
        "tool_name": "Agent",
        "tool_input": tool_input,
    }
    if cwd is not None:
        event["cwd"] = str(cwd)
    return event


def resume(cwd: Path | str, message: str) -> dict[str, Any]:
    """A `SendMessage` PreToolUse event."""
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": "SendMessage",
        "tool_input": {"to": "a1b2c3d4", "message": message},
        "cwd": str(cwd),
    }


def verdict(result: subprocess.CompletedProcess[str]) -> tuple[str | None, str, str]:
    """(permissionDecision or None, reason, additionalContext); ('', '', '') when silent."""
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    if not result.stdout.strip():
        return None, "", ""
    output = json.loads(result.stdout)["hookSpecificOutput"]
    assert output["hookEventName"] == "PreToolUse"
    return (
        output.get("permissionDecision"),
        output.get("permissionDecisionReason", ""),
        output.get("additionalContext", ""),
    )


def silent(result: subprocess.CompletedProcess[str]) -> None:
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    assert result.stdout.strip() == "", f"expected no output, got {result.stdout!r}"


def denied(result: subprocess.CompletedProcess[str], starts: str) -> str:
    decision, reason, _ = verdict(result)
    assert decision == "deny", f"not denied: {result.stdout!r}"
    assert reason.startswith(starts), f"reason {reason!r} does not start with {starts!r}"
    return reason


def task(block: str) -> str:
    return f"{block}\n\nThe task: write the red tests for x."


# ---------------------------------------------------------------- 10


def test_10_a_fresh_block_on_its_persona_passes(project: Path) -> None:
    silent(run_hook(project, agent(project, task(make_block("tester", project))), project))


def test_10_brief_pys_own_output_passes(project: Path, tmp_path: Path) -> None:
    """End to end: the block brief.py prints is the block the hook accepts."""
    home = tmp_path.resolve() / "home"
    write_transcript(home, project, [human("Go.", "2026-10-04T06:00:00.000Z")])
    (project / "docs" / "increments").mkdir(parents=True)
    (project / "docs" / "increments" / "x.md").write_text("# x\n\nStatus: design.\n")
    printed = run_brief(project, home, "tester", "--worktree", str(project), "--beside", "none",
                        "--increment", "docs/increments/x.md")  # fmt: skip
    assert printed.returncode == 0, printed.stderr
    silent(run_hook(project, agent(project, task(printed.stdout)), project))


# ---------------------------------------------------------------- 11, rule 5


def test_11_a_persona_spawn_without_a_block_is_denied(project: Path) -> None:
    reason = denied(run_hook(project, agent(project, "Write the tests."), project), NO_BLOCK)
    assert "python3 tools/brief.py" in reason


@pytest.mark.parametrize("persona", ["architect", "developer", "orchestrator", "perf", "reviewer"])
def test_11_every_persona_needs_a_block(project: Path, persona: str) -> None:
    denied(run_hook(project, agent(project, "Do it.", persona), project), NO_BLOCK)


@pytest.mark.parametrize("kind", ["general-purpose", "Explore", None])
def test_11_a_helper_that_is_not_a_persona_needs_no_block(project: Path, kind: str | None) -> None:
    silent(run_hook(project, agent(project, "Find the file.", kind), project))


# ---------------------------------------------------------------- 12, rule 6


def test_12_one_changed_body_character_is_denied_as_edited(project: Path) -> None:
    block = make_block("tester", project, body="Read the files.\nThen write the tests.")
    edited = block.replace("Then write", "Then wrote")
    reason = denied(run_hook(project, agent(project, task(edited)), project), "brief: edited")
    assert "paste the output of tools/brief.py unchanged" in reason


def test_12_a_deleted_body_line_is_denied_as_edited(project: Path) -> None:
    block = make_block("tester", project, body="Read the files.\nRun the mutation round.\nEnd.")
    edited = block.replace("Run the mutation round.\n", "")
    denied(run_hook(project, agent(project, task(edited)), project), "brief: edited")


def test_12_trailing_spaces_and_crlf_pass(project: Path) -> None:
    block = make_block("tester", project, body="Read the files.\nThen write the tests.")
    pasted = "\r\n".join(line + "   " for line in block.split("\n"))
    silent(run_hook(project, agent(project, task(pasted)), project))


# ---------------------------------------------------------------- 13, rule 4


def malformed_variants(project: Path) -> dict[str, str]:
    block = make_block("tester", project)
    header, *body, end = block.split("\n")
    match = HEADER.match(header)
    assert match is not None and END.match(end) is not None
    head = match.group(3)
    other = make_block("tester", project, body="Another block.")
    return {
        "two blocks": f"{block}\n\n{other}",
        "39-digit head": block.replace(f"head={head}", f"head={head[:39]}", 1),
        "no end line": "\n".join([header, *body]),
        "end hash differs": "\n".join([header, *body, "<<<END BRIEF 000000000000>>>"]),
    }


@pytest.mark.parametrize(
    "variant", ["two blocks", "39-digit head", "no end line", "end hash differs"]
)
def test_13_a_malformed_or_doubled_block_is_denied(project: Path, variant: str) -> None:
    text = malformed_variants(project)[variant]
    decision, reason, _ = verdict(run_hook(project, agent(project, task(text)), project))
    assert decision == "deny"
    expected = "brief: two blocks" if variant == "two blocks" else "brief: malformed"
    assert reason.startswith(expected), reason


# ---------------------------------------------------------------- 14, rule 7


def test_14_a_block_for_another_persona_is_denied_naming_both(project: Path) -> None:
    event = agent(project, task(make_block("tester", project)), "developer")
    reason = denied(run_hook(project, event, project), "brief: persona")
    assert "tester" in reason and "developer" in reason


# ---------------------------------------------------------------- 15, rule 8


def test_15_a_commit_after_the_block_makes_it_stale(project: Path, tree: Path) -> None:
    block = make_block("tester", tree)
    git(tree, "commit", "-q", "--allow-empty", "-m", "a later commit")
    reason = denied(run_hook(project, agent(project, task(block)), project), "brief: stale")
    assert "run brief.py again" in reason


def test_15_a_removed_worktree_is_denied(project: Path, tree: Path) -> None:
    block = make_block("tester", tree)
    git(project, "worktree", "remove", "--force", str(tree))
    denied(run_hook(project, agent(project, task(block)), project), "brief: worktree")


def test_15_a_block_for_a_linked_worktree_at_its_head_passes(project: Path, tree: Path) -> None:
    git(tree, "commit", "-q", "--allow-empty", "-m", "ahead of master")
    silent(run_hook(project, agent(project, task(make_block("tester", tree))), project))


# ---------------------------------------------------------------- 16, rule 9


def test_16_a_block_with_isolation_is_denied(project: Path) -> None:
    event = agent(project, task(make_block("tester", project)), isolation="worktree")
    denied(run_hook(project, event, project), "brief: isolation")


# ---------------------------------------------------------------- 17, rule 2


def wrong_places(project: Path, tree: Path, tmp_path: Path) -> dict[str, Path]:
    gone = tmp_path.resolve() / "gone"
    gone.mkdir()
    gone.rmdir()
    return {"worktree": tree, "subdirectory": project / "docs", "deleted": gone}


@pytest.mark.parametrize("place", ["worktree", "subdirectory", "deleted"])
@pytest.mark.parametrize("tool", ["Agent", "SendMessage"])
def test_17_a_spawn_or_resume_outside_the_project_directory_is_denied(
    project: Path, tree: Path, tmp_path: Path, place: str, tool: str
) -> None:
    where = wrong_places(project, tree, tmp_path)[place]
    block = make_block("tester", project)
    event = agent(where, task(block)) if tool == "Agent" else resume(where, "Carry on.")
    reason = denied(run_hook(project, event, project), "cwd:")
    assert str(where) in reason
    assert str(project) in reason
    assert f"run: cd {project}" in reason


def test_17_a_symlink_to_the_project_is_the_project(project: Path, tmp_path: Path) -> None:
    link = tmp_path.resolve() / "link"
    link.symlink_to(project, target_is_directory=True)
    silent(run_hook(project, agent(link, task(make_block("tester", project))), project))


def test_17_the_directory_is_checked_before_the_block(project: Path, tree: Path) -> None:
    """A valid block from the wrong directory is refused for the directory."""
    denied(run_hook(project, agent(tree, task(make_block("tester", project))), project), "cwd:")


# ---------------------------------------------------------------- 18, resumes


def test_18_a_resume_without_a_block_passes(project: Path) -> None:
    silent(run_hook(project, resume(project, "Re-verify after the ruling."), project))


def test_18_a_resume_with_a_fresh_block_passes(project: Path) -> None:
    silent(run_hook(project, resume(project, task(make_block("tester", project))), project))


def test_18_a_resume_with_an_edited_block_is_denied(project: Path) -> None:
    block = make_block("tester", project, body="Read the files.\nThen write the tests.")
    edited = block.replace("Then write", "Then wrote")
    denied(run_hook(project, resume(project, task(edited)), project), "brief: edited")


def test_18_a_resume_with_a_stale_block_is_denied(project: Path, tree: Path) -> None:
    block = make_block("tester", tree)
    git(tree, "commit", "-q", "--allow-empty", "-m", "later")
    denied(run_hook(project, resume(project, task(block)), project), "brief: stale")


@pytest.mark.parametrize("variant", ["bad header", "end without header"])
def test_18_a_resume_with_a_malformed_marker_is_denied(project: Path, variant: str) -> None:
    text = (
        "<<<BRIEF persona=tester worktree=relative/path>>>\nbody\n<<<END BRIEF 000000000000>>>"
        if variant == "bad header"
        else "Some words.\n<<<END BRIEF 0123456789ab>>>"
    )
    denied(run_hook(project, resume(project, text), project), "brief: malformed")


# ---------------------------------------------------------------- 19, rule 1


def test_19_an_event_from_inside_a_subagent_passes(project: Path, tree: Path) -> None:
    event = {**agent(tree, "Do it."), "agent_id": "a1b2c3", "agent_type": "tester"}
    silent(run_hook(project, event, project))


def test_19_another_tool_passes(project: Path, tree: Path) -> None:
    event = {"hook_event_name": "PreToolUse", "tool_name": "Bash",
             "tool_input": {"command": "ls"}, "cwd": str(tree)}  # fmt: skip
    silent(run_hook(project, event, project))


@pytest.mark.parametrize("stdin", ["not json", "[1, 2]", ""])
def test_19_stdin_that_is_not_a_json_object_passes(project: Path, stdin: str) -> None:
    silent(run_hook(project, stdin, project))


# ---------------------------------------------------------------- 20, unavailable


def break_brief(project: Path, how: str) -> None:
    path = project / BRIEF
    if how == "removed":
        path.unlink(missing_ok=True)
    elif how == "raises on import":
        path.write_text('raise RuntimeError("planted import failure")\n')
    elif how == "check raises":
        assert path.exists(), f"{BRIEF} is missing from the copy"
        with path.open("a") as out:
            out.write(
                '\n\ndef check(*args, **kwargs):\n    raise RuntimeError("planted check failure")\n'
            )
    else:  # pragma: no cover
        raise ValueError(how)


@pytest.mark.parametrize("how", ["removed", "raises on import", "check raises"])
def test_20_without_a_working_brief_a_spawn_passes_with_a_notice(project: Path, how: str) -> None:
    break_brief(project, how)
    text = "Write the tests." if how != "check raises" else task(make_block("tester", project))
    decision, _, context = verdict(run_hook(project, agent(project, text), project))
    assert decision in (None, "allow")
    assert "brief check unavailable" in context
    if how == "check raises":
        assert "planted check failure" in context


@pytest.mark.parametrize("how", ["removed", "raises on import", "check raises"])
def test_20_without_a_working_brief_the_directory_is_still_checked(
    project: Path, tree: Path, how: str
) -> None:
    break_brief(project, how)
    denied(run_hook(project, agent(tree, "Write the tests."), project), "cwd:")


# ---------------------------------------------------------------- 21, unattended


def test_21_unattended_mode_changes_no_verdict_and_queues_nothing(project: Path) -> None:
    flag_on(project)
    denied(run_hook(project, agent(project, "Write the tests."), project), NO_BLOCK)
    silent(run_hook(project, agent(project, "Find it.", "general-purpose"), project))
    silent(run_hook(project, agent(project, task(make_block("tester", project))), project))
    assert not queue_path(project).exists()
    assert sorted(p.name for p in state_dir(project).iterdir()) == ["unattended.json"]


# ---------------------------------------------------------------- 22, no comparison


@pytest.mark.parametrize("missing", ["cwd", "CLAUDE_PROJECT_DIR"])
def test_22_without_both_directories_the_comparison_is_skipped_and_said(
    project: Path, missing: str
) -> None:
    cwd = None if missing == "cwd" else project
    env_project = None if missing == "CLAUDE_PROJECT_DIR" else project
    decision, _, context = verdict(
        run_hook(project, agent(cwd, task(make_block("tester", project))), env_project)
    )
    assert decision in (None, "allow")
    assert context, "the skipped comparison is not said"
    assert "cwd" in context or "directory" in context


@pytest.mark.parametrize("missing", ["cwd", "CLAUDE_PROJECT_DIR"])
def test_22_without_both_directories_the_block_checks_still_run(
    project: Path, missing: str
) -> None:
    cwd = None if missing == "cwd" else project
    env_project = None if missing == "CLAUDE_PROJECT_DIR" else project
    denied(run_hook(project, agent(cwd, "Write the tests."), env_project), NO_BLOCK)


# ---------------------------------------------------------------- the fixture itself


def test_the_hook_is_copied_executable(project: Path) -> None:
    """The copy keeps the real file's mode, so a by-path run is a real one."""
    hook = project / ".claude" / "hooks" / "guard_spawn.py"
    assert hook.exists(), ".claude/hooks/guard_spawn.py is missing from the copy"
    assert os.access(hook, os.X_OK)
