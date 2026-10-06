"""`.claude/hooks/guard_governance.py` in and out of unattended mode (h3 U1, T6, T7).

The spec is `docs/increments/h3-unattended-u1.md` §3.6 and §4. Two changes are
pinned: the self-protecting set joins the governed files (asked by day, denied
and queued at night), and the harness state and `away.py` are denied to every
agent in both modes, with nothing queued.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

from harness_fixtures import (
    GUARD_GOVERNANCE,
    REAL,
    bash_event,
    file_event,
    make_repo,
    point_scratchpad,
    pretool_decision,
    queue_lines,
    run_script,
    set_mode,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

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
    # h8 §5: the size table runs inside the SessionStart hook.
    "tools/rule_sizes.py",
    # h9 §5: guard_spawn.py imports brief.py, and the brief template is rule text.
    "tools/brief.py",
    ".claude/briefs/common.md",
    ".claude/briefs/x.md",
    ".claude/profile.toml",
    ".claude/skills/x/SKILL.md",
    ".git/hooks/pre-commit",
    ".git/config",
    # h16 §3: the counter computes CLAUDE.md §2's arithmetic, and
    # guard_governance imports the scratchpad test.
    "tools/count_loc.py",
    "tools/scratchpad.py",
    # h16 G4: a tools/ file named after a stdlib module shadows it for any
    # script run as `python3 tools/x.py`.
    "tools/ast.py",
    "tools/json/__init__.py",
    "tools/subprocess.cpython-314-darwin.so",
    # h16 G1, review round 7: git 2.55 still reads the legacy remote files, so a
    # plain `git fetch s` writes the named refs they give.
    ".git/remotes/s",
    ".git/branches/s",
)

#: Files named `config` that are not a git config: GOVERNED_SUFFIXES is matched
#: with endswith, never on the basename.
NOT_GOVERNED = (
    "src/app/config",
    "docs/config",
    # h16 G4: only a component after `tools/` is judged by the stdlib's names,
    # and a tools/ file that is not named after one stays ordinary.
    "docs/x/ast.py",
    "tools/scratch_copy.py",
    # h16 G1, review round 7: only the legacy remote files under `.git/` are.
    "docs/remotes/s",
)

#: T7: Bash commands that write the harness state or run away.py.
DENIED_COMMANDS = (
    "echo {} > .git/harness/unattended.json",
    "rm .git/harness/unattended.json",
    "python3 tools/away.py 8h",
    "tools/away.py --back",
    "script -q /dev/null python3 tools/away.py 8h",
    # Review 1: a path to the interpreter, a wrapper, or a reader chained to a
    # run is still a run, not a read.
    ".venv/bin/python3 tools/away.py --back",
    "env python3 tools/away.py 8h",
    "python3 -m ruff check x && python3 tools/away.py --back",
)

#: T7: reads of the same paths, which must stay silent.
SILENT_COMMANDS = (
    "cat .git/harness/queue.jsonl",
    "cat tools/away.py",
    "git diff tools/away.py",
    "pytest tests/python/test_away.py",
    # Review 1: a reader is recognised by the basename of its first token, or
    # as `python -m <reader>`, so the gates the green and review steps run on
    # away.py are not refused.
    "../../../.venv/bin/mypy tools/away.py",
    "/usr/bin/grep -n x tools/away.py",
    "python3 -m ruff check tools/away.py",
    "python -m mypy tools/away.py",
)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


def load_hook(name: str) -> ModuleType:
    """Import `.claude/hooks/<name>.py` from the real checkout, for a unit test of one function."""
    spec = importlib.util.spec_from_file_location(
        f"h16_{name}", REAL / ".claude" / "hooks" / f"{name}.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
    # h4 §3: the refusal names the file it judged, not a generic "a rule file changes".
    assert f"it writes {path}" in reason
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"], line["why"]) == (
        "guard_governance",
        f"{tool} {path}",
        f"it writes {path}",
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


# ---------------------------------------------------------------- h16 G4
#
# docs/increments/h16-harness-fixes.md §2 G4, part 2: a path whose component
# after a `tools` component, cut at the first `.`, is a stdlib module name is
# governed, wherever the checkout is.


@pytest.mark.parametrize("command", ["echo x > tools/ast.py", "cp a.py tools/typing.py"])
def test_a_shell_write_of_a_stdlib_named_tools_file_asks(repo: Path, command: str) -> None:
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert "This changes a file that states rules: tools/" in reason


def test_a_stdlib_named_tools_file_in_another_checkout_asks(repo: Path) -> None:
    path = "/abs/wt/tools/typing.py"
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, file_event(repo, "Write", path)))
    assert found is not None, f"{path} is not governed"
    assert found[0] == "ask"
    assert path in found[1]


# ---------------------------------------------------------------- h16 G3a
#
# §2 G3a: a path whose real path lies under a session scratchpad is not
# governed, whatever its name; a symlink out of the scratchpad is followed.


@pytest.fixture
def pad(repo: Path, tmp_path: Path) -> Path:
    """The copy's scratchpad, holding a git repository `r` and a symlink `link/.git` to `repo`'s."""
    pad = point_scratchpad(repo, tmp_path / "faketmp")
    (pad / "r" / ".git").mkdir(parents=True)
    (pad / "r" / ".git" / "config").write_text("[core]\n")
    (pad / "link").mkdir()
    (pad / "link" / ".git").symlink_to(repo / ".git")
    return pad


#: Files in the scratchpad that would be governed anywhere else.
IN_SCRATCHPAD = (
    "r/.git/config",
    "copy/CLAUDE.md",
    "copy/.claude/agents/tester.md",
    "copy/.claude/hooks/guard_push.py",
    "copy/tools/shell_scan.py",
    "copy/tools/ast.py",  # G3a's exemption comes before G4's stdlib names
)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("tool", ["Edit", "Write"])
@pytest.mark.parametrize("relative", IN_SCRATCHPAD)
def test_a_file_under_a_scratchpad_is_not_governed(
    repo: Path, pad: Path, mode: str, tool: str, relative: str
) -> None:
    set_mode(repo, mode)
    event = file_event(repo, tool, str(pad / relative))
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
def test_a_shell_write_under_a_scratchpad_is_not_governed(repo: Path, pad: Path, mode: str) -> None:
    set_mode(repo, mode)
    command = f"cat >> {pad}/r/.git/config <<'EOF'\n[url \"x\"]\n\tinsteadOf = y\nEOF"
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize(
    "where",
    ["link/.git/config", "relative"],
    ids=["symlink-to-a-real-git-dir", "relative-git-config"],
)
def test_a_git_config_that_resolves_outside_the_scratchpad_asks(
    repo: Path, pad: Path, where: str
) -> None:
    path = ".git/config" if where == "relative" else str(pad / where)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, file_event(repo, "Write", path)))
    assert found is not None, f"{path} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert path in reason


# ---------------------------------------------------------------- h16 G3a, option C
#
# §2 G3a, amended after review round 6 (Ola's option C): for a Bash line the
# exemption holds only when the line is one plain command (one simple command,
# no interpreter program or script, nothing the parser cannot read). The guard
# judges targets before the line runs, so a second command could make a link
# and then write through it; any line with more than one command is judged as
# if there were no exemption.

#: Lines that are not one plain command, each writing a governed name under the
#: scratchpad; `{target}` is the path the ask must name.
NOT_PLAIN: dict[str, str] = {
    # Finding 2's route: `{pad}/x` does not exist when the hook runs.
    "ln -sfn {repo} {pad}/x && echo hi > {pad}/x/CLAUDE.md": "{pad}/x/CLAUDE.md",
    "echo hi > {pad}/copy/CLAUDE.md; true": "{pad}/copy/CLAUDE.md",
    "printf x | tee {pad}/copy/CLAUDE.md": "{pad}/copy/CLAUDE.md",
    "echo $(true) > {pad}/copy/CLAUDE.md": "{pad}/copy/CLAUDE.md",
    "sh -c 'echo hi > {pad}/copy/CLAUDE.md'": "{pad}/copy/CLAUDE.md",
    "python3 -c \"open('{pad}/copy/CLAUDE.md', 'w')\"": "{pad}/copy/CLAUDE.md",
}

#: One plain command each: the exemption holds.
PLAIN = (
    "echo hi > {pad}/copy/CLAUDE.md",
    "cp {repo}/notes.txt {pad}/copy/CLAUDE.md",
)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("template", NOT_PLAIN)
def test_a_scratchpad_write_in_a_line_of_more_than_one_command_asks(
    repo: Path, pad: Path, mode: str, template: str
) -> None:
    set_mode(repo, mode)
    command = template.format(repo=repo, pad=pad)
    target = NOT_PLAIN[template].format(pad=pad)
    assert not (pad / "x").exists()
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    if mode == "on":
        assert kind == "deny"
        assert f"it writes {target}" in reason
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_governance", command)
    else:
        assert kind == "ask"
        assert f"This changes a file that states rules: {target}" in reason
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("template", PLAIN)
def test_a_scratchpad_write_in_one_plain_command_is_not_governed(
    repo: Path, pad: Path, mode: str, template: str
) -> None:
    set_mode(repo, mode)
    command = template.format(repo=repo, pad=pad)
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


def test_governed_without_the_scratchpad_exemption_judges_by_the_rules() -> None:
    """`governed(p, scratch_exempt=False)` skips G3a; the default keeps it (§2 G3a's interface)."""
    hook = load_hook("guard_governance")
    path = "/private/tmp/claude-501/-Users-x-project/0f1e2d3c-session/scratchpad/copy/CLAUDE.md"
    assert hook.governed(path) is False
    assert hook.governed(path, scratch_exempt=False) is True


# ---------------------------------------------------------------- h16 G1, round 7
#
# §2 G1, amendment after review round 7: `GOVERNED_PREFIXES` gains
# `.git/remotes/` and `.git/branches/`, as `.git/hooks/` is there.


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        (".git/remotes/s", True),
        (".git/branches/s", True),
        ("/abs/repo/.git/remotes/origin", True),
        ("docs/remotes/s", False),
    ],
)
def test_the_legacy_remote_files_are_governed(path: str, expected: bool) -> None:
    assert load_hook("guard_governance").governed(path) is expected


# ---------------------------------------------------------------- h16 G7, round 8
#
# §2 G7, amended after review round 8: the widened shell set reaches this guard
# too, since a shell's `-c` program is parsed as a nested line.
#
# ---------------------------------------------------------------- h16 G8
#
# §2 G8: a copy or move into a directory named without a trailing `/` is judged
# both as the name itself and as `<dir>/<basename of each source>`; a
# `-t`/`--target-directory` option makes every operand a source written into it.

#: {command: the governed path the ask must name}.
G8_ASKED: dict[str, str] = {
    "dash -c 'cp x CLAUDE.md'": "CLAUDE.md",
    "ksh -c 'echo x > CLAUDE.md'": "CLAUDE.md",
    "cp json.py tools": "tools/json.py",
    "mv ast.py tools": "tools/ast.py",
    "cp a.py b/json.py tools": "tools/json.py",
    "cp s .git/remotes": ".git/remotes/s",
    "cp s .git/branches": ".git/branches/s",
    "cp x.py .claude/hooks": ".claude/hooks/x.py",
    "cp a.md .claude/agents": ".claude/agents/a.md",
    "cp settings.json .claude": ".claude/settings.json",
    "cp config .git": ".git/config",
    "ln -s /x/y.py .claude/hooks": ".claude/hooks/y.py",
    "install x.py .claude/hooks": ".claude/hooks/x.py",
    "cp -t .claude/hooks x.py": ".claude/hooks/x.py",
    "cp -t.claude/hooks x.py": ".claude/hooks/x.py",
    "cp --target-directory .claude/hooks x.py": ".claude/hooks/x.py",
    "cp --target-directory=.claude/hooks x.py": ".claude/hooks/x.py",
    "mv -t tools ast.py": "tools/ast.py",
    # Controls that ask today and still must.
    "cp notes.txt CLAUDE.md": "CLAUDE.md",
    "cp x.py .claude/hooks/": ".claude/hooks/x.py",
    # Pinned false positive: `backup` meant as a new file, read as a directory.
    "cp CLAUDE.md backup": "backup/CLAUDE.md",
}

#: A last operand with an extension is a file; `src`, `/tmp/inc` and `tools`
#: as directories give names that are not governed.
G8_PASSED = (
    "cp CLAUDE.md /tmp/x.md",
    "cp notes.txt notes.bak",
    "cp json.py src",
    "cp x.py tools",
    "mv ast.py tools/ast_helpers.py",
    "mv tools/old.py tools/new.py",
    "cp -r docs/increments /tmp/inc",
)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G8_ASKED)
def test_a_copy_into_a_directory_or_a_shell_write_asks_naming_the_path(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    if mode == "on":
        assert kind == "deny"
        assert "it writes " in reason and G8_ASKED[command] in reason
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_governance", command)
    else:
        assert kind == "ask"
        assert "This changes a file that states rules: " in reason
        assert G8_ASKED[command] in reason
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G8_PASSED)
def test_a_copy_to_a_file_or_an_ungoverned_directory_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h16 G4, round 9
#
# §2 G4, amended after review round 9 (Ola's option A): `governed()` judges a
# path also in its normalised form (`posixpath.normpath`: `.`, `..` and a
# doubled `/` resolved as text), and in its form as written, so a trailing `/`
# or `/.` that the prefix rules rely on still counts.

RULES = "This changes a file that states rules: "

#: {command: text the ask must name after its lead-in}: the governed file's
#: basename, since the form of the path named (as written or normalised) is
#: not fixed by the design.
G4_ROUND9_ASKED: dict[str, str] = {
    "cp json.py tools/.": "json.py",
    "cp json.py tools/./": "json.py",
    "mv json.py tools/./json.py": "json.py",
    "cp json.py tools//json.py": "json.py",
    "cp json.py tools/x/../json.py": "json.py",
    "cp json.py ./tools/.": "json.py",
    "cp -t tools/. json.py": "json.py",
    "install json.py tools/.": "json.py",
    "echo x > tools/./ast.py": "ast.py",
    "echo x > .claude/./agents/x.md": "x.md",
    "cp a.md .claude/./agents/a.md": "a.md",
    # Controls that ask today by the `.claude/hooks/` prefix on the form as
    # written, and would pass if only the normalised form were judged.
    "rm -r .claude/hooks/": ".claude/hooks/",
    "rm -r .claude/hooks/.": ".claude/hooks/",
}

#: Normalised, each names a file no rule governs.
G4_ROUND9_PASSED = (
    "cp x.py tools/.",
    "cp json.py src/.",
    "mv tools/./old.py tools/./new.py",
    "cp json.py tools/../json.py",
    "cp CLAUDE.md /tmp/x.md",
)

#: Written into a checkout through Write or Edit.
G4_ROUND9_FILES = (".claude/./agents/x.md", "tools/./json.py")


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G4_ROUND9_ASKED)
def test_a_governed_path_written_with_dot_segments_asks(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    lead = "it writes " if mode == "on" else RULES
    assert kind == ("deny" if mode == "on" else "ask")
    assert lead in reason
    assert G4_ROUND9_ASKED[command] in reason.split(lead, 1)[1]
    if mode == "on":
        [line] = queue_lines(repo)
        assert (line["hook"], line["act"]) == ("guard_governance", command)
    else:
        assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", G4_ROUND9_PASSED)
def test_an_ungoverned_path_written_with_dot_segments_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("tool", ["Edit", "Write"])
@pytest.mark.parametrize("relative", G4_ROUND9_FILES)
def test_a_file_tool_write_to_a_governed_path_with_dot_segments_asks(
    repo: Path, mode: str, tool: str, relative: str
) -> None:
    set_mode(repo, mode)
    path = f"{repo}/{relative}"
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, file_event(repo, tool, path)))
    assert found is not None, f"{path} is not governed"
    kind, reason = found
    assert kind == ("deny" if mode == "on" else "ask")
    assert path in reason
    assert len(queue_lines(repo)) == (1 if mode == "on" else 0)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("tool", ["Edit", "Write"])
def test_a_dot_segment_path_under_a_scratchpad_is_not_governed(
    repo: Path, pad: Path, mode: str, tool: str
) -> None:
    """G3a's scratchpad check stays first."""
    set_mode(repo, mode)
    event = file_event(repo, tool, f"{pad}/tools/./json.py")
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("tools/./json.py", True),
        ("tools//ast.py", True),
        ("a/tools/x/../typing.py", True),
        ("tools/./scratch_copy.py", False),
        ("docs/x/./ast.py", False),
        ("tools/../json.py", False),
    ],
)
def test_governed_judges_the_normalised_path_too(path: str, expected: bool) -> None:
    assert load_hook("guard_governance").governed(path) is expected
