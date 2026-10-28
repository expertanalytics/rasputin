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
    # h20 §2 row 6: guard_push and record_review import these two, so an edit
    # to either is an edit to a hook.
    "tools/push_grant.py",
    "tools/ci_changes.py",
    # h31 §10: session_state.py imports it, so it is harness code.
    "tools/work_pipes.py",
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


# ---------------------------------------------------------------- h24 item 5b
#
# docs/increments/h24-harness-fixes.md §3 item 5b, at the governance guard: a
# Python program that only reads a governed file, or the harness state, by
# one of the exempt shapes writes nothing, so the guard is silent. A literal
# that can reach a write still counts, as before.

SETTINGS_READ = "python3 -c 'print(open(\".claude/settings.json\").read())'"
SETTINGS_REUSED = 'python3 -c \'n = open(".claude/settings.json").name; open(n, "w").write("x")\''


def state_read(repo: Path) -> str:
    """Section 3 item 5's reproduction: a json.load of the flag, by its absolute path."""
    flag = repo / ".git" / "harness" / "unattended.json"
    return f"python3 -c 'import json; print(json.load(open(\"{flag}\")))'"


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h24_5b_a_read_of_a_governed_file_is_silent(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    event = bash_event(repo, SETTINGS_READ)
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []


# A Python program is not on the short list of programs a line naming the
# record folder may run, so even a read of the state by one is denied; `cat`
# and `rg` read it.


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h24_5b_a_json_load_of_the_harness_state_is_denied(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    event = bash_event(repo, state_read(repo))
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


def test_h24_5b_a_governed_name_reused_for_a_write_still_asks(repo: Path) -> None:
    set_mode(repo, "off")
    event = bash_event(repo, SETTINGS_REUSED)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found is not None, "the write through a reused name passed silently"
    kind, reason = found
    assert kind == "ask"
    assert "This changes a file that states rules: " in reason
    assert ".claude/settings.json" in reason


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h24_5b_a_read_of_the_harness_state_with_an_opener_is_denied(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    flag = repo / ".git" / "harness" / "unattended.json"
    program = f'def W(*a):\n    return 0\nprint(open("{flag}", opener=W).read())'
    event = bash_event(repo, f"python3 -c '{program}'")
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


# ------------------------------------------------- h24 item 5b, code review 1
#
# The scanner reads the program text before the shell expands it. In a
# double-quoted `-c` or an unquoted heredoc, `$(<q)` is replaced
# by the shell first: with q holding `','w').write('WRITTEN'))#`, the read
# `open('<path>$(<q)').read()` runs as `open('<path>','w').write(...)`. The fix
# proposed in review: no exemption when the program text holds `$` or a
# backtick. Single-quoted `-c` text and a quoted heredoc are not expanded.

#: The file the expansion reads; its text turns the read into a write.
SPLICE = "','w').write('WRITTEN'))#"


def expanding_reads(target: str) -> dict[str, str]:
    """A read of `target` whose program text the shell rewrites before Python sees it."""
    return {
        "double-quoted-c": f"python3 -c \"print(open('{target}$(<q)').read())\"",
        "unquoted-heredoc": f"python3 - <<EOF\nprint(open('{target}$(<q)').read())\nEOF",
    }


def quoted_reads(target: str) -> dict[str, str]:
    """The same read with no expansion: single-quoted `-c` text and a quoted heredoc."""
    program = f'print(open("{target}").read())'
    return {
        "single-quoted-c": f"python3 -c '{program}'",
        "quoted-heredoc": f"python3 - <<'EOF'\n{program}\nEOF",
    }


EXPANSIONS = list(expanding_reads("T"))
QUOTED = list(quoted_reads("T"))


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("form", EXPANSIONS)
def test_h24_5b_an_expanded_read_of_the_harness_state_is_denied(
    repo: Path, mode: str, form: str
) -> None:
    set_mode(repo, mode)
    (repo / "q").write_text(SPLICE)
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, expanding_reads(str(flag))[form])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


@pytest.mark.parametrize("form", EXPANSIONS)
def test_h24_5b_an_expanded_read_of_a_governed_file_asks(repo: Path, form: str) -> None:
    set_mode(repo, "off")
    (repo / "q").write_text(SPLICE)
    event = bash_event(repo, expanding_reads(".claude/settings.json")[form])
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found is not None, "a read the shell turns into a write passed silently"
    kind, reason = found
    assert kind == "ask"
    assert f"{RULES}.claude/settings.json" in reason, reason


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("form", QUOTED)
def test_h24_5b_an_unexpanded_read_of_the_harness_state_is_denied(
    repo: Path, mode: str, form: str
) -> None:
    set_mode(repo, mode)
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, quoted_reads(str(flag))[form])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


@pytest.mark.parametrize("form", QUOTED)
def test_h24_5b_an_unexpanded_read_of_a_governed_file_is_silent(repo: Path, form: str) -> None:
    set_mode(repo, "off")
    event = bash_event(repo, quoted_reads(".claude/settings.json")[form])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h24 item 14
#
# docs/increments/h24-harness-fixes.md §3 item 14: where the shell may rewrite
# a Python program before Python reads it, the guard cannot know what runs. A
# literal in a program holding `$` or a backtick is a candidate whatever its
# spaces, judged by its static form; and the interpreter is `unknown`, so the
# line falls to the deny rule (it names `.git/harness`) or the text rule (it
# names a governed file). In an unquoted heredoc a backslash counts too: the
# shell turns `\\` into `\`, so the scanner's comment is Python's write.


def backtick_reads(target: str) -> dict[str, str]:
    """A read of `target` with a backtick substitution spliced into the literal."""
    return {
        "double-quoted-c": f"python3 -c \"print(open('{target}`cat q`').read())\"",
        "unquoted-heredoc": f"python3 - <<EOF\nprint(open('{target}`cat q`').read())\nEOF",
    }


BACKTICKS = list(backtick_reads("T"))


def backslash_body(target: str) -> str:
    """Design review round 3's B(t): as written a comment hides the write; after the
    shell's `\\` -> `\\` in an unquoted heredoc, Python runs it."""
    return rf"x = 'a\\'; y = 1 #'; open('{target}','w').write('X')"


def heredoc(body: str, opener: str) -> str:
    return f"python3 - {opener}\n{body}\nEOF"


def asked_naming(repo: Path, command: str, name: str) -> None:
    """The hook asks, by day, and names `name` in its reason."""
    set_mode(repo, "off")
    (repo / "q").write_text(SPLICE)
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found is not None, f"{command!r} passed silently"
    kind, reason = found
    assert kind == "ask"
    assert RULES in reason and name in reason, reason


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("form", BACKTICKS)
def test_h24_14_a_backtick_read_of_the_harness_state_is_denied(
    repo: Path, mode: str, form: str
) -> None:
    set_mode(repo, mode)
    (repo / "q").write_text(SPLICE)
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, backtick_reads(str(flag))[form])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


SUBSTITUTED_GOVERNED: dict[str, tuple[str, str]] = {
    "backtick-read-of-settings": (
        "python3 -c \"print(open('.claude/settings.json`cat q`').read())\"",
        ".claude/settings.json"),
    "dollar-cat-read-of-settings": (
        "python3 -c \"print(open('.claude/settings.json$(cat q)').read())\"",
        ".claude/settings.json"),
    "a-write-whose-path-is-all-substitution": (
        "python3 -c \"open('$(echo CLAUDE.md)','w').write('x')\"", "CLAUDE.md"),
    "backslash-in-an-unquoted-heredoc": (
        heredoc(backslash_body("CLAUDE.md"), "<<EOF"), "CLAUDE.md"),
}  # fmt: skip


@pytest.mark.parametrize(
    ("command", "name"), SUBSTITUTED_GOVERNED.values(), ids=SUBSTITUTED_GOVERNED.keys()
)
def test_h24_14_a_program_the_shell_rewrites_naming_a_governed_file_asks(
    repo: Path, command: str, name: str
) -> None:
    asked_naming(repo, command, name)


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h24_14_a_backslash_in_an_unquoted_heredoc_writing_the_harness_state_is_denied(
    repo: Path, mode: str
) -> None:
    set_mode(repo, mode)
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, heredoc(backslash_body(str(flag)), "<<EOF"))
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


def item_14_silent(flag: Path) -> dict[str, str]:
    """Lines item 14 must leave silent: nothing substituted, or a quoted heredoc."""
    return {
        # Quoted: Python too reads the write as a comment.
        "backslash-in-a-quoted-heredoc": heredoc(backslash_body("CLAUDE.md"), "<<'EOF'"),
        "a-dollar-in-a-prose-literal": "python3 -c 'print(\"a $ b\")'",
        "a-governed-name-in-prose": "python3 -c 'print(\"see CLAUDE.md for rules\")'",
        "a-dollar-beside-a-plain-read": (
            "python3 -c 'print(open(\"/x/notes.txt\").read().count(\"$\"))'"),
    }  # fmt: skip


def item_14_state_reads(flag: Path) -> dict[str, str]:
    """Python reads of the harness state by a quoted heredoc: denied, as no Python
    program is on the record folder's short list."""
    return {
        "a-json-read-by-quoted-heredoc": heredoc(
            rf'import json; print("\n".join(json.load(open("{flag}"))))', "<<'EOF'"),
        "a-dollar-read-by-quoted-heredoc": heredoc(
            f'print(open("{flag}").read().count("$"))', "<<'EOF'"),
    }  # fmt: skip


@pytest.mark.parametrize("label", list(item_14_state_reads(Path("/x"))))
def test_h24_14_a_python_read_of_the_harness_state_is_denied(repo: Path, label: str) -> None:
    set_mode(repo, "off")
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, item_14_state_reads(flag)[label])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


ITEM_14_SILENT = list(item_14_silent(Path("/x")))


@pytest.mark.parametrize("label", ITEM_14_SILENT)
def test_h24_14_a_line_nothing_rewrites_or_that_only_reads_is_silent(
    repo: Path, label: str
) -> None:
    set_mode(repo, "off")
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, item_14_silent(flag)[label])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h24 item 15
#
# docs/increments/h24-harness-fixes.md §3 item 15: the shell hands a
# here-string (`<<< word`), plus a newline, to the command's stdin, so for an
# interpreter or shell with no inline program it is the program. zsh joins
# every stdin redirect of one command (MULTIOS), so each body counts. A
# here-string holding `$` or a backtick is unknown whatever its quoting.

#: A Python write of CLAUDE.md, double-quoted for the shell.
WRITE_CLAUDE = "\"open('CLAUDE.md','w').write('x')\""

HERE_STRING_ASKED: dict[str, tuple[str, str]] = {
    "python3-double-quoted": (f"python3 <<< {WRITE_CLAUDE}", "CLAUDE.md"),
    "python3-single-quoted": ("python3 <<< 'open(\"CLAUDE.md\",\"w\").write(\"x\")'", "CLAUDE.md"),
    "bash": ('bash <<< "echo x > CLAUDE.md"', "CLAUDE.md"),
    "zsh": ("zsh <<< 'echo x > .claude/settings.json'", ".claude/settings.json"),
    "perl": ("perl <<< 'open(F, \">CLAUDE.md\")'", "CLAUDE.md"),
    "env": (f"env python3 <<< {WRITE_CLAUDE}", "CLAUDE.md"),
    "env-with-assignment": (f"env A=1 python3 <<< {WRITE_CLAUDE}", "CLAUDE.md"),
    "nice": (f"nice -n 5 python3 <<< {WRITE_CLAUDE}", "CLAUDE.md"),
    "timeout": (f"timeout 5 python3 <<< {WRITE_CLAUDE}", "CLAUDE.md"),
    "timeout-bash": ('timeout 5 bash <<< "echo x > CLAUDE.md"', "CLAUDE.md"),
    # Two stdin redirects on one command: zsh runs both, so each counts.
    "join-write-second": (f"python3 <<< 'print(1)' <<< {WRITE_CLAUDE}", "CLAUDE.md"),
    "join-write-first": (f"python3 <<< {WRITE_CLAUDE} <<< 'print(2)'", "CLAUDE.md"),
    "here-string-write-beside-a-heredoc": (
        f"python3 <<< {WRITE_CLAUDE} <<EOF\nprint(1)\nEOF", "CLAUDE.md"),
    # Unchanged: a heredoc's write beside a here-string still counts.
    "heredoc-write-beside-a-here-string": (
        "python3 <<< 'print(1)' <<EOF\nopen('CLAUDE.md','w').write('x')\nEOF", "CLAUDE.md"),
    # `$` or a backtick: unknown whatever the quoting, so the line's text rule
    # judges it; the governed name sits only in a comment, outside any literal.
    "dollar-single-quoted": ("python3 <<< 'x = \"$y\"  # CLAUDE.md'", "CLAUDE.md"),
    "dollar-double-quoted": ("python3 <<< \"x = '$(cat q)'  # CLAUDE.md\"", "CLAUDE.md"),
    "backtick-single-quoted": ("python3 <<< 'x = \"`cat q`\"  # CLAUDE.md'", "CLAUDE.md"),
    # Item 14's write whose path is all substitution, by here-string.
    "a-write-whose-path-is-all-substitution": (
        "python3 <<< \"open('$(echo CLAUDE.md)','w').write('x')\"", "CLAUDE.md"),
    # A shell parses every stdin body whatever its -c program.
    "shell-c-program-that-is-a-shell": (
        "bash -c 'bash' <<'EOF'\necho x > CLAUDE.md\nEOF", "CLAUDE.md"),
    # A heredoc on a shell that has any argument feeds what the script runs
    # (docs/increments/h25-stdin-routes.md §3.2, line feeds).
    "shell-script-beside-a-heredoc": (
        "bash script.sh <<'EOF'\necho x > CLAUDE.md\nEOF", "CLAUDE.md"),
}  # fmt: skip


@pytest.mark.parametrize(
    ("command", "name"), HERE_STRING_ASKED.values(), ids=HERE_STRING_ASKED.keys()
)
def test_h24_15_a_here_string_program_naming_a_governed_file_asks(
    repo: Path, command: str, name: str
) -> None:
    asked_naming(repo, command, name)
    assert queue_lines(repo) == []


def here_string_denied(flag: Path) -> dict[str, str]:
    """Here-string programs that write, or may write, the harness state."""
    return {
        "python3-write": f"python3 <<< \"open('{flag}','w').write('x')\"",
        "bash-write": f"bash <<< 'echo x > {flag}'",
        "python3-backtick-read": f"python3 <<< \"print(open('{flag}`cat q`').read())\"",
    }


HERE_STRING_DENIED = list(here_string_denied(Path("/x")))


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("label", HERE_STRING_DENIED)
def test_h24_15_a_here_string_program_writing_the_harness_state_is_denied(
    repo: Path, mode: str, label: str
) -> None:
    set_mode(repo, mode)
    (repo / "q").write_text(SPLICE)
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, here_string_denied(flag)[label])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


def here_string_silent(flag: Path) -> dict[str, str]:
    """Here-strings item 15 must leave silent: data, an ordinary file, or a read."""
    return {
        "cat-plain": "cat <<< x",
        "cat-of-a-python-write": "cat <<< 'open(\"CLAUDE.md\",\"w\")'",
        "grep-of-a-governed-name": "grep -c x <<< CLAUDE.md",
        "python3-writes-an-ordinary-file": "python3 <<< \"open('notes/a.txt','w').write('x')\"",
        "bash-writes-an-ordinary-file": 'bash <<< "echo x > notes/a.txt"',
        "two-here-strings-that-print": "python3 <<< 'print(1)' <<< 'print(2)'",
        # A -c program keeps stdin as data.
        "c-program-beside-a-here-string": "python3 -c 'print(1)' <<< \"open('CLAUDE.md','w')\"",
        # A data reader's inline program keeps a heredoc as data.
        "c-program-beside-a-heredoc": (
            "python3 -c 'print(1)' <<'EOF'\nopen('CLAUDE.md','w').write('x')\nEOF"),
    }  # fmt: skip


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h24_15_a_here_string_json_read_of_the_harness_state_is_denied(
    repo: Path, mode: str
) -> None:
    """No Python program is on the record folder's short list, a read included."""
    set_mode(repo, mode)
    flag = repo / ".git" / "harness" / "unattended.json"
    command = f"python3 <<< 'import json; print(json.load(open(\"{flag}\")))'"
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


HERE_STRING_SILENT = list(here_string_silent(Path("/x")))


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("label", HERE_STRING_SILENT)
def test_h24_15_a_here_string_that_writes_nothing_governed_is_silent(
    repo: Path, mode: str, label: str
) -> None:
    set_mode(repo, mode)
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, here_string_silent(flag)[label])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- h25
#
# docs/increments/h25-stdin-routes.md §3.2: a runner (an interpreter, a shell,
# or a command that could be one) whose stdin or script could carry text the
# scanner did not read is unread, and so unknown: the line is judged by its
# text. A Python data reader (§3.4, an allowlist) keeps its stdin as data.

#: A Python write of CLAUDE.md, the program fed on stdin.
W = "open('CLAUDE.md','w').write('x')"
#: A shell write of CLAUDE.md, the commands fed on stdin.
SH_W = "echo x > CLAUDE.md"
#: A Python program a data reader may run: only json and sys, read only.
READ_STDIN = "import sys; print(sys.stdin.read())"

H25_REVIEW_ESCAPES = (
    f'echo "{W}" |\npython3',
    f'echo "{W}" | {{ true; python3; }}',
    f'echo "{W}" | (true; python3)',
    f'echo "{W}" | if true; then python3; fi',
    f'echo "{W}" | bash -c python3',
    f'bash -c python3 <<< "{W}"',
    f"echo '{SH_W}' | bash -s foo",
    f"bash -s foo <<< '{SH_W}'",
    f"echo '{SH_W}' | bash -o errexit",
    f'echo "{W}" | sudo -E python3 -',
    f'python3 -i -c pass <<< "{W}"',
    f'echo "{W}" | python3 -i -c pass',
    f'python3 =(echo "{W}")',
    f'python3 3< <(echo "{W}") < /dev/fd/3',
    f'python3 3< <(echo "{W}") <&3',
    f'echo "{W}" | xargs -I{{}} python3 -c {{}}',
    f"echo \"{W}\" | env -S 'python3 -'",
    f"echo '{SH_W}' | $SHELL",
    f'python3 //dev/stdin <<< "{W}"',
    f'exec 0< <(echo "{W}"); python3',
    f"echo \"{W}\" | python3 -c 'exec(input())'",
    f"python3 -c 'exec(input())' <<< \"{W}\"",
    f'echo "{W}" | python3 -m cProfile /dev/stdin',
    f'echo "{W}" | arch -arm64 python3',
)

H25_FURTHER_ROUTES = (
    f'echo "{W}" |  # c\npython3',
    f'echo "{W}" | (python3)',
    f'echo "{W}" |& python3',
    f'echo "{W}" | python3 -m code',
    f'echo "{W}" | uv run -',
    f'echo "{W}" | uv run python3',
    f'uv run python3 - <<< "{W}"',
    f'echo "{W}" | tee >(python3)',
    f'python3 <(echo "{W}")',
    f"source <(echo '{SH_W}')",
    f"echo '{SH_W}' | while read l; do eval \"$l\"; done",
    f'echo "{W}" | echo "$(python3)"',
    f'echo "{W}" | {{ cat; }} | python3',
    f'f() {{ python3; }}; echo "{W}" | f',
    f"function g {{ sh; }}; echo '{SH_W}' | g",
    f'coproc python3; echo "{W}" >&"${{COPROC[1]}}"',
    f"$SHELL <<< '{SH_W}'",
    f"bash -c 'echo \"{W}\" | python3'",
    f"echo \"{W}\" | sh -c 'cat | python3'",
    f'{{ python3; }} <<< "{W}"',
    f'{{ python3; }} < <(echo "{W}")',
    f'( python3 ) < <(echo "{W}")',
    f'while read l; do python3; done < <(echo "{W}")',
    f"echo '{SH_W}' | parallel",
    'cat CLAUDE.md | python3 -c "import os,sys; os.system(sys.stdin.read())"',
    # Asks today too: a process substitution on an output redirect.
    f'echo "{W}" > >(python3)',
)

#: Python programs that run their stdin, each double-quoted for the shell.
H25_ROUND2_PROGRAMS = (
    "import sys; sys.modules['posix'].system(sys.stdin.read())",
    "import sys; __loader__.load_module('posix').system(sys.stdin.read())",
    "import operator,sys; operator.methodcaller('system', sys.stdin.read())(sys.modules['posix'])",
    "import string,sys; string.Formatter().get_field('0.system', [sys.modules['posix']], {})[0]"
    "(sys.stdin.read())",
    "import sys; g=(x for x in [1]); g.gi_frame.f_builtins['exec'](sys.stdin.read())",
    "import sys; [0 for exec in []]; exec(sys.stdin.read())",
)

#: Design review round 2's routes. Each program passes a list of what is
#: dangerous and fails the allowlist of §3.4; a bound builtin name in a
#: comprehension leaves the builtin in place.
H25_ROUND2 = (
    *(f'echo "{W}" | python3 -c "{program}"' for program in H25_ROUND2_PROGRAMS),
    f"echo '{SH_W}' | /bin/s?",
    f"echo '{SH_W}' | =sh",
    f"read -r l <<< '{SH_W}'; eval \"$l\"",
    f"l=$(cat <<< '{SH_W}'); eval \"$l\"",
    f'echo "{W}" | env -S"python3 -"',
    f"echo \"{W}\" | python3 -W ignore -c '{READ_STDIN}'",
)  # fmt: skip


@pytest.mark.parametrize("command", [*H25_REVIEW_ESCAPES, *H25_FURTHER_ROUTES, *H25_ROUND2])
def test_h25_a_runner_fed_text_naming_a_governed_file_asks(repo: Path, command: str) -> None:
    asked_naming(repo, command, "CLAUDE.md")
    assert queue_lines(repo) == []


def h25_denied(flag: Path) -> dict[str, str]:
    """Fed programs that write the harness state."""
    write = f"open('{flag}','w').write('x')"
    return {
        "piped-python3": f'echo "{write}" | python3',
        "shell-with-an-operand": f"bash -s foo <<< 'echo x > {flag}'",
        "compound-fed-by-process-substitution": f'{{ python3; }} < <(echo "{write}")',
    }


H25_DENIED = list(h25_denied(Path("/x")))


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("label", H25_DENIED)
def test_h25_a_fed_runner_writing_the_harness_state_is_denied(
    repo: Path, mode: str, label: str
) -> None:
    set_mode(repo, mode)
    flag = repo / ".git" / "harness" / "unattended.json"
    event = bash_event(repo, h25_denied(flag)[label])
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


H25_SILENT: dict[str, str] = {
    "stdin-from-a-file": "python3 < prog.py",
    "after-or": "echo x || python3",
    "after-and": "echo x && python3",
    "named-script-is-a-data-reader": f'echo "{W}" | python3 x.py',
    "governed-tool-as-a-data-reader": "git diff | python3 tools/count_loc.py",
    "c-data-reader-of-a-governed-file": (
        "cat CLAUDE.md | python3 -c 'import sys; print(len(sys.stdin.read()))'"),
    "json-reader-fed-a-governed-file": (
        "python3 -c 'import json,sys; print(json.load(sys.stdin))' < CLAUDE.md"),
    "heredoc-to-a-named-script": "python3 tools/x.py <<'EOF'\n{\"path\": \"CLAUDE.md\"}\nEOF",
    "stdin-from-dev-null": "python3 x.py < /dev/null",
    "process-substitution-on-a-reader": f'diff <(echo "{W}") x; python3 x.py',
    "json-tool-fed-data": "echo '{}' | python3 -m json.tool",
    "piped-program-naming-nothing-governed": "printf 'print(1)' | python3",
    "source-of-a-file": "source .venv/bin/activate",
    "uv-run-of-a-tool": "uv run pytest tests/python",
    "loop-fed-a-listing": "while read l; do echo $l; done < <(ls CLAUDE.md)",
    # An interpreter's name among the words of a command that runs nothing.
    "grep-for-an-interpreter-name": "grep -i python CLAUDE.md | head",
}  # fmt: skip


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", H25_SILENT.values(), ids=H25_SILENT.keys())
def test_h25_data_on_stdin_or_a_runner_not_fed_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# ------------------------------------------------- h25 code review round 1
#
# docs/increments/h25-stdin-routes.md §3.4: a data reader's script operand is
# "not under /dev/", judged after `.` and `..` are resolved. A path that
# reaches /dev/stdin through them runs stdin, as `python3 /dev/stdin` and
# `python3 -` do (the first two, the controls), so the program is fed.

H25_STDIN_SCRIPTS = (
    "-",
    "/dev/stdin",
    "/./dev/stdin",
    "../../../../../../../../dev/stdin",
    "/usr/../dev/stdin",
    "/usr/lib/../../dev/fd/0",
)


@pytest.mark.parametrize("script", H25_STDIN_SCRIPTS)
def test_h25_a_piped_script_path_resolving_to_stdin_asks(repo: Path, script: str) -> None:
    asked_naming(repo, f'echo "{W}" | python3 {script}', "CLAUDE.md")
    assert queue_lines(repo) == []


@pytest.mark.parametrize("script", H25_STDIN_SCRIPTS)
def test_h25_a_here_string_on_a_script_path_resolving_to_stdin_asks(
    repo: Path, script: str
) -> None:
    asked_naming(repo, f'python3 {script} <<< "{W}"', "CLAUDE.md")
    assert queue_lines(repo) == []


# ------------------------------------------------- h25 code review round 2
#
# docs/increments/h25-stdin-routes.md §3.4: a data reader's script operand is
# "not under /dev/". Each line below runs its stdin as the script, checked on
# macOS under zsh and bash with `print('ran')` on 2026-10-20:
# - the macOS root is case-insensitive, so `/DEV/stdin` and `/Dev/fd/0` are
#   /dev/stdin and /dev/fd/0 (devfs below it is not: `/dev/STDIN` and
#   `/DEV/FD/0` fail, and are not pinned);
# - a relative script on a line that changes directory resolves against a
#   directory the guard cannot know (`cd "$D"`), so it counts as fed.

H25_CASE_SCRIPTS = ("/DEV/stdin", "/Dev/fd/0")

H25_CD_ROUTES: dict[str, str] = {
    "cd-and": f'cd /dev && echo "{W}" | python3 stdin',
    "cd-semicolon-fd": f'cd /dev/fd; echo "{W}" | python3 0',
    "cd-newline": f'cd /dev\necho "{W}" | python3 stdin',
    "cd-variable": f'D=/dev; cd "$D" && echo "{W}" | python3 stdin',
    "pushd-and": f'pushd /dev && echo "{W}" | python3 stdin',
    "pushd-semicolon-fd": f'pushd /dev/fd; echo "{W}" | python3 0',
}  # fmt: skip

#: The controls: a relative script with no cd or pushd on the line is a data
#: reader as today, and a cd before a runner nothing feeds stays silent.
H25_CD_CONTROLS: dict[str, str] = {
    "relative-stdin-no-cd": f'echo "{W}" | python3 stdin',
    "relative-0-no-cd": f'echo "{W}" | python3 0',
    "cd-before-an-unfed-runner": "cd docs && python3 x.py",
}  # fmt: skip


@pytest.mark.parametrize("script", H25_CASE_SCRIPTS)
def test_h25_a_piped_script_path_naming_dev_in_another_case_asks(repo: Path, script: str) -> None:
    asked_naming(repo, f'echo "{W}" | python3 {script}', "CLAUDE.md")
    assert queue_lines(repo) == []


@pytest.mark.parametrize("command", H25_CD_ROUTES.values(), ids=H25_CD_ROUTES.keys())
def test_h25_a_relative_script_on_a_line_that_changes_directory_asks(
    repo: Path, command: str
) -> None:
    asked_naming(repo, command, "CLAUDE.md")
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("command", H25_CD_CONTROLS.values(), ids=H25_CD_CONTROLS.keys())
def test_h25_a_relative_script_with_no_directory_change_is_silent(
    repo: Path, mode: str, command: str
) -> None:
    set_mode(repo, mode)
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- the record folder
#
# A line that names the harness record folder (in any letter case, also after
# quotes are removed), or that changes directory into `.git`, passes only when
# every simple command is on a short list of programs that cannot write, every
# write target is a plain path outside the folder, no `cd` leads into it, and
# no expansion runs code the parser did not read. Every other line is judged
# as before. The folder's name is built from parts here, so a line that writes
# this file never names it.

H = ".git/" + "harness"
DENIED = ("deny", [], "")


def record_folder_passes() -> dict[str, str]:
    """Lines that only mention the folder or read from it, each with an id."""
    return {
        "quoted-heredoc-into-session": (
            "cat >> .claude/current-task/session.md <<'EOF'\n"
            f"NOW: see {H}/x and $(git rev-parse --git-dir)/harness\nEOF"
        ),
        "unquoted-heredoc-into-session": (
            f"cat >> .claude/current-task/session.md <<EOF\nNOW: see {H}/x\nEOF"
        ),
        "echo-into-a-note": f'echo "see {H}" > note.md',
        "printf-appending": f"printf '%s\\n' {H} >> notes.md",
        "commit-m": f'git commit -m "guard: {H}"',
        "commit-qm": f'git commit -qm "guard: {H}"',
        "commit-F": f"echo {H} && git commit -F msg.txt",
        "commit-message-long": f'git commit --message "guard: {H}"',
        "brief-task": f'python3 tools/brief.py a --task "the {H}/x file"',
        "brief-task-to-a-file": f'python3 tools/brief.py a --task "the {H}/x" > /tmp/brief.txt',
        "cat": f"cat {H}/queue.jsonl",
        "grep": f"grep -n x {H}/queue.jsonl",
        "rg": f"rg x {H}",
        "rg-to-a-file": f"rg x {H} > /tmp/hits",
        "wc": f"wc -l {H}/queue.jsonl",
        "cut": f"cut -c1-9 {H}/queue.jsonl",
        "head-cut": f"head -3 {H}/queue.jsonl | cut -c1-9",
        "tail": f"tail -n 2 {H}/queue.jsonl",
        "ls": f"ls -la {H}",
        "stat": f"stat {H}/unattended.json",
        "diff": f"diff {H}/queue.jsonl /tmp/q",
        "git-log": f"git log --oneline && echo {H}",
        "git-show": f"git show HEAD && echo {H}",
        "git-status": f"git status && echo {H}",
        "git-diff": f"git diff --stat && echo {H}",
        "cd-absolute-then-grep": f"cd /abs/checkout && grep -n x {H}/queue.jsonl",
        "cd-git-status-echo": f"cd .git && git status && echo {H}",
        "plain-expansions": f"echo $HOME {H} ${{HOME}} $1",
        "cd-to-toplevel": "cd $(git rev-parse --show-toplevel) && ls",
        "cd-hooks-ls": "cd .git/hooks && ls",
        "upper-case-read": f"cat {H.upper()}/queue.jsonl",
        # A reader on the list reads the folder, its output filtered by readers.
        "grep-review-log": f"grep -c APPROVED {H}/reviews.jsonl | tail -3",
    }


def record_folder_denials() -> dict[str, str]:
    """Lines that write, or may write, the folder, or run code unread; each with an id."""
    lines = {
        "cd-git-touch": "cd .git && touch harness/x",
        "cd-folder-python-c": f"cd {H} && python3 -c \"open('x','w')\"",
        "cd-folder-rsync": f"cd {H} && rsync -a /s/ . > /tmp/log",
        "cd-folder-tar": f"cd {H} && tar xf /tmp/a.tar > /tmp/log",
        "cd-folder-vim": f"cd {H} && vim x > /tmp/log",
        "cd-folder-sed-w": f"cd {H} && sed -n 'w x' in > /tmp/log",
        "cd-folder-file-C": f"cd {H} && file -C -m m",
        "cd-folder-cat": f"cd {H} && cat queue.jsonl",
        "cd-git-rsync": "cd .git && rsync -a /s/ harness/",
        "python-argv-open": f"python3 -c \"import sys; open(sys.argv[1],'w')\" {H}/x",
        "bash-script": f"bash s.sh {H}/x",
        "python-tool": f"python3 tools/x.py {H}/x",
        "python-m": f"python3 -m mod {H}/x",
        "sed-w": f"sed -n 'w {H}/x' in",
        "awk-print": f"awk '{{print > \"{H}/x\"}}' in",
        "git-apply-then-cat": f"git apply p && cat {H}/x",
        "assigned-target": f'N={H}/x; echo x > "$N"',
        "ed-on-a-pipe": f"printf 'a\\nx\\n.\\nw {H}/x\\n' | ed -s",
        "ed-quoted-heredoc": f"ed -s <<'EOF'\na\nx\n.\nw {H}/x\nEOF",
        "clobber": f"echo x >! {H}/y",
        "append-clobber": f"echo x >>! {H}/y",
        "both-clobber": f"echo x &>! {H}/y",
        "append-pipe-clobber": f"echo x >>| {H}/y",
        "glob-question": f"echo {H} > .git/harnes?/q",
        "glob-bracket": f"echo {H} > .git/harnes[s]/q",
        "e-flag": f"echo ${{(e):-'$(touch x)'}} {H}",
        "brace-default": f"echo ${{HOME:-x}} {H}",
        "arithmetic": f"echo $((1+1)) {H}",
        "old-arithmetic": f"echo $[1] {H}",
        "unquoted-heredoc-runs": f"cat > note.md <<EOF\n$(touch {H}/x)\nEOF",
        "unquoted-heredoc-runs-elsewhere": f"cat > note.md <<EOF\n$(touch x) {H}\nEOF",
        "unquoted-heredoc-backtick": f"cat > note.md <<EOF\n`touch x` {H}\nEOF",
        "here-string-quoted": f'cat <<< "$(touch {H}/x)"',
        "here-string-unquoted": f"cat <<< $(touch {H}/x)",
        "glob-qualifier": f"echo {H} *(e:'touch x':)",
        "rg-pre": f"rg --pre ./x pat {H}",
        "rg-pre-equals": f"rg --pre=./x pat {H}",
        "rg-pre-glob": f"rg --pre-glob '*' --pre ./x pat {H}",
        "rg-hostname-bin": f"rg --hostname-bin=./x pat {H}",
        "git-log-output": f"git log --output=/tmp/o && echo {H}",
        "git-log-output-in": f"git log --output={H}/x",
        "git-c": f"git -c core.pager=x log && echo {H}",
        "git-C": f"git -C /tmp/r status && echo {H}",
        "git-dir": f"git --git-dir=/tmp/r status && echo {H}",
        "ripgrep-config": f"RIPGREP_CONFIG_PATH=/tmp/c rg x {H}",
        "pythonpath": f'PYTHONPATH=/tmp python3 tools/brief.py a --task "{H}"',
        "read-write-open": f"cat <> {H}/x",
        "cd-git-echo": "cd .git && echo x > harness/y",
        "cd-hooks-up": "cd .git/hooks && echo x > ../harness/y",
        "cd-elsewhere-brief": f"cd /tmp/e && python3 tools/brief.py a --task {H}/x",
        "pushd-elsewhere-brief": f"pushd /tmp/e && python3 tools/brief.py a --task {H}/x",
        "env-rg": f"env rg x {H}",
        "xargs-touch": f"ls {H} | xargs touch",
        "pipe-to-sh": f"echo touch {H}/x | sh",
        "home-target": f"echo x > ~/r/{H}/y",
        "cd-two-arguments": f"cd a b && echo {H}",
        "cd-option": f"cd -P .git && echo {H}",
        "sort": f"sort {H}/queue.jsonl",
        "jq": f"jq . {H}/queue.jsonl",
        "rm-then-heredoc": f"rm n.md && cat >> session.md <<'EOF'\n{H}\nEOF",
        # A short option cluster holding `e` opens the editor.
        "commit-em": f'git commit -em "x {H}"',
        "commit-e-m": f'git commit -e -m "x {H}"',
        "commit-edit": f'git commit --edit -m "x {H}"',
        "commit-no-message": f"git commit && echo {H}",
        # A commit message spares only the words inside it.
        "commit-m-then-python": f"git commit -m \"{H}\" && python3 -c \"open('{H}/q','w')\"",
        "python-then-commit-m": f"python3 -c \"open('{H}/q','w')\"; git commit -m \"{H}\"",
        "add-then-commit-em": f'git add -A && git commit -em "x {H}"',
        "add-then-commit-e-m": f'git add -A && git commit -e -m "x {H}"',
        "add-then-commit-m-runs": f'git add -A && git commit -m "$(touch x) {H}"',
        "add-then-commit-m-backtick": f'git add -A && git commit -m "`touch x` {H}"',
        "add-then-commit-m-pathspec": f'git add -A && git commit -m "{H}" -- {H}/x',
        "add-then-tag-m": f'git add -A && git tag -m "x {H}" v1',
        "readers-and-python": f"grep x {H}/r.jsonl | cut -c1-9; python3 -c \"open('{H}/q')\"",
        # The disk is case-insensitive: `.GIT/HARNESS` is the folder.
        "upper-case-target": "echo x > " + H.upper() + "/y",
        "mixed-case-cd": "cd .Git && touch Harness/y",
        "mixed-case-mention": "touch x && echo " + ".Git/" + "Harness",
        # The trigger also reads the words with quotes removed.
        "quoted-inside-the-name": "rsync -a s/ .git/har''ness/ > log",
        "double-quoted-name": 'rsync -a s/ .git/"harness"/ > log',
    }
    for program in (
        "rsync -a /s/ {}/",
        "tar -C {} -xf /tmp/a.tar",
        "find {} -delete",
        "ex -sc wq {}/x",
        "vim {}/x",
        "sort -o {}/x y",
        "uniq y {}/x",
    ):
        word = program.split()[0]
        lines[f"{word}-alone"] = program.format(H)
        lines[f"{word}-to-a-log"] = program.format(H) + " > log"
    return lines


@pytest.fixture(scope="module")
def guard() -> ModuleType:
    return load_hook("guard_governance")


@pytest.mark.parametrize("command", record_folder_passes().values(),
                         ids=record_folder_passes().keys())  # fmt: skip
def test_a_line_that_only_mentions_or_reads_the_record_folder_passes(
    guard: ModuleType, command: str
) -> None:
    assert guard.judge_bash(command) is None


@pytest.mark.parametrize("command", record_folder_denials().values(),
                         ids=record_folder_denials().keys())  # fmt: skip
def test_a_line_that_may_write_the_record_folder_is_denied(guard: ModuleType, command: str) -> None:
    assert guard.judge_bash(command) == DENIED


#: The folder named only inside a commit message, beside a command off the list.
COMMIT_MESSAGE_MENTIONS = {
    "add-then-commit-m": f'git add -A && git commit -m "guard denies grep over {H}/r.jsonl"',
    "add-then-commit-qm": f'git add -A && git commit -qm "x {H}"',
    "add-then-commit-message-equals": f'git add -A && git commit --message="x {H}"',
    "add-then-commit-glued-m": f'git add -A && git commit -m"x {H}"',
}


@pytest.mark.xfail(strict=True, reason="the guard edit waits for Ola's approval")
@pytest.mark.parametrize("command", COMMIT_MESSAGE_MENTIONS.values(),
                         ids=COMMIT_MESSAGE_MENTIONS.keys())  # fmt: skip
def test_a_mention_only_inside_a_commit_message_passes(guard: ModuleType, command: str) -> None:
    assert guard.judge_bash(command) is None


def test_the_folder_lists_hold_each_line_once() -> None:
    passes, denials = record_folder_passes(), record_folder_denials()
    assert not set(passes.values()) & set(denials.values())
    assert len(set(denials.values())) == len(denials)


@pytest.mark.parametrize("mode", ["off", "on"])
def test_a_mention_in_a_heredoc_passes_through_the_hook(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    command = record_folder_passes()["quoted-heredoc-into-session"]
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command))) is None
    assert queue_lines(repo) == []


@pytest.mark.parametrize("mode", ["off", "on"])
def test_a_case_changed_write_is_denied_through_the_hook(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    command = record_folder_denials()["upper-case-target"]
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, bash_event(repo, command)))
    assert found == ("deny", ALWAYS_DENIED)
    assert queue_lines(repo) == []


def test_shell_scan_marks_a_leading_assignment() -> None:
    shell_scan = load_tool("shell_scan")
    [assigned] = shell_scan.parse("RIPGREP_CONFIG_PATH=/tmp/c rg x y")
    [bare] = shell_scan.parse("rg x y")
    assert (assigned.assigns, bare.assigns) == (True, False)


@pytest.mark.parametrize(
    ("command", "expected"),
    [("cat > n.md <<EOF\n$(touch x)\nEOF", True),
     ("cat > n.md <<EOF\n`touch x`\nEOF", True),
     ("cat > n.md <<'EOF'\n$(touch x)\nEOF", False),
     ("cat > n.md <<EOF\nplain $HOME\nEOF", False)],
    ids=["unquoted-dollar-paren", "unquoted-backtick", "quoted", "no-command"],
)  # fmt: skip
def test_shell_scan_marks_a_heredoc_body_the_shell_runs(command: str, expected: bool) -> None:
    shell_scan = load_tool("shell_scan")
    owner = [s for s in shell_scan.parse(command) if s.argv[:1] == ["cat"]]
    assert [s.runs_body for s in owner] == [expected]


# ---------------------------------------------------------------- guard follow-ups (h37)
#
# docs/increments/h37-guard-follow-ups.md §4: every line §1 found passing is
# refused (asked for the zsh redirects and the letter-case rows, denied for the
# record folder), except the rows §2 leaves open; the lines §4 keeps, and every
# line §1 found refused, are judged as they are today. The folder name is built
# from parts, as above.

HN = "har" + "ness"
GLOB_HN = HN[:-1] + "?"

#: §1 (1): zsh's clobbering and combined redirects, each writing the file named.
ZSH_OPERATORS = (">!", ">>!", ">>|", "&>!", ">&!", ">&|", ">>&", "&>>!", "&>>|", ">>&!",
                 ">>&|", "2>!", "2>>|")  # fmt: skip


def zsh_redirect_asks() -> dict[str, tuple[str, str]]:
    """§1 (1), passing today: id -> (line, the governed target it must be asked for)."""
    rows = {f"echo{op}": (f"echo x {op} CLAUDE.md", "CLAUDE.md") for op in ZSH_OPERATORS}
    return rows | {
        "no-space-append-clobber": ("echo x >>|CLAUDE.md", "CLAUDE.md"),
        "clobber-alone": (">! CLAUDE.md", "CLAUDE.md"),
        "exec-fd-clobber": ("exec {fd}>! CLAUDE.md", "CLAUDE.md"),
        "clobber-settings": ("echo x >! .claude/settings.json", ".claude/settings.json"),
        "clobber-away": ("echo x >! tools/away.py", "tools/away.py"),
    }


#: §1 (2), case: the macOS disk ignores letter case, so each writes a governed file.
CASE_ASKS = {
    "lower-claude": ("echo x > claude.md", "claude.md"),
    "agents-capital": ("echo x > .Claude/agents/x.md", ".Claude/agents/x.md"),
    "away-capitals": ("echo x > Tools/Away.py", "Tools/Away.py"),
    "principles-lower": ("echo x > docs/principles.md", "docs/principles.md"),
}

#: §1 (2): the folder written with `.` and doubled `/` segments.
FOLDER_SPELLINGS = (f".git//{HN}/m", f".git/./{HN}/m", f"./.git//{HN}/m", f".git///{HN}/m",
                    f"/r/.git//{HN}/m")  # fmt: skip


def folder_denials() -> dict[str, str]:
    """§1 (2) and (3), passing today; each a write of the record folder, to be denied."""
    lines = {f"{form.split()[0]}-{path}": form.format(path)
             for path in FOLDER_SPELLINGS
             for form in ("echo x > {}", "touch {}", "rm {}")}  # fmt: skip
    return lines | {
        "cd-git-dir-folder-touch": f"cd $(git rev-parse --git-dir)/{HN} && touch x",
        "rsync-into-git-dir-folder": f"rsync -a s/ $(git rev-parse --git-dir)/{HN}/",
        "cd-git-dir-touch-folder": f"cd $(git rev-parse --git-dir) && touch {HN}/x",
        "touch-glob-question": f"touch .git/{GLOB_HN}/x",
        "redirect-glob-bracket": f"echo x > .git/{HN[:-1]}[s]/x",
        "touch-glob-git": f"touch .gi?/{HN}/x",
        "redirect-dot-star": f"echo x > .*/{HN}/m",
        # §2 rule 4a: a link to the git folder, made and written through on one line.
        "link-git-then-write": f"ln -s .git g && echo x > g/{HN}/m",
        "link-common-dir": "ln -s $(git rev-parse --git-common-dir) g",
        "link-absolute-git-dir": "ln -s $(git rev-parse --absolute-git-dir) g",
        # §3's rule 4a row: the main checkout's git folder linked into a scratchpad.
        "link-assigned-git": "ln -sfn $R/.git $SP/lk/.git",
        # §3's rows 2 and 3: a non-reader on a folder path held in an expansion.
        "sed-on-common-dir-folder": f"G=$(git rev-parse --git-common-dir); sed -n 1p $G/{HN}/x",
    }


#: §4: lines kept passing (beside h33's list above, which reruns as it is).
KEPT = {
    "cd-job-tmp-pytest": "cd $CLAUDE_JOB_DIR/tmp/x && pytest",
    "cd-toplevel-ls": "cd $(git rev-parse --show-toplevel) && ls",
    "echo-git-dir-folder": f"echo $(git rev-parse --git-dir)/{HN}",
    "make-stderr-grep": "make 2>&1 | grep x",
    "echo-to-stderr": "echo x >&2",
    "exec-close-fd": "exec 3>&-",
    "git-grep-conflict-markers": "git grep -nE '^(<<<<<<<|>>>>>>>)'",
    "link-build": "ln -s ../build b",
    "link-data": "ln -s $DATA data",
    "rm-glob-in-expansion": "rm $D/*",
    "find-not-git": "find . -not -path '*/.git/*' -name x",
    "known-scratch-cd": f"S=/tmp/s; cd $S/{HN} && echo x > y",
    # §3: reads of the folder through an expansion stay open to the readers.
    "cat-common-dir-folder": f"G=$(git rev-parse --git-common-dir); cat $G/{HN}/queue.jsonl",
    # §2 rule 5: the standard-library-module rule stays exact in case.
    "stdlib-name-other-case": "echo x > tools/Ast.py",
}

#: §2 *Stays open*: pinned as passing; only the operating system's sandbox closes them.
STAYS_OPEN = {
    "script-shell": "bash s.sh",
    "script-python": "python3 tools/x.py",
    "cd-expansion-unknown-writer": f"cd $X && rsync s/ {HN}/",
    # §2's example `H=$G/harness; touch $H/x` is denied today (the scanner puts
    # the assignment into the target), so it is pinned under STAYS_REFUSED; a
    # variable set to a command substitution is the form that stays open.
    "variable-set-to-substitution": f"H=$(echo $G/{HN}); touch $H/x",
}

#: §1's refused rows: each stays refused (asked or denied).
STAYS_REFUSED = {
    "redirect": "echo x > CLAUDE.md",
    "redirect-pipe": "echo x >| CLAUDE.md",
    "both-pipe": "echo x &>| CLAUDE.md",
    "both-append": "echo x &>> CLAUDE.md",
    "clobber-no-space": "echo x >!CLAUDE.md",
    "exec-fd": "exec {fd}> CLAUDE.md",
    "hooks-up-to-folder": f"echo x > .git/hooks/../{HN}/m",
    "folder-double-slash-inside": f"echo x > .git/{HN}//m",
    "folder-up-and-back": f"echo x > .git/{HN}/../{HN}/m",
    "backslash-dot": f"echo x > \\.git/{HN}/m",
    "quoted-inside-git": f'echo x > .g"it"/{HN}/m',
    "cd-double-slash": f"cd .git//{HN} && touch x",
    "cd-git-dot": f"cd .git/. && touch {HN}/x",
    "cp-common-dir": f"cp x $(git rev-parse --git-common-dir)/{HN}/m",
    "redirect-common-dir": f"echo x > $(git rev-parse --git-common-dir)/{HN}/m",
    "variable-set-to-expansion": f"H=$G/{HN}; touch $H/x",
}


@pytest.mark.parametrize(("command", "target"), zsh_redirect_asks().values(),
                         ids=zsh_redirect_asks().keys())  # fmt: skip
def test_h37_a_zsh_redirect_to_a_governed_file_asks(
    guard: ModuleType, command: str, target: str
) -> None:
    found = guard.judge_bash(command)
    assert found is not None, f"{command!r} passed"
    assert found[0] == "ask" and target in found[1], found


def test_h37_an_unreadable_line_with_a_zsh_redirect_asks_by_its_text(guard: ModuleType) -> None:
    cannot_read = "the guard cannot read this command's targets"
    assert guard.judge_bash("echo 'u >! CLAUDE.md") == ("ask", ["CLAUDE.md"], cannot_read)


@pytest.mark.parametrize(("command", "target"), CASE_ASKS.values(), ids=CASE_ASKS.keys())
def test_h37_a_governed_name_in_another_letter_case_asks(
    guard: ModuleType, command: str, target: str
) -> None:
    found = guard.judge_bash(command)
    assert found is not None, f"{command!r} passed"
    assert found[0] == "ask" and target in found[1], found


@pytest.mark.parametrize("command", folder_denials().values(), ids=folder_denials().keys())
def test_h37_a_write_of_the_folder_by_any_spelling_is_denied(
    guard: ModuleType, command: str
) -> None:
    assert guard.judge_bash(command) == DENIED


@pytest.mark.parametrize("command", KEPT.values(), ids=KEPT.keys())
def test_h37_a_kept_line_passes(guard: ModuleType, command: str) -> None:
    assert guard.judge_bash(command) is None


@pytest.mark.parametrize("command", STAYS_OPEN.values(), ids=STAYS_OPEN.keys())
def test_h37_a_line_the_design_leaves_open_passes(guard: ModuleType, command: str) -> None:
    # §2 *Stays open*: what a program writes is fixed when it runs; the full
    # close is Claude Code's sandbox, a question for Ola. Pinned so a change is seen.
    assert guard.judge_bash(command) is None


@pytest.mark.parametrize("command", STAYS_REFUSED.values(), ids=STAYS_REFUSED.keys())
def test_h37_a_line_refused_today_stays_refused(guard: ModuleType, command: str) -> None:
    assert guard.judge_bash(command) is not None


def test_h37_the_lists_hold_each_line_once() -> None:
    refused = [c for c, _ in (*zsh_redirect_asks().values(), *CASE_ASKS.values())]
    refused += [*folder_denials().values(), *STAYS_REFUSED.values()]
    passed = [*KEPT.values(), *STAYS_OPEN.values(), *record_folder_passes().values()]
    assert len(set(refused)) == len(refused) and not set(refused) & set(passed)


@pytest.mark.parametrize("mode", ["off", "on"])
@pytest.mark.parametrize("spelling", [f".git//{HN}/m", f".git/./{HN}/m"])
def test_h37_a_file_tool_write_to_the_folder_by_any_spelling_is_denied(
    repo: Path, mode: str, spelling: str
) -> None:
    set_mode(repo, mode)
    event = file_event(repo, "Write", f"{repo}/{spelling}")
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found == ("deny", ALWAYS_DENIED)


def test_h37_a_file_tool_write_through_hooks_up_to_the_folder_stays_refused(repo: Path) -> None:
    event = file_event(repo, "Write", f"{repo}/.git/hooks/../{HN}/m")
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is not None


def test_h37_a_file_tool_write_to_a_governed_name_in_another_case_asks(repo: Path) -> None:
    event = file_event(repo, "Write", f"{repo}/claude.md")
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found is not None and found[0] == "ask", found


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h37_a_file_tool_write_to_a_note_named_like_claude_passes(repo: Path, mode: str) -> None:
    set_mode(repo, mode)
    event = file_event(repo, "Write", f"{repo}/notes/claude-notes.md")
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None


@pytest.fixture
def git_link(repo: Path) -> Path:
    """`g` in the repository, a link an earlier line made to its git folder (§2 rule 4b)."""
    (repo / "g").symlink_to(repo / ".git")
    return repo / "g"


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h37_a_shell_write_through_a_link_made_earlier_is_denied(
    repo: Path, git_link: Path, mode: str
) -> None:
    set_mode(repo, mode)
    event = bash_event(repo, f"echo x > g/{HN}/m")
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)


def test_h37_judge_bash_follows_a_link_from_the_given_directory(
    guard: ModuleType, repo: Path, git_link: Path
) -> None:
    assert guard.judge_bash(f"echo x > g/{HN}/m", cwd=str(repo)) == DENIED


def test_h37_a_file_tool_write_through_a_link_made_earlier_is_denied(
    repo: Path, git_link: Path
) -> None:
    event = file_event(repo, "Write", f"{git_link}/{HN}/m")
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) == ("deny", ALWAYS_DENIED)


@pytest.mark.parametrize("tool", ["Bash", "Write"])
def test_h37_a_write_through_a_link_to_a_governed_file_asks(repo: Path, tool: str) -> None:
    (repo / "n").symlink_to(repo / "CLAUDE.md")
    event = (bash_event(repo, "echo x > n") if tool == "Bash"
             else file_event(repo, "Write", str(repo / "n")))  # fmt: skip
    found = pretool_decision(run_script(repo, GUARD_GOVERNANCE, event))
    assert found is not None and found[0] == "ask", found


@pytest.mark.parametrize("mode", ["off", "on"])
def test_h37_a_link_to_a_git_folder_under_a_scratchpad_passes(
    repo: Path, pad: Path, mode: str
) -> None:
    set_mode(repo, mode)
    event = bash_event(repo, f"ln -s {pad}/r/.git g")
    assert pretool_decision(run_script(repo, GUARD_GOVERNANCE, event)) is None
    assert queue_lines(repo) == []
