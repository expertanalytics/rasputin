"""The guards judge the files a shell line writes, not the words in it (h4).

The spec is `docs/increments/h4-guard-fixes.md` §3 and §4, with Ola's two
rulings of §6: a line the guards cannot read is judged by today's text rule,
and a refusal names the targets it judged.

The commands are the probe's `CASES`, loaded from
`docs/increments/h4-probes/h4_probe.py` so they stay verbatim; §4's further
commands are added under labels of their own. Every guard case runs by day
(flag off) and at night (flag on), in the U1 fixture's temporary repository.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

from harness_fixtures import (
    GUARD_GOVERNANCE,
    GUARD_PUSH,
    REAL,
    Tool,
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


def _load_probe() -> ModuleType:
    path = REAL / "docs" / "increments" / "h4-probes" / "h4_probe.py"
    spec = importlib.util.spec_from_file_location("h4_probe", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PROBE = _load_probe()
AWAY: str = PROBE.AWAY

#: label -> (tool, payload, agent_type). The probe's cases, then §4's own.
CASES: dict[str, tuple[str, str, str | None]] = {
    label: (tool, payload, agent) for label, tool, payload, agent in PROBE.CASES
} | {
    "gh-pr-body": ("Bash", "gh pr create --body 'x git push'", None),
    "rm-governed": ("Bash", "rm CLAUDE.md", None),
    "expansion-tail": ("Bash", "echo x > $D/CLAUDE.md", None),
    "bash-c-sed": ("Bash", 'bash -c "sed -i s/a/b/ CLAUDE.md"', None),
    "env-away": ("Bash", f"env python3 tools/{AWAY} 8h", None),
    "subst-away": ("Bash", f'echo "$(python3 tools/{AWAY} --back)"', None),
    "unbalanced": ("Bash", "echo 'unbalanced > CLAUDE.md", None),
    "for-loop": ("Bash", "for f in CLAUDE.md; do sed -i s/a/b/ $f; done", None),
    "rm-expansion": ("Bash", "rm $f", None),
    "git-apply": ("Bash", "git apply x.patch", None),
    # The adversarial comparison's three weakenings (W1-W3): each was judged
    # by text before h4 and must not slip past the target-based reading.
    "eval-push": ("Bash", 'eval "git push"', None),
    "git-mv-governed": ("Bash", "git mv notes.txt CLAUDE.md", None),
    "assigned-harness": ("Bash", "H=.git/harness; echo x > $H/unattended.json", None),
    "assigned-governed": ("Bash", "G=CLAUDE; echo x > $G.md", None),
    # The re-run comparison's further routes to the harness state, each denied
    # by the pre-h4 text rule: a variable set by export, read, a substitution or
    # declare, and a relative write after cd into the state dir.
    "export-harness": ("Bash", "export H=.git/harness && echo x > $H/unattended.json", None),
    "read-harness": ("Bash", "read H <<< .git/harness; echo x > $H/unattended.json", None),
    "subst-harness": ("Bash", "H=$(echo .git/harness); echo x > $H/unattended.json", None),
    "declare-harness": ("Bash", "declare H=.git/harness; echo x > $H/unattended.json", None),
    "cd-harness": ("Bash", "cd .git/harness && echo x > unattended.json", None),
    # A here-string commit message naming git push: a commit, not a push.
    "commit-herestring": ("Bash", 'git commit -F - <<< "git push"', None),
}

#: §2's five false positives, and a here-string commit message naming git push:
#: neither guard speaks, by day or at night.
FALSE_POSITIVES = (
    "fp1-write",
    "fp1-bash",
    "fp2-heredoc",
    "fp2-body-file",
    "fp2-commit-msg",
    "fp3-sed",
    "fp3-heredoc",
    "fp3-commit-msg",
    "fp4-git-show",
    "fp4-diff-to-tmp",
    "fp5-grep",
    "fp5-echo-memory",
    "fp5-printf-task",
    "commit-herestring",
)

#: Publishing acts and the one reason each gives, as today. `gh-pr-body` names
#: `git push` only in its body, which must not add a second reason.
PUSH_ASKS = {
    "tp-push": "git push writes to the remote",
    "tp-rebase": "this rewrites history, which is destructive once anything is published",
    "gh-pr-body": "gh pr changes a pull request",
    "eval-push": "git push writes to the remote",  # W1: eval runs its words
}

#: Governed writes and the governed targets the ask and the refusal must name.
GOVERNED_WRITES = {
    "tp-claude-md-sed": ("CLAUDE.md",),
    "tp-readme-heredoc": ("docs/increments/README.md",),
    "tp-settings-py": (".claude/settings.json",),
    "tp-check-cp": ("tools/check_citations.py",),
    "tp-edit-claude": ("CLAUDE.md",),
    "rm-governed": ("CLAUDE.md",),
    "expansion-tail": ("CLAUDE.md",),
    "bash-c-sed": ("CLAUDE.md",),
    "git-mv-governed": ("CLAUDE.md",),  # W2: git mv writes its destination
}

#: W3: a governed target built from a variable assigned on the same line. What
#: the refusal says is not pinned; that it asks by day and is refused at night is.
ASSIGNED_GOVERNED = "assigned-governed"

#: Runs of the away script and writes of the harness state: denied in both modes.
#: `assigned-harness` (W3) reaches the state through a variable set on the line.
ALWAYS_DENIED_CASES = (
    "tp-run-away",
    "tp-harness-write",
    "env-away",
    "subst-away",
    "assigned-harness",
    "export-harness",
    "read-harness",
    "subst-harness",
    "declare-harness",
    "cd-harness",
)
ALWAYS_DENIED = (
    "Only Ola enters or leaves unattended mode, and only hooks and away.py write the "
    "harness state. Nothing is queued: this act is not an agent's to wait for."
)

#: Cannot judge (Ola's ruling 1): today's text rule decides, so a governed name
#: plus a write construct asks by day and is refused and queued at night ...
CANNOT_JUDGE_GOVERNED = ("unbalanced", "for-loop")
#: ... and without a governed name the line passes.
CANNOT_JUDGE_PLAIN = ("rm-expansion", "git-apply")
CANNOT_READ = "the guard cannot read this command's targets"

MODES = ("off", "on")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


def decisions(repo: Path, label: str) -> dict[str, tuple[str, str]]:
    """{guard: (decision, reason)} for each guard that speaks on the case, as the probe runs it."""
    tool, payload, agent = CASES[label]
    extra = {"agent_id": "p1", "agent_type": agent} if agent else {}
    if tool == "Bash":
        event = bash_event(repo, payload, **extra)
        hooks = {"governance": GUARD_GOVERNANCE, "push": GUARD_PUSH}
    else:
        event = file_event(repo, tool, payload, **extra)
        hooks = {"governance": GUARD_GOVERNANCE}
    found = {name: pretool_decision(run_script(repo, hook, event)) for name, hook in hooks.items()}
    return {name: verdict for name, verdict in found.items() if verdict is not None}


def act_of(label: str) -> str:
    tool, payload, _ = CASES[label]
    return payload if tool == "Bash" else f"{tool} {payload}"


# ---------------------------------------------------------------- the false positives


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("label", FALSE_POSITIVES)
def test_a_false_positive_passes_both_guards(repo: Path, mode: str, label: str) -> None:
    set_mode(repo, mode)
    assert decisions(repo, label) == {}
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- publishing acts


@pytest.mark.parametrize("label", sorted(PUSH_ASKS))
def test_a_publishing_act_asks_by_day_with_one_reason(repo: Path, label: str) -> None:
    found = decisions(repo, label)
    assert set(found) == {"push"}
    kind, reason = found["push"]
    assert kind == "ask"
    assert f"This reaches beyond the working tree: {PUSH_ASKS[label]}.\n" in reason
    assert queue_lines(repo) == []


@pytest.mark.parametrize("label", sorted(PUSH_ASKS))
def test_a_publishing_act_is_refused_and_queued_at_night(repo: Path, label: str) -> None:
    set_mode(repo, "on")
    found = decisions(repo, label)
    assert set(found) == {"push"}
    assert found["push"][0] == "deny"
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"], line["why"]) == (
        "guard_push",
        act_of(label),
        PUSH_ASKS[label],
    )


# ---------------------------------------------------------------- governed writes


@pytest.mark.parametrize("label", sorted(GOVERNED_WRITES))
def test_a_governed_write_asks_by_day_naming_its_targets(repo: Path, label: str) -> None:
    found = decisions(repo, label)
    assert set(found) == {"governance"}
    kind, reason = found["governance"]
    assert kind == "ask"
    for target in GOVERNED_WRITES[label]:
        assert target in reason, f"the ask does not name {target}: {reason!r}"
    assert queue_lines(repo) == []


@pytest.mark.parametrize("label", sorted(GOVERNED_WRITES))
def test_a_governed_write_is_refused_at_night_naming_its_targets(repo: Path, label: str) -> None:
    set_mode(repo, "on")
    found = decisions(repo, label)
    assert set(found) == {"governance"}
    kind, reason = found["governance"]
    assert kind == "deny"
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_governance", act_of(label))
    why = line["why"]
    assert why.startswith("it writes "), why
    for target in GOVERNED_WRITES[label]:
        assert target in why, f"the refusal does not name {target}: {why!r}"
    assert why in reason


def test_a_governed_target_built_on_the_line_asks_by_day(repo: Path) -> None:
    """W3: `G=CLAUDE; echo x > $G.md` writes CLAUDE.md, whatever the static tail says."""
    found = decisions(repo, ASSIGNED_GOVERNED)
    assert set(found) == {"governance"}
    assert found["governance"][0] == "ask"
    assert queue_lines(repo) == []


def test_a_governed_target_built_on_the_line_is_refused_and_queued_at_night(repo: Path) -> None:
    set_mode(repo, "on")
    found = decisions(repo, ASSIGNED_GOVERNED)
    assert set(found) == {"governance"}
    assert found["governance"][0] == "deny"
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"]) == ("guard_governance", act_of(ASSIGNED_GOVERNED))


# ---------------------------------------------------------------- always denied


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("label", ALWAYS_DENIED_CASES)
def test_a_run_of_away_or_a_state_write_is_always_denied(repo: Path, mode: str, label: str) -> None:
    set_mode(repo, mode)
    assert decisions(repo, label) == {"governance": ("deny", ALWAYS_DENIED)}
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- cannot judge


@pytest.mark.parametrize("label", CANNOT_JUDGE_GOVERNED)
def test_an_unreadable_governed_line_asks_by_day(repo: Path, label: str) -> None:
    found = decisions(repo, label)
    assert set(found) == {"governance"}
    kind, reason = found["governance"]
    assert kind == "ask"
    assert "CLAUDE.md" in reason
    assert queue_lines(repo) == []


@pytest.mark.parametrize("label", CANNOT_JUDGE_GOVERNED)
def test_an_unreadable_governed_line_is_refused_at_night_saying_so(repo: Path, label: str) -> None:
    set_mode(repo, "on")
    found = decisions(repo, label)
    assert set(found) == {"governance"}
    kind, reason = found["governance"]
    assert kind == "deny"
    assert CANNOT_READ in reason
    [line] = queue_lines(repo)
    assert (line["hook"], line["act"], line["why"]) == (
        "guard_governance",
        act_of(label),
        CANNOT_READ,
    )


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("label", CANNOT_JUDGE_PLAIN)
def test_an_unreadable_line_naming_no_governed_file_passes(
    repo: Path, mode: str, label: str
) -> None:
    set_mode(repo, mode)
    assert decisions(repo, label) == {}
    assert queue_lines(repo) == []


# ---------------------------------------------------------------- shell_scan.parse

shell_scan = Tool("shell_scan")

#: label -> every simple command's (argv, writes), in order. Interpreter
#: commands whose program has string literals are left to PROGRAMS: whether a
#: literal is a candidate target is the guard's reading, not pinned here.
PARSED: dict[str, list[tuple[list[str], list[str]]]] = {
    "fp1-bash": [
        (["printf", "ask: x\\n"], [".claude/current-task/reviewer-231000.md"]),
        (["cat", ".claude/REQUIRED-READING.md"], []),
    ],
    "fp2-heredoc": [(["cat"], [".claude/current-task/session.md"])],
    "fp2-body-file": [(["cat"], ["$CLAUDE_JOB_DIR/tmp/body.md"])],
    "fp2-commit-msg": [
        (["git", "commit", "-q", "-F", "$CLAUDE_JOB_DIR/tmp/msg.txt"], []),
        (["echo", "next: git push once Ola says yes"], []),
    ],
    "fp3-sed": [
        (
            ["sed", "-i", "", "s/a/b/", "docs/increments/23-basin-scale.md"],
            ["docs/increments/23-basin-scale.md"],
        ),
        (["python3", "tools/check_citations.py"], []),
        (["tail", "-3"], []),
    ],
    "fp3-commit-msg": [
        (["printf", "docs: cite CLAUDE.md section 2\\n"], ["$CLAUDE_JOB_DIR/tmp/m.txt"]),
        (["git", "commit", "-q", "-F", "$CLAUDE_JOB_DIR/tmp/m.txt"], []),
    ],
    "fp4-git-show": [
        (["git", "show", "agents-web-tools:.claude/settings.json"], []),
        (["python3", "-c", "import json,sys; print(json.load(sys.stdin).keys())"], []),
    ],
    "fp4-diff-to-tmp": [
        (["python3", "tools/check_citations.py"], []),
        (["git", "diff", "-U0", "HEAD"], ["$CLAUDE_JOB_DIR/tmp/diff.txt"]),
        (["python3", "-"], []),
    ],
    "fp5-grep": [
        (["grep", "-n", "-B3", "-A12", "def _open_tty\\|/dev/tty", f"tools/{AWAY}"], []),
        (["head", "-60"], []),
        (["which", "-a", "python3"], []),
    ],
    "fp5-echo-memory": [
        (["echo", f"- each night is a trial of {AWAY} + guards; report idle time"], ["MEMORY.md"]),
    ],
    "fp5-printf-task": [
        (
            ["printf", f"Ask: spec amendment for {AWAY} /dev/tty defect; commit, no push.\\n"],
            [".claude/current-task/architect-224000.md"],
        ),
    ],
    "tp-push": [(["git", "push", "-u", "origin", "x"], []), (["tail", "-2"], [])],
    "tp-rebase": [(["cd", "w"], []), (["git", "rebase", "master"], []), (["tail", "-3"], [])],
    "tp-claude-md-sed": [(["sed", "-i", "", "s/a/b/", "CLAUDE.md"], ["CLAUDE.md"])],
    "tp-readme-heredoc": [(["cat"], ["docs/increments/README.md"])],
    "tp-check-cp": [
        (["cp", "/tmp/x.py", "tools/check_citations.py"], ["tools/check_citations.py"])
    ],
    "tp-run-away": [(["python3", f"tools/{AWAY}", "--back"], [])],
    "tp-harness-write": [(["echo", "{}"], [".git/harness/unattended.json"])],
    "gh-pr-body": [(["gh", "pr", "create", "--body", "x git push"], [])],
    "rm-governed": [(["rm", "CLAUDE.md"], ["CLAUDE.md"])],
    "expansion-tail": [(["echo", "x"], ["$D/CLAUDE.md"])],
    "env-away": [(["python3", f"tools/{AWAY}", "8h"], [])],
    "rm-expansion": [(["rm", "$f"], ["$f"])],
}

#: label -> (argv of the interpreter command, its program: the -c text or heredoc body).
PROGRAMS: dict[str, tuple[list[str], str]] = {
    "fp3-heredoc": (
        ["python3", "-"],
        "p='docs/increments/23-basin-scale.md'; s=open(p).read()\n"
        "s=s.replace('see CLAUDE.md', 'see docs/increments/README.md')\n"
        "open(p,'w').write(s)",
    ),
    "fp4-git-show": (
        ["python3", "-c", "import json,sys; print(json.load(sys.stdin).keys())"],
        "import json,sys; print(json.load(sys.stdin).keys())",
    ),
    "fp4-diff-to-tmp": (["python3", "-"], "print(1)"),
    "tp-settings-py": (
        ["python3", "-c", "open('.claude/settings.json','w').write('{}')"],
        "open('.claude/settings.json','w').write('{}')",
    ),
}

#: label -> one (argv, writes) that must be among the simple commands, where
#: §3 fixes the inner command but not how the enclosing one is reported.
CONTAINS: dict[str, tuple[list[str], list[str] | None]] = {
    "bash-c-sed": (["sed", "-i", "s/a/b/", "CLAUDE.md"], ["CLAUDE.md"]),
    "subst-away": (["python3", f"tools/{AWAY}", "--back"], []),
    "for-loop": (["sed", "-i", "s/a/b/", "$f"], ["$f"]),
    "git-apply": (["git", "apply", "x.patch"], None),  # its targets are unknown
    "eval-push": (["git", "push"], []),  # W1: eval's words are a command
}

#: Unparseable (§3): unbalanced quotes, an unclosed `$(`, an unterminated heredoc.
UNPARSEABLE = (
    CASES["unbalanced"][1],
    "echo $(ls tools",
    "cat > notes.txt <<'EOF'\nno terminator",
)

#: Lines with no heredoc: a program, when present, is the whole -c text.
_NO_HEREDOC = {"fp4-git-show", "tp-settings-py"}


@pytest.mark.parametrize("label", sorted(PARSED))
def test_parse_gives_each_simple_command_with_its_writes(label: str) -> None:
    simples = shell_scan.parse(CASES[label][1])
    assert simples is not None, "parse gave up on a readable line"
    assert [(s.argv, s.writes) for s in simples] == PARSED[label]


@pytest.mark.parametrize("label", sorted(PARSED))
def test_parse_reads_no_heredoc_or_quoted_text_as_a_program_of_a_non_interpreter(
    label: str,
) -> None:
    simples = shell_scan.parse(CASES[label][1])
    assert simples is not None
    for simple in simples:
        if simple.argv[0] != "python3":
            assert simple.program is None, f"{simple.argv} has a program"


@pytest.mark.parametrize("label", sorted(PROGRAMS))
def test_parse_gives_an_interpreter_its_program(label: str) -> None:
    simples = shell_scan.parse(CASES[label][1])
    assert simples is not None
    argv, program = PROGRAMS[label]
    [simple] = [s for s in simples if s.argv == argv]
    assert simple.program is not None
    if label in _NO_HEREDOC:
        assert simple.program == program
    else:  # a heredoc body: with or without its last newline
        assert simple.program.rstrip("\n") == program


@pytest.mark.parametrize("label", sorted(CONTAINS))
def test_parse_reaches_into_wrappers_loops_and_substitutions(label: str) -> None:
    simples = shell_scan.parse(CASES[label][1])
    assert simples is not None, "parse gave up on a readable line"
    argv, writes = CONTAINS[label]
    matching = [s for s in simples if s.argv == argv]
    assert matching, f"no simple command {argv} in {[s.argv for s in simples]}"
    if writes is not None:
        assert [s.writes for s in matching] == [writes]


def test_parse_counts_the_destination_of_git_mv_as_written() -> None:
    """W2: `git mv` writes its destination, like `mv`."""
    simples = shell_scan.parse(CASES["git-mv-governed"][1])
    assert simples is not None, "parse gave up on a readable line"
    assert "CLAUDE.md" in simples[0].writes, simples[0].writes


@pytest.mark.parametrize("command", UNPARSEABLE)
def test_parse_gives_none_for_an_unparseable_line(command: str) -> None:
    assert shell_scan.parse(command) is None


def test_a_move_into_a_directory_named_without_a_slash_names_both_readings() -> None:
    """h16 §2 G8: the last operand as itself and as a directory; `mv` keeps its sources."""
    found, unnamed = shell_scan.writer_targets("mv", ["a", "b", "tools"])
    assert {"a", "b", "tools", "tools/a", "tools/b"} <= set(found), found
    assert unnamed is False
