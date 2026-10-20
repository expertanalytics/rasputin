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

import ast
import importlib.util
import shlex
from pathlib import Path
from types import ModuleType
from typing import Any

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


# ---------------------------------------------------------------- h24 item 5b
#
# docs/increments/h24-harness-fixes.md §3 item 5b: a Python string literal
# that only flows into a read (`json.load(open(LIT))`, `open(LIT).read()`,
# `Path(LIT).read_text()`, ...) is not a write candidate, provided the whole
# program imports only json and pathlib, rebinds none of the names the read
# relies on, and patches or reflects on nothing. Otherwise every literal is a
# candidate, as before. Each program runs through `parse` as Python's `-c`
# text; the first shape also as a heredoc.

#: The harness folder's `unattended.json`, as an absolute path.
H = "/x/.git/harness/unattended.json"
#: A function the program defines itself: no import, no rebinding of a read name.
W = "def W(*args, **kwargs):\n    return 0\n"


def candidates_of(program: str, heredoc: bool = False) -> list[str]:
    """`candidates` of a Python program that parses, so the AST rules judge it, not the regex."""
    ast.parse(program)
    line = f"python3 - <<'X'\n{program}\nX" if heredoc else f"python3 -c {shlex.quote(program)}"
    simples = shell_scan.parse(line)
    assert simples is not None and len(simples) == 1, simples
    assert simples[0].program is not None and simples[0].program.strip() == program.strip()
    found: list[str] = shell_scan.candidates(simples[0])
    return found


READS: dict[str, tuple[str, bool]] = {
    "json-load-by-c": (f"import json; print(json.load(open({H!r})))", False),
    "json-load-by-heredoc": (f"import json\nprint(json.load(open({H!r})))", True),
    "path-read-text": (f"from pathlib import Path; print(Path({H!r}).read_text())", False),
    "open-r-encoding-readlines": (
        f'print(len(open({H!r}, "r", encoding="utf-8").readlines()))', False),
    "open-rb-read": (f'print(open({H!r}, "rb").read())', False),
    "pathlib-path-read-bytes": (f"import pathlib; print(pathlib.Path({H!r}).read_bytes())", False),
}  # fmt: skip


@pytest.mark.parametrize(("program", "heredoc"), READS.values(), ids=READS.keys())
def test_h24_5b_a_literal_only_read_is_not_a_candidate(program: str, heredoc: bool) -> None:
    assert H not in candidates_of(program, heredoc)


STILL_CANDIDATES: dict[str, str] = {
    # The design review's writes, each its own case.
    "reused-name": f'n = open({H!r}).name; open(n, "w").write("x")',
    "opener": f"{W}print(open({H!r}, opener=W).read())",
    "builtins-open-patched": f"{W}import builtins; builtins.open = W; print(open({H!r}).read())",
    "path-read-text-patched": (
        f"{W}from pathlib import Path; Path.read_text = W; print(Path({H!r}).read_text())"),
    "setattr-before-the-read": (
        f'{W}from pathlib import Path; setattr(Path, "read_text", W); '
        f"print(Path({H!r}).read_text())"),
    "import-as-json": f"import shutil as json; json.load(open({H!r}))",
    "path-from-elsewhere": f"from x import Path; Path({H!r}).read_text()",
    "star-import": f"from os import *; print(open({H!r}).read())",
    "json-rebound": f"{W}import json; json = W; json.load(open({H!r}))",
    # The reflection route, and an import outside the allowlist.
    "code-object": f"print(open({H!r}).read()); (lambda: 0).__code__",
    "import-sys": f"import sys; print(open({H!r}).read())",
    # Unchanged: writes, unread modes, named file objects, and built paths.
    "open-w": f'open({H!r}, "w")',
    "open-a": f'open({H!r}, "a")',
    "open-mode-variable": f'mode = "r"; open({H!r}, mode)',
    "mode-keyword": f'open({H!r}, mode="r").read()',
    "path-write-text": f'from pathlib import Path; Path({H!r}).write_text("x")',
    "path-open-w": f'from pathlib import Path; Path({H!r}).open("w")',
    "named-path": f"from pathlib import Path; p = Path({H!r}); p.read_text()",
    "named-file": f"f = open({H!r}); f.read()",
    "for-loop-over-open": f"for l in open({H!r}): print(l)",
    "with-open-as": f"with open({H!r}) as f: print(f.read())",
    "shutil-copy": f'import shutil; shutil.copy("a.json", {H!r})',
    "open-redefined": f"def open(p):\n    return 0\nopen({H!r}).read()",
}  # fmt: skip


@pytest.mark.parametrize("program", STILL_CANDIDATES.values(), ids=STILL_CANDIDATES.keys())
def test_h24_5b_a_literal_outside_the_exempt_reads_stays_a_candidate(program: str) -> None:
    assert H in candidates_of(program)


def test_h24_5b_the_harness_folder_joined_with_slash_stays_a_candidate() -> None:
    folder = H.rsplit("/", 1)[0]
    program = f'from pathlib import Path; print((Path({folder!r}) / "unattended.json").read_text())'
    assert folder in candidates_of(program)


# ---------------------------------------------------------------- h24 item 14
#
# docs/increments/h24-harness-fixes.md §3 item 14: an interpreter (not a
# shell) whose program the shell may rewrite is `unknown`. For `-c` that is a
# program holding `$` or a backtick; for a heredoc, an unquoted delimiter and a
# body holding `$`, a backtick or a backslash. Any quoting of the delimiter
# (`'EOF'`, `"EOF"`, `\EOF`) turns the shell's expansion off.


def the_simple(line: str, name: str) -> Any:
    simples = shell_scan.parse(line)
    assert simples is not None, line
    [simple] = [s for s in simples if s.argv[:1] == [name]]
    return simple


def test_h24_14_a_double_quoted_c_program_holding_a_substitution_is_unknown() -> None:
    simple = the_simple("python3 -c \"print(open('a$(cat q)').read())\"", "python3")
    assert simple.unknown is True


def test_h24_14_a_shell_c_program_holding_a_substitution_is_parsed_and_the_shell_unread() -> None:
    # The program is still parsed, but a shell whose -c program holds an
    # expansion is fed through a variable, so the shell itself is unread and
    # unknown (docs/increments/h25-stdin-routes.md §3.2, "Through a variable").
    simples = shell_scan.parse('bash -c "echo $(cat q)"')
    assert simples is not None
    assert len(simples) == 4, [s.argv for s in simples]
    assert [s.argv for s in simples if s.unknown] == [["bash", "-c", "echo $(cat q)"]]
    assert [s.argv for s in simples if s.unread] == [["bash", "-c", "echo $(cat q)"]]


#: Design review round 3's B(t), here with an ordinary target.
BACKSLASH_BODY = r"x = 'a\\'; y = 1 #'; open('notes/plain.txt','w').write('X')"


def test_h24_14_an_unquoted_heredoc_with_a_backslash_is_unknown() -> None:
    simple = the_simple(f"python3 - <<EOF\n{BACKSLASH_BODY}\nEOF", "python3")
    assert simple.unknown is True


@pytest.mark.parametrize("opener", ["<<'EOF'", '<<"EOF"', "<<\\EOF"])
def test_h24_14_a_quoted_heredoc_delimiter_keeps_the_program_readable(opener: str) -> None:
    simple = the_simple(f"python3 - {opener}\n{BACKSLASH_BODY}\nEOF", "python3")
    assert simple.program is not None and BACKSLASH_BODY in simple.program
    assert simple.unknown is False


# ---------------------------------------------------------------- h24 item 15
#
# docs/increments/h24-harness-fixes.md §3 item 15: a here-string is its word,
# quotes removed as `word()` removes them, plus a newline, handed to stdin. An
# interpreter with no `-c`/`-e` text and no script runs it; a shell parses it
# as commands. Every stdin body of one command counts (zsh's MULTIOS). `$` or a
# backtick in the word marks the interpreter unknown whatever the quoting.

#: Design review round 4's two writes of an ordinary file.
PYTHON_WRITE = "open('notes/a.txt','w').write('x')"


def test_h24_15_a_python_here_string_is_its_program() -> None:
    simple = the_simple(f'python3 <<< "{PYTHON_WRITE}"', "python3")
    assert simple.program == PYTHON_WRITE + "\n"
    assert "notes/a.txt" in shell_scan.candidates(simple)
    assert simple.unknown is False


def test_h24_15_a_shell_here_string_is_parsed_as_commands() -> None:
    simples = shell_scan.parse('bash <<< "echo x > notes/a.txt"')
    assert simples is not None
    assert [(s.argv, s.writes) for s in simples] == [
        (["bash"], []),
        (["echo", "x"], ["notes/a.txt"]),
    ]


@pytest.mark.parametrize(
    "line",
    [
        f'env python3 <<< "{PYTHON_WRITE}"',
        f'env A=1 python3 <<< "{PYTHON_WRITE}"',
        f'nice -n 5 python3 <<< "{PYTHON_WRITE}"',
        f'timeout 5 python3 <<< "{PYTHON_WRITE}"',
        f'python3 - <<< "{PYTHON_WRITE}"',
    ],
    ids=["env", "env-with-assignment", "nice", "timeout", "dash-script"],
)
def test_h24_15_a_here_string_after_a_wrapper_or_dash_is_the_program(line: str) -> None:
    simple = the_simple(line, "python3")
    assert simple.program == PYTHON_WRITE + "\n"
    assert "notes/a.txt" in shell_scan.candidates(simple)


@pytest.mark.parametrize("line", ['env bash <<< "echo x > notes/a.txt"',
                                  'timeout 5 bash <<< "echo x > notes/a.txt"'])  # fmt: skip
def test_h24_15_a_shell_here_string_after_a_wrapper_is_parsed(line: str) -> None:
    simples = shell_scan.parse(line)
    assert simples is not None
    assert (["echo", "x"], ["notes/a.txt"]) in [(s.argv, s.writes) for s in simples]


@pytest.mark.parametrize(
    "line",
    [
        f"python3 <<< 'print(1)' <<< \"{PYTHON_WRITE}\"",
        f"python3 <<< \"{PYTHON_WRITE}\" <<< 'print(1)'",
        f'python3 <<< "{PYTHON_WRITE}" <<EOF\nprint(1)\nEOF',
        f"python3 <<< 'print(1)' <<EOF\n{PYTHON_WRITE}\nEOF",
    ],
    ids=["write-second", "write-first", "here-string-before-heredoc", "heredoc-write"],
)
def test_h24_15_every_stdin_body_of_one_command_is_in_its_program(line: str) -> None:
    simple = the_simple(line, "python3")
    assert simple.program is not None
    assert "print(1)" in simple.program and PYTHON_WRITE in simple.program
    assert "notes/a.txt" in shell_scan.candidates(simple)


def test_h24_15_two_shell_here_strings_are_both_parsed() -> None:
    simples = shell_scan.parse("zsh <<< 'echo a' <<< 'echo x > notes/a.txt'")
    assert simples is not None
    found = [(s.argv, s.writes) for s in simples]
    assert (["echo", "a"], []) in found and (["echo", "x"], ["notes/a.txt"]) in found


@pytest.mark.parametrize(
    "line",
    [
        "python3 <<< \"print('$(cat q)')\"",
        "python3 <<< 'print(\"$(cat q)\")'",
        "python3 <<< 'print(\"$HOME\")'",
        "python3 <<< 'print(\"`cat q`\")'",
        "python3 <<< \"print('`cat q`')\"",
    ],
    ids=["dollar-paren-double", "dollar-paren-single", "parameter-single", "backtick-single",
         "backtick-double"],
)  # fmt: skip
def test_h24_15_a_here_string_holding_a_substitution_is_unknown(line: str) -> None:
    assert the_simple(line, "python3").unknown is True


def test_h24_15_a_substitution_in_a_second_here_string_is_unknown() -> None:
    simple = the_simple("python3 <<< 'print(1)' <<< 'print(\"$HOME\")'", "python3")
    assert simple.unknown is True


def test_h24_15_a_shell_here_string_holding_a_substitution_is_parsed_not_unknown() -> None:
    simples = shell_scan.parse('bash <<< "echo $(cat q)"')
    assert simples is not None
    assert ["echo", "$(cat q)"] in [s.argv for s in simples]
    assert [s.argv for s in simples if s.unknown] == []


def test_h24_15_a_backslash_in_a_here_string_is_not_a_sign_of_rewriting() -> None:
    simple = the_simple("python3 <<< 'print(\"a\\\\nb\")'", "python3")
    assert simple.program == 'print("a\\\\nb")\n'
    assert simple.unknown is False


def test_h24_15_a_here_string_to_a_non_interpreter_is_data() -> None:
    simple = the_simple("cat <<< x", "cat")
    assert simple.program is None
    assert simple.unknown is False


@pytest.mark.parametrize(
    ("line", "program"),
    [("python3 -c 'print(1)' <<< \"open('CLAUDE.md','w')\"", "print(1)"),
     ("python3 run.py <<< \"open('CLAUDE.md','w')\"", None)],
    ids=["c-program", "script"],
)  # fmt: skip
def test_h24_15_a_c_program_or_a_script_keeps_a_here_string_as_data(
    line: str, program: str | None
) -> None:
    simple = the_simple(line, "python3")
    assert simple.program == program
    assert "CLAUDE.md" not in shell_scan.candidates(simple)


# ---------------------------------------------------------------- h25
#
# docs/increments/h25-stdin-routes.md §3.2-3.4: a runner fed text the scanner
# did not read is `unread`, and so `unknown`; a heredoc or here-string on a
# literal interpreter is its program whatever its options, unless the command
# is a data reader; a process substitution is never a written file.

#: A Python write of CLAUDE.md, the program fed on stdin.
H25_W = "open('CLAUDE.md','w').write('x')"


def test_h25_a_runner_piped_across_a_newline_is_unread() -> None:
    simple = the_simple("echo S |\npython3", "python3")
    assert (simple.unread, simple.unknown) == (True, True)


def test_h25_a_runner_after_a_pipe_and_a_semicolon_is_not_unread() -> None:
    simples = shell_scan.parse("python3 x | tail -1; python3 y")
    assert simples is not None
    pythons = [s for s in simples if s.argv[:1] == ["python3"]]
    assert len(pythons) == 2
    assert [s.unread for s in pythons] == [False, False]


def test_h25_a_runner_inside_a_piped_compound_is_unread() -> None:
    assert the_simple("echo S | if true; then python3; fi", "python3").unread is True


def test_h25_an_interpreter_name_among_a_readers_words_is_not_a_runner() -> None:
    simples = shell_scan.parse("grep -i python CLAUDE.md | head")
    assert simples is not None
    assert [s.argv for s in simples if s.unread] == []


def test_h25_a_here_string_is_read_whatever_the_interpreters_options() -> None:
    simple = the_simple(f'python3 -i -c pass <<< "{H25_W}"', "python3")
    assert simple.program == f"pass\n{H25_W}\n"
    assert "CLAUDE.md" in shell_scan.candidates(simple)


def test_h25_a_process_substitution_target_is_not_a_write() -> None:
    simples = shell_scan.parse("echo x > >(sh)")
    assert simples is not None
    assert [s.writes for s in simples if s.writes] == []


#: §3.4: programs a data reader may run, measured as passing the allowlist.
DATA_PROGRAMS = (
    "import json,sys; print(json.load(sys.stdin)['a'])",
    "import sys; print(len(sys.stdin.read()))",
    "print(1)",
)

#: §3.4: programs outside the allowlist, one rule each.
NOT_DATA_PROGRAMS: dict[str, str] = {
    "sys-modules": "import sys; sys.modules['posix'].system(sys.stdin.read())",
    "loader": "import sys; __loader__.load_module('posix').system(sys.stdin.read())",
    "methodcaller": "import operator; operator.methodcaller('system', 'x')",
    "formatter": "import string; string.Formatter().get_field('0', [1], {})",
    "generator-frame": "g=(x for x in [1]); g.gi_frame.f_builtins['exec']('print(1)')",
    "builtin-bound-in-a-comprehension": "import sys; [0 for exec in []]; exec(sys.stdin.read())",
    "builtin-bound-in-a-lambda": "import sys; f = lambda exec: 0; exec(sys.stdin.read())",
    "builtin-bound-in-dead-code": "import sys\nif 0:\n    exec = print\nexec(sys.stdin.read())",
    "exec-of-input": "exec(input())",
    "import-as": "import json as j; print(j.dumps(1))",
    "from-import": "from json import loads; print(loads('1'))",
    "a-third-module": "import os; print(1)",
    "underscore-name": "import sys; print(_)",
    "attribute-assigned": "import sys; sys.stdout = 1",
    "attribute-deleted": "import sys; del sys.argv",
    "sys-attribute-off-the-list": "import sys; print(sys.path)",
    "method-off-the-list": "x = 'a'; x.encode()",
    "def": "def f():\n    return 1\nprint(f())",
    "class": "class C:\n    pass",
    "lambda-calling-a-builtin-off-the-list": "f = lambda: open('CLAUDE.md', 'w'); f()",
    "syntax-error": "print(",
}  # fmt: skip


@pytest.mark.parametrize("program", DATA_PROGRAMS)
def test_h25_a_program_on_the_allowlist_is_a_data_program(program: str) -> None:
    assert shell_scan.data_program(program) is True


@pytest.mark.parametrize("program", NOT_DATA_PROGRAMS.values(), ids=NOT_DATA_PROGRAMS.keys())
def test_h25_a_program_off_the_allowlist_is_not_a_data_program(program: str) -> None:
    assert shell_scan.data_program(program) is False


#: §3.4: the command shapes around a data program; each is fed by a pipe.
NOT_DATA_READERS: dict[str, str] = {
    "a-wrapper-before-it": "echo S | env python3 x.py",
    "w-option": "echo S | python3 -W ignore x.py",
    "x-option": "echo S | python3 -X dev x.py",
    "i-option": "echo S | python3 -i x.py",
    "m-option": "echo S | python3 -m json.tool",
    "dash-script": "echo S | python3 -",
    "a-long-option": "echo S | python3 --check-hash-based-pycs always x.py",
    "a-script-under-dev": "echo S | python3 /dev/stdin",
    "a-script-starting-with-two-slashes": "echo S | python3 //dev/stdin",
    "a-script-holding-a-dollar": "echo S | python3 $X",
    "a-script-holding-a-glob": "echo S | python3 x*.py",
}

DATA_READERS: dict[str, str] = {
    "a-plain-script": "echo S | python3 x.py",
    "plain-flags-in-a-cluster": "echo S | python3 -BEs x.py",
    "a-data-program": "echo S | python3 -u -c 'import sys; print(sys.stdin.read())'",
    "python3-minor-version": "echo S | python3.12 tools/count_loc.py",
}


@pytest.mark.parametrize("line", NOT_DATA_READERS.values(), ids=NOT_DATA_READERS.keys())
def test_h25_a_fed_python_outside_the_reader_shape_is_unread(line: str) -> None:
    simples = shell_scan.parse(line)
    assert simples is not None
    assert [s.unread for s in simples if s.argv[:1] != ["echo"]] == [True]


@pytest.mark.parametrize("line", DATA_READERS.values(), ids=DATA_READERS.keys())
def test_h25_a_fed_data_reader_is_not_unread(line: str) -> None:
    simples = shell_scan.parse(line)
    assert simples is not None
    assert [s.unread for s in simples] == [False, False]
