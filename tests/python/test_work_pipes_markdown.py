"""`tools/work_pipes.py show --markdown`: the verbose pipe view, in GitHub-flavoured Markdown.

The view is pasted into a reply where Markdown renders: one `###` section per
pipe with a table of its open items, nothing truncated, then the questions for
Ola. Commands run as the harness runs them, a copied script in a temporary
repository, so the real live lists are never read.
"""

from __future__ import annotations

import json
import re
import subprocess
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import git, make_repo, run_script

pytestmark = pytest.mark.harness

HEAD = "| # | Item | Step | Waiting on | Next |"


def ok(repo: Path, script: str, *args: str) -> str:
    result: subprocess.CompletedProcess[str] = run_script(repo, script, "", *args)
    assert result.returncode == 0, f"{script} {' '.join(args)}: {result.stderr}"
    return result.stdout


def pipes(repo: Path, *args: str) -> str:
    return ok(repo, "tools/work_pipes.py", *args)


def rec(item: str, state: str, pipe: str = "core", **extra: Any) -> dict[str, Any]:
    since = datetime.now(UTC) - timedelta(hours=2)
    record: dict[str, Any] = {
        "id": item, "title": f"title of {item}", "state": state, "pipe": pipe,
        "next": "", "branch": f"worktree-{item}", "pr": None, "plain": False,
        "fallback": False, "since": since.isoformat(timespec="seconds"),
    }  # fmt: skip
    return record | extra


def append(repo: Path, *records: dict[str, Any]) -> None:
    common = Path(git(repo, "rev-parse", "--path-format=absolute", "--git-common-dir").strip())
    for record in records:
        path = common / "pipes" / f"{record['pipe']}.jsonl"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a") as out:
            out.write(json.dumps(record) + "\n")


def note(repo: Path, name: str) -> None:
    path = repo / ".claude" / "current-task" / f"developer-{name}-190000.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ask\n")


def ask(repo: Path, text: str, *more: str) -> None:
    ok(repo, "tools/questions.py", "add", text, "--default", "wait", *more)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    root = make_repo(tmp_path / "repo")
    (root / ".claude" / "pipes").mkdir(parents=True)
    (root / ".claude" / "pipes" / "core.toml").write_text('name = "core"\npurpose = "p"\n')
    (root / ".claude" / "pipes" / "lsm.toml").write_text('name = "lsm"\npurpose = "q"\n')
    pipes(root, "priority", "core=1", "lsm=2", "--words", "core first")
    return root


def section(lines: list[str], heading: str) -> list[str]:
    """The lines under `heading` up to the next `###` heading."""
    start = lines.index(heading) + 1
    end = next((i for i in range(start, len(lines)) if lines[i].startswith("###")), len(lines))
    return lines[start:end]


def rows(lines: list[str]) -> list[list[str]]:
    """The table rows of a section, split on unescaped `|`, the head, rule and number dropped."""
    found = [line for line in lines if line.startswith("|")]
    assert found[0] == HEAD and set(found[1]) <= set("|- :"), found
    split = [line.replace("\\|", "\0").strip("|").split("|") for line in found[2:]]
    cells = [[cell.strip().replace("\0", "\\|") for cell in row] for row in split]
    assert all(re.fullmatch(r"\d+\.\d+", row[0]) for row in cells), cells
    return [row[1:] for row in cells]


def test_one_section_per_pipe_with_open_items_in_flight_first(repo: Path) -> None:
    long_title = "a title well past thirty characters, " * 4 + "kept whole"
    append(repo, rec("q1", "queued", title=long_title), rec("r1", "red"),
           rec("m1", "merged"), rec("d1", "dropped"), rec("p1", "parked"),
           rec("g1", "green", next="run the suite", pr=12),
           rec("l1", "queued", pipe="lsm"))  # fmt: skip
    lines = pipes(repo, "show", "--markdown").splitlines()
    core, lsm = rows(section(lines, "### core")), rows(section(lines, "### lsm"))
    assert [r[0] for r in core] == ["**g1** title of g1", "**r1** title of r1",
                                     f"**q1** {long_title.strip()}"]  # fmt: skip
    assert [r[1] for r in core] == ["green", "red", "queued"]
    assert "run the suite" in core[0][3] and "#12" in core[0][3]
    assert [r[0] for r in lsm] == ["**l1** title of l1"]
    assert lines.index("### core") < lines.index("### lsm")


def test_waiting_on_names_the_question_the_agent_or_a_free_writer(repo: Path) -> None:
    append(repo, rec("a", "queued"), rec("b", "designing"), rec("c", "red"), rec("e", "queued"))
    ask(repo, "Which name for the command?", "--choices", "show|list")
    ask(repo, "Split it in two?", "--item", "b")
    ask(repo, "Merge it now?", "--item", "e")  # a second yes/no question: held
    note(repo, "c")
    lines = pipes(repo, "show", "--markdown").splitlines()
    waits = {r[0].split()[0]: r[2] for r in rows(section(lines, "### core"))}
    assert waits["**a**"] == "a free writer"
    assert waits["**b**"] == "Ola (question 2)"
    assert waits["**c**"] == "a running agent"
    assert waits["**e**"].startswith("Ola (")


def test_the_first_line_and_the_questions_section(repo: Path) -> None:
    append(repo, rec("r1", "red"), rec("q1", "queued"))
    ask(repo, "Which name for the command?", "--choices", "show|list")
    ask(repo, "Split it in two?", "--item", "q1")
    ask(repo, "Merge it now?")  # held behind the yes/no question above
    lines = pipes(repo, "show", "--markdown").splitlines()
    assert lines[0].startswith("In flight 1 of 3") and "**r1**" in lines[0], lines[0]
    questions = section(lines, "### Questions for Ola")
    listed = run_script(repo, "tools/questions.py", "", "list").stdout.splitlines()
    numbered = [line for line in listed if line[0].isdigit()]
    assert numbered == [line for line in questions if line[:1].isdigit()]
    assert len(numbered) == 2
    assert any("1 more" in line and "held" in line for line in questions), questions


def test_a_pipe_in_a_cell_is_escaped_and_plain_show_is_unchanged(repo: Path) -> None:
    append(repo, rec("x", "queued", title="this | that", next="a|b"))
    plain = pipes(repo, "show")
    lines = pipes(repo, "show", "--markdown").splitlines()
    [row] = rows(section(lines, "### core"))
    assert row[0] == "**x** this \\| that" and row[3] == "a\\|b"
    assert pipes(repo, "show") == plain and not plain.startswith("In flight")
    assert "###" not in plain
