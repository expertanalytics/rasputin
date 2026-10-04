"""`tools/rule_sizes.py`: rule text in words, against the last retrospective (h8, test 9).

Spec: `docs/increments/h8-window-and-recap.md` §3.6. The reference is derived
from git, never stored: the newest commit reachable from HEAD that added a
dated retrospective file. Files removed since are found at that commit with
`git ls-tree`, so they show as `removed (was n)`. The pure parts (`words`,
`table`) are tested on strings; the git half end to end, by running the
fixture's copy of the script, so its root is the fixture repository.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from harness_fixtures import REAL, Tool, clean_env, git, make_repo, run_script

rule_sizes = Tool("rule_sizes")

LABEL = "2026-10-01-first.md (abc1234)"


def row(words: int, change: str, path: str) -> str:
    """One file line of §3.6: `  <words:5> <change:>6> <path>`."""
    return f"  {words:5} {change:>6} {path}"


# ---------------------------------------------------------------- words


ASCII_FIXTURE = "  # Rules\n\nOne  two\tthree\n- four, five.\n\n\tsix\r\nseven-eight `nine`\n"


def test_words_counts_what_wc_w_counts_on_ascii_text() -> None:
    wc = shutil.which("wc")
    if wc is None:  # pragma: no cover - wc is POSIX
        pytest.skip("wc is not on PATH")
    counted = subprocess.run(
        [wc, "-w"], input=ASCII_FIXTURE, capture_output=True, text=True, check=True
    )
    assert rule_sizes.words(ASCII_FIXTURE) == int(counted.stdout.split()[0])
    assert rule_sizes.words(ASCII_FIXTURE) == 11


def test_words_of_empty_or_blank_text_is_zero() -> None:
    assert rule_sizes.words("") == 0
    assert rule_sizes.words(" \n\t\n") == 0


def test_rule_files_are_the_four_governed_rule_files_in_order() -> None:
    assert rule_sizes.RULE_FILES == (
        "CLAUDE.md",
        ".claude/REQUIRED-READING.md",
        "docs/increments/README.md",
        "docs/PRINCIPLES.md",
    )


# ---------------------------------------------------------------- table, pure

NOW = {"CLAUDE.md": 112, ".claude/REQUIRED-READING.md": 97, ".claude/agents/b.md": 50,
       ".claude/agents/new.md": 5}  # fmt: skip
THEN = {"CLAUDE.md": 100, ".claude/REQUIRED-READING.md": 100, ".claude/agents/b.md": 50,
        ".claude/agents/old.md": 7}  # fmt: skip


def test_table_shows_growth_cuts_new_unchanged_and_removed_files() -> None:
    lines = rule_sizes.table(NOW, THEN, LABEL)
    assert lines[0] == f"Rule text in words, change since {LABEL}:"
    body = lines[1:-1]
    assert row(112, "+12", "CLAUDE.md") in body
    assert row(97, "-3", ".claude/REQUIRED-READING.md") in body
    assert row(50, "0", ".claude/agents/b.md") in body
    assert row(5, "new", ".claude/agents/new.md") in body
    assert len(body) == 5
    # §8.13: words column 0, change text unpadded, after every file present now.
    assert body[4] == f"  {0:5} removed (was 7) .claude/agents/old.md"


def test_table_total_is_the_change_over_both_sides_including_removed_files() -> None:
    # now 264, then 257 (the removed file's 7 words count at the reference).
    assert rule_sizes.table(NOW, THEN, LABEL)[-1] == row(264, "+7", "total")


def test_table_lists_files_present_now_in_the_order_given() -> None:
    lines = rule_sizes.table(NOW, THEN, LABEL)
    shown = [line.rsplit(" ", 1)[1] for line in lines[1:-1]]
    assert [p for p in shown if p in NOW] == list(NOW)


def test_table_with_a_smaller_total_shows_a_negative_total_change() -> None:
    lines = rule_sizes.table({"CLAUDE.md": 90}, {"CLAUDE.md": 100}, LABEL)
    assert lines == [
        f"Rule text in words, change since {LABEL}:",
        row(90, "-10", "CLAUDE.md"),
        row(90, "-10", "total"),
    ]


def test_table_without_a_reference_has_no_change_column() -> None:
    assert rule_sizes.table({"CLAUDE.md": 112, "docs/PRINCIPLES.md": 40}, None, LABEL) == [
        "Rule text in words (no retrospective found):",
        f"  {112:5} CLAUDE.md",
        f"  {40:5} docs/PRINCIPLES.md",
        f"  {152:5} total",
    ]


# ---------------------------------------------------------------- end to end
#
# The fixture: the four rule files, two agents, a nested skill file, committed;
# a retrospective added in its own commit (the reference); then CLAUDE.md grown
# by 12 words, REQUIRED-READING cut by 3, one agent added and one deleted.

RETRO = "docs/retrospectives/2026-10-01-first.md"
LATER_RETRO = "docs/retrospectives/2026-10-02-second.md"

STARTING = {
    "CLAUDE.md": "# rules\n" + "word " * 18,  # 20 words
    ".claude/REQUIRED-READING.md": "read " * 30,
    "docs/increments/README.md": "loop " * 15,
    "docs/PRINCIPLES.md": "principle " * 10,
    ".claude/agents/tester.md": "test " * 8,
    ".claude/agents/old.md": "old " * 7,
    ".claude/skills/py/SKILL.md": "skill " * 6,
    ".claude/skills/py/ref/deep.md": "deep " * 4,
}


def write(repo: Path, relative: str, text: str) -> None:
    (repo / relative).parent.mkdir(parents=True, exist_ok=True)
    (repo / relative).write_text(text)


def commit_all(repo: Path, message: str) -> str:
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message)
    return git(repo, "log", "-1", "--format=%h").strip()


@pytest.fixture
def sized(tmp_path: Path) -> tuple[Path, str]:
    """The fixture repository after the edits; and the reference commit's hash."""
    repo = make_repo(tmp_path.resolve() / "repo")
    for relative, text in STARTING.items():
        write(repo, relative, text)
    commit_all(repo, "rule files")
    write(repo, RETRO, "# A retrospective\n")
    reference = commit_all(repo, "retrospective")
    write(repo, "CLAUDE.md", STARTING["CLAUDE.md"] + "more " * 12)
    write(repo, ".claude/REQUIRED-READING.md", "read " * 27)
    write(repo, ".claude/agents/new.md", "new " * 5)
    (repo / ".claude/agents/old.md").unlink()
    # An undated file under docs/retrospectives/ is not a retrospective.
    write(repo, "docs/retrospectives/next.md", "next items\n")
    commit_all(repo, "edits since the retrospective")
    return repo, reference


def sizes(repo: Path) -> list[str]:
    result = run_script(repo, "tools/rule_sizes.py", "")
    assert result.returncode == 0, result.stderr
    return result.stdout.splitlines()


def test_the_script_prints_the_table_against_the_last_retrospective(
    sized: tuple[Path, str],
) -> None:
    repo, reference = sized
    lines = sizes(repo)
    assert lines[0] == f"Rule text in words, change since 2026-10-01-first.md ({reference}):"
    assert row(32, "+12", "CLAUDE.md") in lines
    assert row(27, "-3", ".claude/REQUIRED-READING.md") in lines
    assert row(15, "0", "docs/increments/README.md") in lines
    assert row(10, "0", "docs/PRINCIPLES.md") in lines
    assert row(5, "new", ".claude/agents/new.md") in lines
    assert row(8, "0", ".claude/agents/tester.md") in lines
    assert row(6, "0", ".claude/skills/py/SKILL.md") in lines
    assert row(4, "0", ".claude/skills/py/ref/deep.md") in lines
    [removed] = [line for line in lines if line.endswith(" .claude/agents/old.md")]
    assert removed.endswith(" removed (was 7) .claude/agents/old.md")
    # now 32+27+15+10+5+8+6+4 = 107; then 20+30+15+10+8+7+6+4 = 100.
    assert lines[-1] == row(107, "+7", "total")


def test_the_script_orders_the_rule_files_then_agents_then_skills(
    sized: tuple[Path, str],
) -> None:
    repo, _ = sized
    paths = [line.rsplit(" ", 1)[1] for line in sizes(repo)[1:-1]]
    present = [p for p in paths if p != ".claude/agents/old.md"]
    assert present == [
        "CLAUDE.md",
        ".claude/REQUIRED-READING.md",
        "docs/increments/README.md",
        "docs/PRINCIPLES.md",
        ".claude/agents/new.md",
        ".claude/agents/tester.md",
        ".claude/skills/py/SKILL.md",
        ".claude/skills/py/ref/deep.md",
    ]


def test_the_script_lists_only_rule_files(sized: tuple[Path, str]) -> None:
    repo, _ = sized
    text = "\n".join(sizes(repo))
    assert "notes.txt" not in text
    assert "retrospectives" not in text.split("\n", 1)[1]
    assert ".claude/hooks/" not in text


def test_a_later_retrospective_moves_the_reference(sized: tuple[Path, str]) -> None:
    repo, _ = sized
    write(repo, LATER_RETRO, "# Another\n")
    later = commit_all(repo, "second retrospective")
    lines = sizes(repo)
    assert lines[0] == f"Rule text in words, change since 2026-10-02-second.md ({later}):"
    assert row(32, "0", "CLAUDE.md") in lines
    assert not any("old.md" in line for line in lines)
    assert lines[-1] == row(107, "0", "total")


def test_without_a_retrospective_the_table_has_no_reference(tmp_path: Path) -> None:
    repo = make_repo(tmp_path.resolve() / "repo")
    write(repo, ".claude/agents/tester.md", "test " * 8)
    commit_all(repo, "an agent")
    lines = sizes(repo)
    assert lines[0] == "Rule text in words (no retrospective found):"
    assert f"  {2:5} CLAUDE.md" in lines
    assert f"  {8:5} .claude/agents/tester.md" in lines
    assert lines[-1] == f"  {10:5} total"


def test_an_uncommitted_edit_counts_because_now_is_the_working_tree(
    sized: tuple[Path, str],
) -> None:
    # §8.14: the committed +12 plus an uncommitted 3-word append.
    repo, _ = sized
    with (repo / "CLAUDE.md").open("a") as out:
        out.write("one two three\n")
    lines = sizes(repo)
    assert row(35, "+15", "CLAUDE.md") in lines
    assert lines[-1] == row(110, "+10", "total")


def test_without_a_repository_the_counts_print_with_a_git_failed_heading(
    tmp_path: Path,
) -> None:
    # §8.19: git fails, so no reference; the working-tree counts still print.
    plain = tmp_path.resolve() / "plain"
    (plain / "tools").mkdir(parents=True)
    source = REAL / "tools" / "rule_sizes.py"
    assert source.exists(), "tools/rule_sizes.py is missing from the checkout"
    shutil.copy2(source, plain / "tools" / "rule_sizes.py")
    (plain / "CLAUDE.md").write_text("a b c\n")
    result = subprocess.run(
        [sys.executable, str(plain / "tools" / "rule_sizes.py")],
        input="",
        capture_output=True,
        text=True,
        cwd=plain,
        env=clean_env(),
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [
        "Rule text in words (no reference: git failed):",
        f"  {3:5} CLAUDE.md",
        f"  {3:5} total",
    ]
