"""`tools/check_citations.py`: pinned, ID and heading-name citations, governed files, scan set.

The design is `docs/increments/h2-citation-pinning.md`; its section 6 is the
checklist this module follows, and each block below names the test numbers it
covers. The existing `test_check_citations.py` keeps the `legacy/` behaviour.

Black-box, like the existing suite: the tool finds its repository from its own
path, so each test copies it into a scratch repository, commits, runs it as a
subprocess and asserts on the exit status, the report section a citation lands
in, and the citing location and cited token inside that section. Nothing here
reads the live tree's citations.

This module's own comments and docstrings are scanned by the tool in the real
tree, so every line reference in them is a placeholder (`testing.md:<n>`);
the citations under test live in string literals, which are not scanned.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import textwrap
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import pytest

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

TOOL = Path(__file__).resolve().parents[2] / "tools" / "check_citations.py"

# A governed target long enough for every planted line reference.
TESTING_MD = "# Testing\n\n## Frameworks [partly live]\n\nline four\nline five\n"
PLANT = "testing.md:3"


@dataclass(frozen=True)
class Run:
    code: int
    out: str

    def section(self, name: str) -> str:
        """The body of the `== <name> (<n>)` section, header included; '' when absent."""
        kept: list[str] = []
        inside = False
        for line in self.out.splitlines():
            stripped = line.strip()
            if stripped.startswith("== "):
                inside = stripped.startswith(f"== {name} (")
            if inside:
                kept.append(line)
        return "\n".join(kept)

    def count(self, name: str) -> int:
        """The `<n>` of the `== <name> (<n>)` header, 0 when the section is absent."""
        match = re.search(rf"^\s*== {re.escape(name)} \((\d+)\)", self.out, re.MULTILINE)
        return int(match.group(1)) if match else 0


class Scratch:
    """A git repository on `master` holding a copy of the tool at `tools/`."""

    def __init__(self, root: Path) -> None:
        self.root = root

    def git(self, *args: str) -> str:
        result = subprocess.run(
            ["git", "-C", str(self.root), *args], check=True, capture_output=True, text=True
        )
        return result.stdout.strip()

    def write(self, rel: str, text: str) -> None:
        path = self.root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def commit(self, message: str = "change") -> str:
        """Commit everything and return the new commit's full sha."""
        self.git("add", "-A")
        self.git("commit", "-q", "--allow-empty", "-m", message)
        return self.git("rev-parse", "HEAD")

    def cite(self, *lines: str) -> None:
        """Commit `docs/record.md`: a heading, a blank line, then `lines` from line 3."""
        self.write("docs/record.md", "# Record\n\n" + "".join(f"{line}\n" for line in lines))
        self.commit("cite")

    def check(self, *paths: str, base: str = "master") -> Run:
        """Run the tool; with no `paths`, its default scan set (no `--paths`)."""
        argv = [sys.executable, str(self.root / "tools" / "check_citations.py"), "--base", base]
        if paths:
            argv += ["--paths", *paths]
        result = subprocess.run(argv, cwd=self.root, capture_output=True, text=True, timeout=60)
        return Run(result.returncode, result.stdout + result.stderr)


def _init(root: Path) -> Scratch:
    root.mkdir(parents=True)
    repo = Scratch(root)
    repo.git("init", "-q", "-b", "master")
    for key, value in (
        ("user.name", "t"),
        ("user.email", "t@t"),
        ("commit.gpgsign", "false"),
        ("tag.gpgsign", "false"),
    ):
        repo.git("config", key, value)
    (root / "tools").mkdir()
    shutil.copy(TOOL, root / "tools" / "check_citations.py")
    repo.write("docs/.keep.md", "# Keep\n")
    repo.commit("tool")
    return repo


@pytest.fixture
def repo(tmp_path: Path) -> Scratch:
    return _init(tmp_path / "repo")


def _lines(n: int) -> str:
    return "".join(f"line {i}\n" for i in range(1, n + 1))


# The citing doc starts its citations at line 3 (`Scratch.cite`).
FIRST = "docs/record.md:3"


# -- Pinned: tests 1-9 --------------------------------------------------------


@pytest.fixture
def pinned(repo: Scratch) -> tuple[Scratch, str]:
    """`notes/target.md` has 5 lines at the returned sha, then a branch deletes it."""
    repo.write("notes/target.md", _lines(5))
    sha = repo.commit("target")
    repo.git("checkout", "-q", "-b", "feature")
    return repo, sha


@pytest.mark.parametrize("change", ["rewritten", "deleted"])
@pytest.mark.parametrize(
    ("rev_of", "lines"),
    [(lambda sha: sha, "4"), (lambda sha: sha[:7], "2-5")],
    ids=["full-sha", "short-sha-range"],
)
def test_a_sha_pin_resolves_at_its_commit_and_is_never_at_risk(
    pinned: tuple[Scratch, str], rev_of: Callable[[str], str], lines: str, change: str
) -> None:
    repo, sha = pinned
    if change == "deleted":
        (repo.root / "notes" / "target.md").unlink()
    else:
        repo.write("notes/target.md", "one line now\n")
    repo.commit("the branch edits the pinned file")
    rev = rev_of(sha)
    repo.cite(f"`notes/target.md@{rev}:{lines}`")

    run = repo.check("docs")

    assert run.code == 0, run.out
    assert run.section("at risk") == "", run.out
    assert run.section("broken") == "", run.out


@pytest.mark.parametrize(
    ("tag", "annotated"),
    [("v1", False), ("rel-1", True), ("archive/2026-09", True), ("light/one", False)],
    ids=["lightweight", "annotated", "annotated-with-slash", "lightweight-with-slash"],
)
def test_a_tag_pin_resolves_at_the_tag(
    pinned: tuple[Scratch, str], tag: str, annotated: bool
) -> None:
    repo, sha = pinned
    if annotated:
        repo.git("tag", "-a", tag, "-m", "pin", sha)
    else:
        repo.git("tag", tag, sha)
    (repo.root / "notes" / "target.md").unlink()
    repo.cite(f"`notes/target.md@{tag}:5`")

    run = repo.check("docs")

    assert run.code == 0, run.out


@pytest.mark.parametrize("lines", ["6", "0", "2-6"], ids=["past-end", "line-0", "range-past"])
def test_a_pin_past_the_end_of_the_file_at_its_rev_is_broken(
    pinned: tuple[Scratch, str], lines: str
) -> None:
    repo, sha = pinned
    # The working tree has plenty of lines: only the file at the rev counts.
    repo.write("notes/target.md", _lines(50))
    repo.cite(f"`notes/target.md@{sha}:{lines}`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    broken = run.section("broken")
    assert FIRST in broken
    assert f"notes/target.md@{sha}:{lines.split('-')[0]}" in broken


def test_a_pin_to_a_path_absent_at_the_rev_is_broken_with_no_working_tree_fallback(
    pinned: tuple[Scratch, str],
) -> None:
    repo, sha = pinned
    repo.write("notes/later.md", _lines(5))
    repo.cite(f"`notes/later.md@{sha}:1`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    broken = run.section("broken")
    assert FIRST in broken
    assert f"notes/later.md@{sha}:1" in broken


def test_a_pinned_basename_is_not_searched_for(pinned: tuple[Scratch, str]) -> None:
    repo, sha = pinned  # at the rev the file is `notes/target.md`, never `target.md`
    repo.cite(f"`target.md@{sha}:1`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    assert f"target.md@{sha}:1" in run.section("broken")


@pytest.mark.parametrize("kind", ["absent", "blob", "tree"])
def test_a_hex_rev_that_is_not_a_present_commit_is_broken(
    pinned: tuple[Scratch, str], kind: str
) -> None:
    repo, sha = pinned
    rev = {
        "absent": "0123456789abcdef0123456789abcdef01234567",
        "blob": repo.git("rev-parse", f"{sha}:notes/target.md"),
        "tree": repo.git("rev-parse", f"{sha}^{{tree}}"),
    }[kind]
    repo.cite(f"`notes/target.md@{rev}:1`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    broken = run.section("broken")
    assert FIRST in broken
    assert rev in broken
    # A full clone: the report must not send the reader after a fetch-depth.
    assert "shallow" not in broken


@pytest.mark.parametrize("rev", ["master", "feature", "HEAD", "topic"])
def test_a_branch_or_head_is_not_a_pin(pinned: tuple[Scratch, str], rev: str) -> None:
    repo, sha = pinned
    # `topic` exists only as a remote-tracking branch.
    repo.git("update-ref", "refs/remotes/origin/topic", sha)
    repo.cite(f"`notes/target.md@{rev}:1`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    broken = run.section("broken")
    assert FIRST in broken
    assert f"notes/target.md@{rev}:1" in broken
    assert "not a pin" in broken


def test_a_name_that_is_both_a_tag_and_a_branch_pins_as_the_tag(
    pinned: tuple[Scratch, str],
) -> None:
    repo, sha = pinned
    repo.git("tag", "both", sha)
    repo.git("branch", "both", sha)
    repo.cite("`notes/target.md@both:5`")

    run = repo.check("docs")

    assert run.code == 0, run.out


@pytest.mark.parametrize("pin", ["sha", "tag"])
def test_in_a_shallow_clone_a_missing_rev_says_the_clone_is_shallow(
    tmp_path: Path, pin: str
) -> None:
    origin = _init(tmp_path / "origin")
    origin.write("notes/target.md", _lines(5))
    old = origin.commit("target")
    origin.git("tag", "-a", "v0", "-m", "old", old)
    rev = old if pin == "sha" else "v0"
    origin.cite(f"`notes/target.md@{rev}:1`")
    clone = tmp_path / "clone"
    subprocess.run(
        ["git", "clone", "-q", "--depth", "1", f"file://{origin.root}", str(clone)],
        check=True,
        capture_output=True,
    )
    shallow = Scratch(clone)
    assert shallow.git("rev-parse", "--is-shallow-repository") == "true"

    run = shallow.check("docs")

    assert run.code == 1, run.out
    broken = run.section("broken")
    assert FIRST in broken
    assert rev in broken
    assert "shallow" in broken


def test_a_pinned_citation_is_one_entry_not_also_an_unpinned_one(
    pinned: tuple[Scratch, str],
) -> None:
    repo, sha = pinned
    # The tag name ends in `.md`, so an unpinned match inside the pinned span
    # would read `release/notes.md` plus the line as a second, broken citation.
    repo.git("tag", "release/notes.md", sha)
    repo.cite("`notes/target.md@release/notes.md:2`", f"`notes/target.md@{sha}:9`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    assert run.count("broken") == 1, run.out
    assert "docs/record.md:4" in run.section("broken")
    assert "docs/record.md:3" not in run.out


# -- Legacy: test 10 ----------------------------------------------------------


@pytest.fixture
def archived(repo: Scratch) -> Scratch:
    repo.write("legacy/x.py", _lines(5))
    repo.commit("with legacy")
    repo.git("tag", "-a", "legacy-archive", "-m", "archive")
    shutil.rmtree(repo.root / "legacy")
    repo.commit("legacy leaves")
    return repo


@pytest.mark.parametrize(
    ("citation", "code"),
    [
        ("legacy/x.py@legacy-archive:5", 0),
        ("legacy/x.py:5", 0),
        ("legacy/x.py@legacy-archive:6", 1),
        ("legacy/x.py:6", 1),
    ],
    ids=["pinned-resolves", "shorthand-resolves", "pinned-past-end", "shorthand-past-end"],
)
def test_the_legacy_shorthand_resolves_as_the_pinned_form(
    archived: Scratch, citation: str, code: int
) -> None:
    archived.cite(f"`{citation}`")

    run = archived.check("docs")

    assert run.code == code, run.out
    if code:
        assert FIRST in run.section("broken")


# -- ID and heading-name: tests 11-15 -----------------------------------------

RULES_MD = """\
# Rules

## 1. Intro

## 2. Constraints

## 3. Testing

### 3.2 Citation stability

### C. Data source

## E1 — Principle one

### QW1. Quick win

## §4 Governed

## Appendix

### A tie point
"""

# ID -> the edit to `RULES_MD` that leaves it undeclared.
UNDECLARE = {
    "2": ("## 2. Constraints", "## Constraints"),
    "3.2": ("### 3.2 Citation stability", "### Citation stability"),
    "3C": ("## 3. Testing", "## Testing"),
    "C": ("### C. Data source", "### Data source"),
    "E1": ("## E1 — Principle one", "## Principle one"),
    "QW1": ("### QW1. Quick win", "### Quick win"),
    "4": ("## §4 Governed", "## Governed"),
}


@pytest.fixture
def rules(repo: Scratch) -> Scratch:
    repo.write("docs/rules.md", RULES_MD)
    repo.commit("rules")
    return repo


@pytest.mark.parametrize(
    "form",
    [
        "`docs/rules.md` §{id}",
        "docs/rules.md §{id}",
        "`docs/rules.md`\t§{id}.",
        "(`rules.md` §{id})",
    ],
    ids=["backticked", "bare", "tab-and-period", "basename-in-parens"],
)
@pytest.mark.parametrize("ident", list(UNDECLARE))
def test_an_id_citation_resolves_when_a_heading_declares_it(
    rules: Scratch, ident: str, form: str
) -> None:
    rules.cite("See " + form.format(id=ident) + " for this.")

    run = rules.check("docs")

    assert run.code == 0, run.out


@pytest.mark.parametrize("ident", list(UNDECLARE))
def test_an_id_citation_is_broken_when_no_heading_declares_it(rules: Scratch, ident: str) -> None:
    old, new = UNDECLARE[ident]
    rules.write("docs/rules.md", RULES_MD.replace(old, new))
    rules.cite(f"See `docs/rules.md` §{ident} for this.")

    run = rules.check("docs")

    assert run.code == 1, run.out
    broken = run.section("broken")
    assert FIRST in broken
    assert "docs/rules.md" in broken
    assert ident in broken


@pytest.mark.parametrize(
    "ident",
    ["A", "5", "9.9"],
    ids=["letter-without-period", "undeclared-number", "undeclared-subsection"],
)
def test_an_id_no_heading_declares_is_broken(rules: Scratch, ident: str) -> None:
    rules.cite(f"See `docs/rules.md` §{ident} here.")

    run = rules.check("docs")

    assert run.code == 1, run.out
    assert FIRST in run.section("broken")


def test_an_id_citation_to_a_missing_file_is_broken(rules: Scratch) -> None:
    rules.cite("See `docs/nowhere.md` §1 here.")

    run = rules.check("docs")

    assert run.code == 1, run.out
    assert FIRST in run.section("broken")


@pytest.mark.parametrize("fence", ["```", "~~~", "````"])
def test_a_heading_inside_a_fenced_block_declares_nothing(rules: Scratch, fence: str) -> None:
    rules.write(
        "docs/rules.md",
        RULES_MD + f"\n{fence}text\n## 7. Fenced\n## Fenced name\n{fence}\n",
    )
    rules.cite("See `docs/rules.md` §7 here.", "See `docs/rules.md`, *Fenced name* here.")

    run = rules.check("docs")

    assert run.code == 1, run.out
    broken = run.section("broken")
    assert "docs/record.md:3" in broken
    assert "docs/record.md:4" in broken


def test_section_words_and_spaced_section_signs_are_not_citations(rules: Scratch) -> None:
    rules.cite(
        "See `docs/rules.md` § Nonexistent heading.",
        "See `docs/rules.md` section 9 for more.",
        "Placeholders are not citations: `docs/rules.md:<n>`, `docs/rules.md@<rev>:<n>`.",
    )

    run = rules.check("docs")

    assert run.code == 0, run.out


@pytest.mark.parametrize(
    ("heading", "resolves"),
    [
        ("## Name", True),
        ("## Name [tag]", True),
        ("## Name: rest", True),
        ("### Name — aside", True),
        ("## Name (planned)", True),
        ("## Named", False),
        ("## name", False),
        ("## The Name", False),
    ],
)
def test_a_heading_name_citation_matches_the_heading_or_its_prefix(
    repo: Scratch, heading: str, resolves: bool
) -> None:
    repo.write("docs/x.md", f"# X\n\n{heading}\n\ntext\n")
    repo.cite("As `docs/x.md`, *Name* says.")

    run = repo.check("docs")

    assert run.code == (0 if resolves else 1), run.out
    if not resolves:
        assert FIRST in run.section("broken")


def test_a_heading_name_citation_to_a_missing_file_is_broken(repo: Scratch) -> None:
    repo.cite("As `docs/nowhere.md`, *Name* says.")

    run = repo.check("docs")

    assert run.code == 1, run.out
    assert FIRST in run.section("broken")


@pytest.fixture
def two_notes(repo: Scratch) -> Scratch:
    """Two `NOTES.md` under different directories, and a directory name that is
    a suffix of another (`a/` and `ba/`)."""
    for rel in ("one/NOTES.md", "two/NOTES.md", "a/N.md", "ba/N.md"):
        repo.write(f"src_python/{rel}", "# Notes\n\n## 1. First\n")
    repo.commit("notes")
    return repo


@pytest.mark.parametrize(
    "citation",
    ["`one/NOTES.md:3`", "`one/NOTES.md` §1", "`one/NOTES.md`, *1. First*", "`a/N.md:3`"],
    ids=["line", "id", "heading-name", "component-not-substring"],
)
def test_a_partial_path_resolves_by_its_path_components(two_notes: Scratch, citation: str) -> None:
    two_notes.cite(citation)

    run = two_notes.check("docs")

    assert run.code == 0, run.out


def test_a_bare_basename_with_several_candidates_stays_ambiguous(two_notes: Scratch) -> None:
    two_notes.cite("`NOTES.md:3`")

    run = two_notes.check("docs")

    assert run.code == 1, run.out
    assert FIRST in run.section("broken")


# -- Governed: tests 16-17 ----------------------------------------------------

GOVERNED_MD = "# Rules\n\n## 1. One\n\n## Frameworks\n\ntext\n"
GOVERNED = (
    "CLAUDE.md",
    "testing.md",
    "docs/PRINCIPLES.md",
    ".claude/REQUIRED-READING.md",
    "docs/increments/README.md",
    ".claude/agents/x.md",
    ".claude/skills/y/SKILL.md",
    ".claude/hooks/z.py",
    "tools/check_w.py",
)


@pytest.fixture
def governed(repo: Scratch) -> tuple[Scratch, str]:
    for rel in GOVERNED:
        repo.write(rel, GOVERNED_MD if rel.endswith(".md") else _lines(7))
    repo.write("docs/free.md", _lines(7))
    sha = repo.commit("rule files")
    return repo, sha


@pytest.mark.parametrize(
    "citation",
    [f"{rel}:3" for rel in GOVERNED] + ["REQUIRED-READING.md:3", "check_w.py:2-3"],
    ids=[*GOVERNED, "basename-into-governed", "basename-range-into-tools-check"],
)
def test_an_unpinned_line_citation_into_a_governed_file_fails(
    governed: tuple[Scratch, str], citation: str
) -> None:
    repo, _ = governed
    repo.cite(f"`{citation}`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    unpinned = run.section("unpinned")
    assert FIRST in unpinned
    assert citation.split("-")[0] in unpinned
    assert FIRST not in run.section("broken")


def test_unpinned_outranks_at_risk_for_a_governed_file_the_branch_edits(
    governed: tuple[Scratch, str],
) -> None:
    repo, _ = governed
    repo.git("checkout", "-q", "-b", "feature")
    repo.write("CLAUDE.md", GOVERNED_MD + "more\n")
    repo.commit("edit a rule file")
    repo.cite("`CLAUDE.md:3`")

    run = repo.check("docs")

    assert run.code == 1, run.out
    assert FIRST in run.section("unpinned")
    assert FIRST not in run.section("at risk")


@pytest.mark.parametrize("rel", GOVERNED)
def test_a_governed_file_cited_pinned_by_id_or_by_heading_passes(
    governed: tuple[Scratch, str], rel: str
) -> None:
    repo, sha = governed
    citations = [f"`{rel}@{sha}:3`"]
    if rel.endswith(".md"):
        citations += [f"`{rel}` §1", f"`{rel}`, *Frameworks*"]
    repo.cite(*citations)

    run = repo.check("docs")

    assert run.code == 0, run.out
    assert run.section("unpinned") == ""


def test_an_unpinned_citation_into_an_edited_ungoverned_file_is_only_at_risk(
    governed: tuple[Scratch, str],
) -> None:
    repo, _ = governed
    repo.git("checkout", "-q", "-b", "feature")
    repo.write("docs/free.md", _lines(8))
    repo.commit("edit an ungoverned file")
    repo.cite("`docs/free.md:3`")

    run = repo.check("docs")

    assert run.code == 0, run.out
    assert FIRST in run.section("at risk")
    assert run.section("unpinned") == ""


def test_an_undiffable_base_still_warns_and_does_not_fail(repo: Scratch) -> None:
    repo.write("docs/free.md", _lines(3))
    repo.cite("`docs/free.md:3`")

    run = repo.check("docs", base="no-such-base")

    assert run.code == 0, run.out
    assert "no-such-base" in run.out


# -- Scan set: tests 18-22 ----------------------------------------------------


@pytest.fixture
def plantable(repo: Scratch) -> Scratch:
    repo.write("testing.md", TESTING_MD)
    repo.commit("governed target")
    return repo


def _dedent(text: str) -> str:
    return textwrap.dedent(text).lstrip("\n")


# (path, text, the line the plant is on)
PLANTS = {
    "cpp-line-comment": ("src/a.cpp", f"int a = 1;\nint b = 2;  // see {PLANT}\n", 2),
    "h-line-comment": ("include/a.h", f"#pragma once\n// see {PLANT}\n", 2),
    "hpp-block-comment": ("include/b.hpp", f"#pragma once\n/* first\n   {PLANT}\n */\n", 3),
    "py-comment": ("tools/plant.py", f"import os\n\nx = os.sep  # {PLANT}\n", 3),
    "py-module-docstring": ("src_python/m.py", f'"""Module.\n\nSee {PLANT}.\n"""\n', 3),
    "py-def-docstring": (
        "src_python/f.py",
        _dedent(f'''
            def f() -> None:
                """Do it.

                See {PLANT}.
                """
            '''),
        4,
    ),
    "py-async-def-docstring": (
        "src_python/g.py",
        f'async def g() -> None:\n    """See {PLANT}."""\n',
        2,
    ),
    "py-class-docstring": ("src_python/c.py", f'class C:\n    """See\n    {PLANT}."""\n', 3),
    "cmakelists-comment": ("CMakeLists.txt", f"project(x)\n# see {PLANT}\n", 2),
    "nested-cmakelists": ("src/CMakeLists.txt", f'add_library(y)  # "{PLANT}"\n', 1),
    "cmake-module": ("cmake/x.cmake", f"# see {PLANT}\n", 1),
    "root-md": ("NOTES.md", f"# Notes\n\nSee {PLANT}.\n", 3),
}


@pytest.mark.parametrize("plant", list(PLANTS))
def test_the_default_scan_finds_an_unpinned_governed_citation_in_every_scanned_kind(
    plantable: Scratch, plant: str
) -> None:
    rel, text, line = PLANTS[plant]
    plantable.write(rel, text)
    plantable.commit("plant")

    run = plantable.check()

    assert run.code == 1, run.out
    assert f"{rel}:{line}" in run.section("unpinned")


# Each file carries one real comment citation (the probe, which must be
# reported, so a scan that saw nothing cannot pass) and decoys that must not.
# (path, text, the probe's line)
DECOYS = {
    "cpp-string-literal": (
        "src/s.cpp",
        f'const char* s = "{PLANT}";\n// {PLANT}\n',
        2,
    ),
    "cpp-raw-string-spanning-lines": (
        "src/r.cpp",
        # A naive lexer closes the string at the inner quote and opens a block
        # comment that swallows the next line.
        f'const char* r = R"x(a " /* b\n{PLANT}\n)x"; // */\nint z = 0;  // {PLANT}\n',
        4,
    ),
    "cpp-string-after-digit-separator": (
        "src/d.cpp",
        f"int v = 0x9513'0000; const char* s = \"x'y // {PLANT}\";\n// {PLANT}\n",
        2,
    ),
    "cpp-escaped-quote-in-string": (
        "src/e.cpp",
        f'const char* s = "a \\" // {PLANT}";\n// {PLANT}\n',
        2,
    ),
    "py-plain-string": (
        "src/p.py",
        f'x = "{PLANT}"\n\n\ndef f() -> None:\n    y = 1\n    "{PLANT}"\n\n\n# {PLANT}\n',
        9,
    ),
    "cmake-quoted-argument": (
        "src/CMakeLists.txt",
        f'set(X "a \\" # {PLANT}")\n# {PLANT}\n',
        2,
    ),
    "md-url-and-absolute-path": (
        "src/u.md",
        f"https://host/{PLANT} and /abs/{PLANT}\n{PLANT}\n",
        2,
    ),
}


@pytest.mark.parametrize("case", list(DECOYS))
def test_a_citation_outside_a_comment_or_docstring_is_not_seen(
    plantable: Scratch, case: str
) -> None:
    rel, text, probe = DECOYS[case]
    plantable.write(rel, text)
    plantable.commit("decoys")

    run = plantable.check("src")

    assert run.code == 1, run.out
    assert run.count("unpinned") == 1, run.out
    assert f"{rel}:{probe}" in run.section("unpinned")


# A comment that a naive lexer would lose, each on line 1 of `src/k.cpp`.
HIDDEN_COMMENTS = {
    "after-string-with-slashes": f'const char* u = "http://x"; // {PLANT}\n',
    "after-string-with-block-opener": f'const char* s = "/*"; // {PLANT}\n',
    "after-char-literal-quote": f"char q = '\"'; // {PLANT}\n",
    "after-digit-separator-and-char": f"int v = 0x9513'0000; char c = 'a'; // {PLANT}\n",
    "block-after-string": f'const char* s = "*/"; /* {PLANT} */\n',
}


@pytest.mark.parametrize("case", list(HIDDEN_COMMENTS))
def test_a_comment_after_a_tricky_literal_is_still_seen(plantable: Scratch, case: str) -> None:
    plantable.write("src/k.cpp", HIDDEN_COMMENTS[case])
    plantable.commit("comment")

    run = plantable.check("src")

    assert run.code == 1, run.out
    assert "src/k.cpp:1" in run.section("unpinned")


@pytest.mark.parametrize(
    ("comment", "code"),
    [
        ("// see nowhere.md:1", 1),
        ("// `testing.md`, *Frameworks* says so", 0),
        ("// `testing.md`, *Nowhere* says so", 1),
        ("// `testing.md` §9", 1),
    ],
    ids=["missing-file", "heading-name-resolves", "heading-name-missing", "id-missing"],
)
def test_every_form_is_checked_inside_a_cpp_comment(
    plantable: Scratch, comment: str, code: int
) -> None:
    plantable.write("src/f.cpp", f"int f();\n{comment}\n")
    plantable.commit("comment")

    run = plantable.check("src")

    assert run.code == code, run.out
    if code:
        assert "src/f.cpp:2" in run.section("broken")


def test_vendored_legacy_worktree_and_ignored_files_are_not_scanned(plantable: Scratch) -> None:
    plantable.write(".gitignore", "scratch/\n")
    for rel in ("lib/v/a.hpp", "legacy/l.py", ".claude/worktrees/w/x.md"):
        plantable.write(rel, f"// {PLANT}\n" if rel.endswith("hpp") else f"# {PLANT}\n")
    plantable.write("docs/gone.md", f"{PLANT}\n")
    plantable.commit("excluded plants")
    (plantable.root / "docs" / "gone.md").unlink()  # tracked, but no longer on disk
    plantable.write("scratch/ignored.md", f"{PLANT}\n")
    plantable.write("docs/untracked.md", f"# New\n\n{PLANT}\n")

    run = plantable.check()

    assert run.code == 1, run.out
    assert "docs/untracked.md:3" in run.section("unpinned")
    for rel in ("lib/v/a.hpp", "legacy/l.py", ".claude/worktrees/w/x.md", "scratch/ignored.md"):
        assert rel not in run.out
    assert "Traceback" not in run.out


def test_a_python_file_that_does_not_parse_is_a_warning_not_a_failure(repo: Scratch) -> None:
    repo.write("src/bad.py", "def broken(:\n    pass\n")
    repo.commit("unparseable")

    run = repo.check("src")

    assert run.code == 0, run.out
    assert "src/bad.py" in run.section("warnings")
    assert "Traceback" not in run.out
