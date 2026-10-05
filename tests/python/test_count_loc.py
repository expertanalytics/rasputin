"""h16 T1: `tools/count_loc.py`, the one counter for `CLAUDE.md` §2.

Spec: `docs/increments/h16-harness-fixes.md` §2, T1. The tests, in the
design's order:

- `kind_of` and `counted_lines`, per kind (Python by `tokenize` and `ast`,
  C++ by a small scanner, CMake and shell by `#`);
- `hunks` on hand-written `git diff -U0` text;
- `main` end to end on a temporary repository whose commits give known
  counts, run as `python3 tools/count_loc.py <base> [<head>]` from a copy of
  the script committed in the repository's root commit, so the copy is never
  part of the diff it counts;
- that user git configuration and replace refs do not change the output;
- the two recorded counts from real history the design names.

The tool is loaded lazily through `harness_fixtures.Tool`, so while it is
missing every test that touches it fails naming `tools/count_loc.py`.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

from harness_fixtures import REAL, Tool, clean_env, git, is_work_tree_top, load_tool

cl = Tool("count_loc")

SCRIPT = "tools/count_loc.py"


def lines_of(text: str, kind_path: str) -> set[int]:
    """`counted_lines` for `text`, of the kind `kind_of(kind_path)` names."""
    kind = cl.kind_of(kind_path)
    assert kind is not None, f"kind_of({kind_path!r}) is None"
    return set(cl.counted_lines(text, kind))


# ---------------------------------------------------------------- kind_of

COUNTED_KINDS = {
    "src_python/tin_engine/gauge.py": "python",
    "src_python/tin_engine/_core.pyi": "python",  # a stub is Python (review round 1)
    "include/terrain/mesh/x.h": "c++",
    "include/terrain/mesh/x.hpp": "c++",
    "src/x.cpp": "c++",
    "src/x.cc": "c++",
    "src/x.cxx": "c++",
    "CMakeLists.txt": "cmake",
    "src/CMakeLists.txt": "cmake",
    "cmake/warnings.cmake": "cmake",
    "tools/run.sh": "shell",
}

NOT_CODE = (
    "README.md",
    "docs/increments/h16-harness-fixes.md",
    ".github/workflows/main.yaml",
    "pyproject.toml",
    "tests/python/data/x.json",
    "notes.txt",
    "LICENSE",
    "data/dem.tif",
)


@pytest.mark.parametrize("path", sorted(COUNTED_KINDS))
def test_kind_of_names_each_counted_kind(path: str) -> None:
    assert cl.kind_of(path) is not None, f"{path} should be counted"


@pytest.mark.parametrize("path", NOT_CODE)
def test_kind_of_is_none_for_files_of_no_kind(path: str) -> None:
    assert cl.kind_of(path) is None


def test_kind_of_groups_by_language() -> None:
    """Files of one language share a kind; Python and C++ do not."""
    by_language: dict[str, set[object]] = {}
    for path, language in COUNTED_KINDS.items():
        by_language.setdefault(language, set()).add(cl.kind_of(path))
    for language, kinds in by_language.items():
        assert len(kinds) == 1, f"{language} files map to {kinds}"
    assert cl.kind_of("a.py") != cl.kind_of("a.cpp")


# ---------------------------------------------------------------- Python

PYTHON_CASES = {
    "blank and comment-only lines": (
        "x = 1\n\n# a comment\n    # an indented comment\ny = 2  # trailing\n",
        {1, 5},
    ),
    "a one-line module docstring": ('"""Module."""\n\nimport os\n', {3}),
    "a three-line module docstring": ('"""Module.\n\nMore.\n"""\nimport os\n', {5}),
    "a one-line function docstring": (
        'def f():\n    """One line."""\n    return 1\n',
        {1, 3},
    ),
    "a three-line function docstring": (
        'def f():\n    """Three\n    lines\n    """\n    return 1\n',
        {1, 5},
    ),
    "a class docstring and a method docstring": (
        'class C:\n    """C.\n    """\n\n    def m(self):\n        """M."""\n        return 2\n',
        {1, 5, 7},
    ),
    "an async function docstring": (
        'async def f():\n    """F.\n    """\n    return 1\n',
        {1, 4},
    ),
    "a multi-line string that is not a docstring counts on its first line": (
        'X = """a\nb\nc"""\nY = 1\n',
        {1, 4},
    ),
    "a second string statement in a function is not a docstring": (
        'def f():\n    """Doc."""\n    """not a\n    docstring"""\n    return 1\n',
        {1, 3, 5},
    ),
    "a multi-line f-string counts on its first line only": (
        'X = f"""a\n{b}\nc"""\n',
        {1},
    ),
    "a bracketed expression counts on every line with a token": (
        "f(\n    1,\n\n    # c\n    2,\n)\n",
        {1, 2, 5, 6},
    ),
}


@pytest.mark.parametrize("name", sorted(PYTHON_CASES))
def test_python_counted_lines(name: str) -> None:
    text, expected = PYTHON_CASES[name]
    assert lines_of(text, "x.py") == expected


def test_python_lines_are_one_based_and_a_frozenset() -> None:
    result = cl.counted_lines("x = 1\n", cl.kind_of("x.py"))
    assert isinstance(result, frozenset)
    assert result == frozenset({1})


def test_python_that_does_not_tokenize_counts_every_non_blank_line() -> None:
    text = 'x = """never closed\n\ny = 1\n'
    assert lines_of(text, "x.py") == {1, 3}


# ---------------------------------------------------------------- C++

CPP_CASES = {
    "blank and // comment-only lines": ("int a;\n\n// c\n   // c\nint b; // c\n", {1, 5}),
    "a raw string over three lines with code after it": (
        'auto s = R"x(first\nmiddle\nlast)x";\n',
        {1, 3},
    ),
    'a raw string whose body holds a )" that is not its end': (
        'auto s = R"d(a)"\n)" still body\n)d";\n',
        {1, 3},
    ),
    "a raw string whose last line holds nothing after its end": (
        'auto s = R"d(\nbody\n)d"\n;\n',
        {1, 4},
    ),
    "/* */ across lines with code after it": (
        "/* a\n   b */ int x = 1;\n/* c\n*/\nint y;\n",
        {2, 5},
    ),
    "a one-line /* */ comment": ("/* c */\nint a; /* c */\n", {2}),
    "comment markers inside an ordinary string literal": (
        'const char* a = "/*";\nint b = 0;\nconst char* u = "//";\nint c;\n',
        {1, 2, 3, 4},
    ),
    # Review round 1: a planted fault that dropped character-literal handling
    # survived. Without it, the `"` inside `'"'` opens a string that ends at the
    # next `"`, and the `/*` after it is taken for a comment swallowing line 2.
    "a character literal holding a double quote, then /* inside a string": (
        'char q = \'"\'; const char* s = "/*";\nint x;\n',
        {1, 2},
    ),
}


@pytest.mark.parametrize("name", sorted(CPP_CASES))
def test_cpp_counted_lines(name: str) -> None:
    text, expected = CPP_CASES[name]
    assert lines_of(text, "x.cpp") == expected


# ---------------------------------------------------------------- CMake and shell


def test_cmake_blank_and_hash_comment_lines_do_not_count() -> None:
    text = "# c\n\nadd_library(x)\n  # c\ntarget_link_libraries(x y)\n"
    assert lines_of(text, "CMakeLists.txt") == {3, 5}


def test_shell_blank_and_hash_comment_lines_do_not_count() -> None:
    """A shebang is a line starting with `#`, so it does not count either."""
    text = "#!/bin/sh\n# c\n\necho hi\n  # c\nexit 0\n"
    assert lines_of(text, "run.sh") == {4, 6}


# ---------------------------------------------------------------- hunks

DIFF = """\
diff --git a/src/new.py b/src/new.py
new file mode 100644
index 0000000..1111111
--- /dev/null
+++ b/src/new.py
@@ -0,0 +1,3 @@
+a = 1
+
+b = 2
diff --git a/src/gone.py b/src/gone.py
deleted file mode 100644
index 2222222..0000000
--- a/src/gone.py
+++ /dev/null
@@ -1,2 +0,0 @@
-x = 1
-y = 2
diff --git a/src/mod.py b/src/mod.py
index 3333333..4444444 100644
--- a/src/mod.py
+++ b/src/mod.py
@@ -2 +2 @@ def f():
-    return 1
+    return 2
@@ -10,0 +11,2 @@ def g():
+    h()
+    i()
@@ -20,3 +22,0 @@ def k():
-p
-q
-r
"""


def test_hunks_reads_a_new_file_a_deletion_and_several_hunks() -> None:
    assert cl.hunks(DIFF) == {
        "src/new.py": ([1, 2, 3], []),
        "src/gone.py": ([], [1, 2]),
        "src/mod.py": ([2, 11, 12], [2, 20, 21, 22]),
    }


def test_hunks_does_not_take_a_changed_line_for_a_header() -> None:
    """`++i;` added and `--j;` removed show as `+++i;` and `---j;` in the diff."""
    diff = """\
diff --git a/src/x.cpp b/src/x.cpp
index 1111111..2222222 100644
--- a/src/x.cpp
+++ b/src/x.cpp
@@ -3,2 +3,2 @@ void f() {
---j;
-int k;
+++i;
+int k = 1;
"""
    assert cl.hunks(diff) == {"src/x.cpp": ([3, 4], [3, 4])}


def test_hunks_of_an_empty_diff_is_empty() -> None:
    assert cl.hunks("") == {}


# ---------------------------------------------------------------- main, end to end

A_BASE = '''\
"""Module doc."""

import os


def f():
    """One line."""
    return os.sep
'''

A_HEAD = '''\
"""Module doc."""

import os


def f():
    """One line."""
    return os.pathsep

# a comment
def g():
    """Three
    lines
    """
    return 1
'''

HPP_HEAD = """\
int a;
// c
/* a
   b */ int b;
auto s = R"d(
x
)d";
"""

#: The rows the fixture's branch must give: path -> (added, removed, net).
ROWS = {
    "src_python/pkg/a.py": (3, 1, 2),
    "src_python/pkg/gone.py": (0, 2, -2),
    "include/x.hpp": (3, 0, 3),
    "CMakeLists.txt": (1, 0, 1),
    "tools/run.sh": (1, 0, 1),
}
TOTAL = (8, 3, 5)
NOT_COUNTED = {
    "tests/python/test_x.py": "tests/",
    "docs/x.md": "docs/",
    "docs/tool.py": "docs/",
    "README.md": "not code",
    "src_python/pkg/data.json": "not code",
}


def commit_files(repo: Path, files: dict[str, str | None], message: str) -> None:
    """Write (or, for None, delete) each file and commit everything."""
    for relative, text in files.items():
        path = repo / relative
        if text is None:
            path.unlink()
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message)


def make_counted_repo(root: Path, base: dict[str, str]) -> Path:
    """A repository on `master` whose root commit holds the script's copy and `base`.

    The copy is committed before the base the counts start from, so it is
    never one of the changed files.
    """
    root.mkdir(parents=True)
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    source = REAL / SCRIPT
    if source.exists():
        (root / SCRIPT).parent.mkdir(parents=True)
        shutil.copy2(source, root / SCRIPT)
    commit_files(root, base, "base")
    return root


@pytest.fixture
def branch_repo(tmp_path: Path) -> Path:
    """`master` at the base, `feature` with the changes of ROWS and NOT_COUNTED,
    then one more commit on `master` that the count must not see."""
    repo = make_counted_repo(
        tmp_path.resolve() / "repo",
        {
            "src_python/pkg/a.py": A_BASE,
            "src_python/pkg/gone.py": "X = 1\n# c\nY = 2\n",
            "include/x.hpp": "int a;\n",
        },
    )
    git(repo, "checkout", "-q", "-b", "feature")
    commit_files(
        repo,
        {
            "src_python/pkg/a.py": A_HEAD,
            "src_python/pkg/gone.py": None,
            "include/x.hpp": HPP_HEAD,
            "CMakeLists.txt": "# c\n\nadd_library(x)\n",
            "tools/run.sh": "#!/bin/sh\n# c\necho hi\n",
            "tests/python/test_x.py": "def test_x():\n    assert True\n",
            "docs/x.md": "# x\n",
            "docs/tool.py": "x = 1\n",
            "README.md": "# r\n",
            "src_python/pkg/data.json": "{}\n",
        },
        "feature",
    )
    git(repo, "checkout", "-q", "master")
    commit_files(repo, {"src_python/pkg/later.py": "Z = 1\nW = 2\n"}, "master moves on")
    git(repo, "checkout", "-q", "feature")
    return repo


def run_count(
    repo: Path, *args: str, env: dict[str, str] | None = None, cwd: Path | None = None
) -> subprocess.CompletedProcess[str]:
    """The copy committed in `repo`, by path, from `cwd` (default `repo`)."""
    script = repo / SCRIPT
    if not script.exists():
        pytest.fail(f"{SCRIPT} is missing from the checkout")
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        cwd=cwd if cwd is not None else repo,
        env=env if env is not None else clean_env(),
        timeout=120,
        check=False,
    )


@dataclass(frozen=True)
class Report:
    rows: dict[str, tuple[int, int, int]]
    total: tuple[int, int, int]
    not_counted: dict[str, str]


def parse_report(result: subprocess.CompletedProcess[str]) -> Report:
    """Exit 0; tab-separated rows, one `total` row, then `not counted:` lines."""
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    rows: dict[str, tuple[int, int, int]] = {}
    totals: list[tuple[int, int, int]] = []
    skipped: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if line.startswith("not counted: "):
            rest = line.removeprefix("not counted: ")
            assert rest.endswith(")") and " (" in rest, f"not a `not counted:` line: {line!r}"
            path, _, reason = rest[:-1].rpartition(" (")
            skipped[path] = reason
            continue
        fields = line.split("\t")
        assert len(fields) == 4, f"not four tab-separated fields: {line!r}"
        numbers = (int(fields[1]), int(fields[2]), int(fields[3]))
        if fields[0] == "total":
            totals.append(numbers)
        else:
            rows[fields[0]] = numbers
    assert len(totals) == 1, f"{len(totals)} total lines in:\n{result.stdout}"
    return Report(rows, totals[0], skipped)


def nonzero(rows: dict[str, tuple[int, int, int]]) -> dict[str, tuple[int, int, int]]:
    """Rows with any line counted; whether a file with none gets a row is open."""
    return {path: row for path, row in rows.items() if row != (0, 0, 0)}


def test_main_counts_the_branch_against_its_merge_base(branch_repo: Path) -> None:
    report = parse_report(run_count(branch_repo, "master", "feature"))
    assert nonzero(report.rows) == ROWS
    assert report.total == TOTAL
    assert "src_python/pkg/later.py" not in report.rows, "master's own commit was counted"


def test_main_lists_each_skipped_file_with_its_reason(branch_repo: Path) -> None:
    report = parse_report(run_count(branch_repo, "master", "feature"))
    assert report.not_counted == NOT_COUNTED


def test_main_prints_rows_then_the_total_then_the_skipped(branch_repo: Path) -> None:
    lines = run_count(branch_repo, "master", "feature").stdout.splitlines()
    kinds = [
        "skip" if line.startswith("not counted: ") else
        "total" if line.startswith("total\t") else "row"
        for line in lines
    ]  # fmt: skip
    assert kinds == sorted(kinds, key=["row", "total", "skip"].index), lines


def test_main_head_defaults_to_head(branch_repo: Path) -> None:
    given = run_count(branch_repo, "master", "feature")
    default = run_count(branch_repo, "master")
    assert default.returncode == 0, default.stderr
    assert default.stdout == given.stdout


def test_main_reads_committed_revisions_not_the_working_tree(branch_repo: Path) -> None:
    before = run_count(branch_repo, "master")
    (branch_repo / "src_python/pkg/a.py").write_text(A_HEAD + "extra = 1\nmore = 2\n")
    (branch_repo / "src_python/pkg/untracked.py").write_text("u = 1\n")
    git(branch_repo, "add", "src_python/pkg/a.py")
    after = run_count(branch_repo, "master")
    assert after.returncode == 0, after.stderr
    assert after.stdout == before.stdout


def test_main_an_unknown_revision_is_exit_2_with_one_stderr_line(branch_repo: Path) -> None:
    result = run_count(branch_repo, "no-such-branch", "feature")
    assert result.returncode == 2, f"exit {result.returncode}; stdout: {result.stdout}"
    lines = [line for line in result.stderr.splitlines() if line.strip()]
    assert len(lines) == 1, result.stderr
    assert lines[0].startswith("count_loc:"), result.stderr


def test_main_a_pure_rename_counts_nothing(tmp_path: Path) -> None:
    body = "".join(f"v{n} = {n}\n" for n in range(20))
    repo = make_counted_repo(tmp_path.resolve() / "repo", {"src/a.py": body})
    git(repo, "mv", "src/a.py", "src/b.py")
    git(repo, "commit", "-q", "-m", "rename")
    report = parse_report(run_count(repo, "HEAD~1", "HEAD"))
    assert report.total == (0, 0, 0)


def test_main_a_python_file_that_does_not_tokenize_counts_non_blank_lines_and_warns(
    tmp_path: Path,
) -> None:
    repo = make_counted_repo(tmp_path.resolve() / "repo", {"src/keep.py": "k = 1\n"})
    commit_files(repo, {"src/bad.py": 'x = """never closed\n\ny = 1\n'}, "bad")
    result = run_count(repo, "HEAD~1", "HEAD")
    report = parse_report(result)
    assert report.rows["src/bad.py"] == (2, 0, 2)
    assert "src/bad.py" in result.stderr, f"no warning naming the file: {result.stderr!r}"


# ---------------------------------------------------------------- what git is told


def test_user_config_does_not_change_the_output(branch_repo: Path, tmp_path: Path) -> None:
    """`diff.noprefix` (and colour, and mnemonic prefixes) in a user config file
    change the headers the parser reads; the tool's git calls ignore that file."""
    plain = run_count(branch_repo, "master", "feature")
    config = tmp_path / "gitconfig"
    config.write_text(
        "[diff]\n\tnoprefix = true\n\tmnemonicPrefix = true\n"
        "[color]\n\tui = always\n\tdiff = always\n"
    )
    hostile = run_count(
        branch_repo, "master", "feature", env={**clean_env(), "GIT_CONFIG_GLOBAL": str(config)}
    )
    assert hostile.returncode == 0, hostile.stderr
    assert hostile.stdout == plain.stdout
    probe = subprocess.run(
        ["git", "-C", str(branch_repo), "diff", "-U0", "master...feature"],
        capture_output=True,
        text=True,
        env={**clean_env(), "GIT_CONFIG_GLOBAL": str(config)},
        check=True,
    )
    assert "+++ b/" not in probe.stdout, "the hostile config did not change git's headers"


def test_replace_refs_do_not_change_the_output(branch_repo: Path) -> None:
    plain = run_count(branch_repo, "master", "feature")
    planted = branch_repo / "planted.py"
    planted.write_text("".join(f"p{n} = {n}\n" for n in range(50)))
    fake = git(branch_repo, "hash-object", "-w", str(planted)).strip()
    planted.unlink()
    real = git(branch_repo, "rev-parse", "feature:src_python/pkg/a.py").strip()
    git(branch_repo, "replace", real, fake)
    shown = git(branch_repo, "show", "feature:src_python/pkg/a.py")
    assert "p49 = 49" in shown, "the replace ref did not take"
    replaced = run_count(branch_repo, "master", "feature")
    assert replaced.returncode == 0, replaced.stderr
    assert replaced.stdout == plain.stdout


def test_repository_diff_prefixes_do_not_change_the_output(branch_repo: Path) -> None:
    """Review round 1: a planted fault that dropped the explicit prefixes
    survived. `diff.srcPrefix`/`diff.dstPrefix` in the repository's own config,
    which `GIT_CONFIG_GLOBAL` does not reach, move the headers the parser reads."""
    plain = run_count(branch_repo, "master", "feature")
    git(branch_repo, "config", "diff.srcPrefix", "SRC/")
    git(branch_repo, "config", "diff.dstPrefix", "DST/")
    probe = git(branch_repo, "diff", "-U0", "master...feature")
    assert "+++ DST/" in probe, "the repository config did not change git's headers"
    hostile = run_count(branch_repo, "master", "feature")
    assert hostile.returncode == 0, hostile.stderr
    assert hostile.stdout == plain.stdout


#: `diff.interHunkContext=50`, three ways that `GIT_CONFIG_GLOBAL` does not
#: reach: the two environment forms and the repository's own config. Fused
#: hunks carry their context lines in the `@@` ranges, so they would count.
INTER_HUNK = {
    "GIT_CONFIG_COUNT": {
        "GIT_CONFIG_COUNT": "1",
        "GIT_CONFIG_KEY_0": "diff.interHunkContext",
        "GIT_CONFIG_VALUE_0": "50",
    },
    "GIT_CONFIG_PARAMETERS": {"GIT_CONFIG_PARAMETERS": "'diff.interHunkContext'='50'"},
    "repository config": {},
}


def inter_hunk_env(repo: Path, how: str) -> dict[str, str]:
    """The environment for `how`; for the repository config, set it in `repo`."""
    if how == "repository config":
        git(repo, "config", "diff.interHunkContext", "50")
    return {**clean_env(), **INTER_HUNK[how]}


@pytest.mark.parametrize("how", sorted(INTER_HUNK))
def test_inter_hunk_context_does_not_change_the_output(tmp_path: Path, how: str) -> None:
    """Two one-line changes with a line of code between them: two hunks, 2 added
    and 2 removed. Fused into one hunk, the middle line would count both ways."""
    repo = make_counted_repo(tmp_path.resolve() / "repo", {"src/m.py": "a = 1\nx = 3\nb = 2\n"})
    commit_files(repo, {"src/m.py": "a = 10\nx = 3\nb = 20\n"}, "two changes")
    plain = parse_report(run_count(repo, "HEAD~1", "HEAD"))
    assert plain.rows["src/m.py"] == (2, 2, 0)
    env = inter_hunk_env(repo, how)
    fused = subprocess.run(
        ["git", "-C", str(repo), "diff", "-U0", "HEAD~1", "HEAD"],
        capture_output=True, text=True, env=env, check=True,
    ).stdout  # fmt: skip
    assert fused.count("\n@@ ") == 1, f"{how} did not fuse the hunks:\n{fused}"
    hostile = parse_report(run_count(repo, "HEAD~1", "HEAD", env=env))
    assert hostile.rows == plain.rows
    assert hostile.total == plain.total


#: `diff.relative=true`, the same three ways as INTER_HUNK. Run from a
#: subdirectory, git would then diff only that subdirectory, with paths
#: relative to it. `--no-relative` (design T1) neutralises all three forms, so
#: this test fails only if both it and the stripping of `GIT_CONFIG_*` go; the
#: stripping alone is pinned by `test_env_drops_config_given_in_the_environment`.
RELATIVE = {
    "GIT_CONFIG_COUNT": {
        "GIT_CONFIG_COUNT": "1",
        "GIT_CONFIG_KEY_0": "diff.relative",
        "GIT_CONFIG_VALUE_0": "true",
    },
    "GIT_CONFIG_PARAMETERS": {"GIT_CONFIG_PARAMETERS": "'diff.relative'='true'"},
    "repository config": {},
}


def relative_env(repo: Path, how: str) -> dict[str, str]:
    """The environment for `how`; for the repository config, set it in `repo`."""
    if how == "repository config":
        git(repo, "config", "diff.relative", "true")
    return {**clean_env(), **RELATIVE[how]}


@pytest.mark.parametrize("how", sorted(RELATIVE))
def test_diff_relative_run_from_a_subdirectory_does_not_change_the_output(
    branch_repo: Path, how: str
) -> None:
    """From `src_python/`, with `diff.relative=true`, git sees only that
    subdirectory: the rows of `include/`, `tools/` and the top level would be
    lost and the remaining paths would lose their `src_python/` prefix."""
    subdir = branch_repo / "src_python"
    plain = run_count(branch_repo, "master", "feature")
    assert parse_report(plain).total == TOTAL
    from_subdir = run_count(branch_repo, "master", "feature", cwd=subdir)
    assert from_subdir.stdout == plain.stdout, "the subdirectory alone changed the output"
    env = relative_env(branch_repo, how)
    probe = subprocess.run(
        ["git", "diff", "--name-only", "master...feature"],
        capture_output=True, text=True, cwd=subdir, env=env, check=True,
    ).stdout  # fmt: skip
    assert "include/x.hpp" not in probe and "pkg/a.py" in probe.splitlines(), (
        f"{how} did not make git's diff relative to {subdir}:\n{probe}"
    )
    hostile = run_count(branch_repo, "master", "feature", env=env, cwd=subdir)
    assert hostile.returncode == 0, hostile.stderr
    assert hostile.stdout == plain.stdout


#: A change on which git's diff algorithms disagree (found by search): myers,
#: git's default, keeps `b = 2` twice and gives 3 added and 1 removed;
#: `histogram` and `patience` anchor on the unique `c = 3` and give 4 and 2.
ALGORITHM_BASE = "b = 2\nb = 2\nc = 3\n"
ALGORITHM_HEAD = "c = 3\nb = 2\nb = 2\nd = 4\na = 1\n"


@pytest.mark.parametrize("algorithm", ["histogram", "patience"])
def test_repository_diff_algorithm_does_not_change_the_output(
    tmp_path: Path, algorithm: str
) -> None:
    """Review round 3, Ola's ruling: the counter pins myers. `diff.algorithm`
    in the repository's own config, which `GIT_CONFIG_GLOBAL` and the
    environment stripping do not reach, must not change the counts."""
    repo = make_counted_repo(tmp_path.resolve() / "repo", {"src/m.py": ALGORITHM_BASE})
    commit_files(repo, {"src/m.py": ALGORITHM_HEAD}, "reorder")
    plain = parse_report(run_count(repo, "HEAD~1", "HEAD"))
    assert plain.rows["src/m.py"] == (3, 1, 2)
    git(repo, "config", "diff.algorithm", algorithm)
    probe = git(repo, "diff", "--numstat", "HEAD~1", "HEAD")
    assert probe.split()[:2] == ["4", "2"], f"{algorithm} did not change git's diff:\n{probe}"
    hostile = parse_report(run_count(repo, "HEAD~1", "HEAD"))
    assert hostile.rows == plain.rows
    assert hostile.total == plain.total


#: Paths the counter does not count, one per reason: `docs/`, `tests/`, a
#: file of no kind.
UNCOUNTED = ("docs/proto.py", "tests/proto.py", "tools/proto.txt")
COUNTED = "tools/proto.py"
PROTO = "".join(f"v{n} = {n}\n" for n in range(40))


def moved_repo(root: Path, source: str, target: str) -> Path:
    """`source` holds PROTO's 40 code lines; the head commit moves it to
    `target` and adds one line, and git reports the move as a rename."""
    repo = make_counted_repo(root, {source: PROTO})
    # `git mv` does not create directories, and the base tree may lack
    # `target`'s parent (`docs/`, `tests/`).
    (repo / target).parent.mkdir(parents=True, exist_ok=True)
    git(repo, "mv", source, target)
    commit_files(repo, {target: PROTO + "extra = 1\n"}, "move")
    status = git(repo, "diff", "-M", "--name-status", "HEAD~1", "HEAD")
    assert status.startswith("R"), f"git did not see a rename:\n{status}"
    return repo


@pytest.mark.parametrize("source", UNCOUNTED)
def test_a_rename_into_a_counted_path_counts_as_an_added_file(tmp_path: Path, source: str) -> None:
    """Review round 1, Ola's ruling: a file moved from a path the counter does
    not count into one it counts counts in full, 41 lines, as a new file would."""
    repo = moved_repo(tmp_path.resolve() / "repo", source, COUNTED)
    report = parse_report(run_count(repo, "HEAD~1", "HEAD"))
    assert report.total == (41, 0, 41)


@pytest.mark.parametrize("target", UNCOUNTED)
def test_a_rename_out_of_a_counted_path_counts_as_a_removed_file(
    tmp_path: Path, target: str
) -> None:
    """The other direction of the same ruling: the 40 lines leave the counted
    tree, so they count in full as removed, as a deletion would."""
    repo = moved_repo(tmp_path.resolve() / "repo", COUNTED, target)
    report = parse_report(run_count(repo, "HEAD~1", "HEAD"))
    assert report.total == (0, 40, -40)


# ---------------------------------------------------------------- recorded counts

#: (base, head) -> the rows with lines counted, and the total, as recorded:
#: 29 PR 2 review round 3 (569 net) and 29 PR 4 (691 net, ef403bd's message
#: and the round-6 review). The per-file rows follow from §2's rule; the
#: totals are the record.
RECORDED = {
    ("529613a", "193079d"): (
        {
            "src_python/tin_engine/burn.py": (175, 0, 175),
            "src_python/tin_engine/catchment.py": (110, 13, 97),
            "src_python/tin_engine/cli.py": (102, 7, 95),
            "src_python/tin_engine/gauge.py": (117, 0, 117),
            "src_python/tin_engine/sensitivity.py": (85, 0, 85),
        },
        (589, 20, 569),
    ),
    ("9e666f4", "9bb1723"): (
        {
            "src_python/tin_engine/catchment.py": (13, 0, 13),
            "src_python/tin_engine/catchment_batch.py": (195, 0, 195),
            "src_python/tin_engine/cli.py": (241, 40, 201),
            "src_python/tin_engine/fetch/nve.py": (19, 9, 10),
            "src_python/tin_engine/gauge.py": (28, 1, 27),
            "src_python/tin_engine/io/geojson.py": (21, 0, 21),
            "src_python/tin_engine/io/station_set.py": (34, 2, 32),
            "src_python/tin_engine/mosaic.py": (6, 1, 5),
            "src_python/tin_engine/reference.py": (187, 0, 187),
        },
        (744, 53, 691),
    ),
}


def needs_history(*revs: str) -> None:
    """Skip, with the reason, outside a git work tree (a `git archive` copy, as
    `tools/scratch_copy.py` makes); fail in a work tree that lacks a commit (a
    shallow clone), so CI's full clone (`fetch-depth: 0`) cannot be lost silently."""
    if not is_work_tree_top(REAL):
        pytest.skip(f"{REAL} is not the top of a git work tree, so it has no history to count")
    for rev in revs:
        found = subprocess.run(
            ["git", "-C", str(REAL), "cat-file", "-e", f"{rev}^{{commit}}"],
            capture_output=True,
            env=clean_env(),
            check=False,
        )
        assert found.returncode == 0, f"{rev} is not in this clone; it needs the full history"


@pytest.mark.parametrize(("base", "head"), sorted(RECORDED))
def test_recorded_counts_from_real_history(base: str, head: str) -> None:
    """Needs the full history; see `needs_history` for when it skips and when it fails."""
    needs_history(base, head)
    rows, total = RECORDED[(base, head)]
    report = parse_report(run_count(REAL, base, head))
    assert nonzero(report.rows) == rows
    assert report.total == total


#: PR #163 (`45acf22`, a merge; its first parent is the base), counted by hand
#: in review round 1: `refine.hpp` 10, 6, 4 and `lattice_mesh.hpp` 31, 14, 17.
PR_163 = ("45acf22^1", "45acf22")
PR_163_ROWS = {
    "include/terrain/mesh/lattice_mesh.hpp": (31, 14, 17),
    "include/terrain/refinement/refine.hpp": (10, 6, 4),
}


@pytest.mark.parametrize("how", ["none", *sorted(set(INTER_HUNK) - {"repository config"})])
def test_pr_163_counts_whatever_the_environment_says_of_hunk_context(how: str) -> None:
    """Review round 1: with `diff.interHunkContext=50` in the environment, git
    fuses PR #163's hunks and the count became 63, 42, 21 instead of 41, 20, 21."""
    needs_history(*PR_163)
    env = clean_env() if how == "none" else {**clean_env(), **INTER_HUNK[how]}
    report = parse_report(run_count(REAL, *PR_163, env=env))
    assert nonzero(report.rows) == PR_163_ROWS
    assert report.total == (41, 20, 21)


#: Config given in the environment, in both of git's forms, with two numbered
#: keys so a pattern that matched only `_0` would be caught.
CONFIG_IN_ENV = {
    "GIT_CONFIG_COUNT": "2",
    "GIT_CONFIG_KEY_0": "diff.algorithm",
    "GIT_CONFIG_VALUE_0": "histogram",
    "GIT_CONFIG_KEY_12": "diff.relative",
    "GIT_CONFIG_VALUE_12": "true",
    "GIT_CONFIG_PARAMETERS": "'diff.algorithm'='patience'",
}


def test_env_drops_config_given_in_the_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    """Review round 3: every flag `DIFF` pins also overrides the environment
    forms end to end, so only this test fails if `_env()` stops stripping them.
    A key no flag pins would get through."""
    for key, value in CONFIG_IN_ENV.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setenv("RASPUTIN_UNRELATED", "kept")
    monkeypatch.setenv("GIT_CONFIG_KEYS", "not git's")  # near miss of the pattern
    env = load_tool("count_loc")._env()  # `Tool` refuses private names
    assert sorted(set(CONFIG_IN_ENV) & set(env)) == []
    assert env["RASPUTIN_UNRELATED"] == "kept"
    assert env["GIT_CONFIG_KEYS"] == "not git's"
    assert env["GIT_CONFIG_GLOBAL"] == os.devnull
