"""`tools/check_code_history.py`: no project history in docstrings, comments or help text.

The rule (`docs/increments/42-readable-docstrings.md`, section 4): docstrings,
comments and the CLI's `help=` strings describe the code as it is, never its
history. The checker reports each match of seven rules (`increment`, `label`,
`ruling-id`, `review`, `person`, `hash`, `pointer`) as `path:line: rule: text`
and exits 1 when it reports anything. Its scope is
`[tool.check_code_history] paths` in `pyproject.toml`, overridden by paths on
the command line. `--same-code BASE` checks that every `.py` file that differs
from `BASE` differs only in docstrings and comment prose: the same AST once
bare strings are removed, and the same tool directives in its comments.

Each test writes its module inline into a scratch repository, copies the tool
in, and runs it there as a subprocess: the tests see only its exit status and
its output, never its internals.
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

# A harness test: CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

TOOL = Path(__file__).resolve().parents[2] / "tools" / "check_code_history.py"

REPORT = re.compile(r"^(?P<path>\S+?):(?P<line>\d+): (?P<rule>[a-z-]+): (?P<text>.+)$")
RULES = ("increment", "label", "ruling-id", "review", "person", "hash", "pointer")


@dataclass(frozen=True)
class Run:
    code: int
    out: str
    err: str

    @property
    def reports(self) -> list[tuple[str, int, str]]:
        """Every report line on stdout, as (path, line, rule)."""
        found = []
        for text in self.out.splitlines():
            match = REPORT.match(text)
            if match:
                found.append((match["path"], int(match["line"]), match["rule"]))
        return found

    @property
    def located(self) -> set[tuple[int, str]]:
        return {(line, rule) for _, line, rule in self.reports}


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _pyproject(repo: Path, paths: list[str]) -> None:
    listed = ", ".join(f'"{p}"' for p in paths)
    (repo / "pyproject.toml").write_text(f"[tool.check_code_history]\npaths = [{listed}]\n")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A git repository with an empty scope and a `pkg/` directory, no commit yet."""
    root = tmp_path / "repo"
    (root / "tools").mkdir(parents=True)
    (root / "pkg").mkdir()
    _git(root, "init", "-q", "-b", "master")
    _pyproject(root, [])
    return root


def _write(repo: Path, name: str, source: str) -> str:
    """Write `source` (dedented) to `pkg/<name>`; return its repository path."""
    (repo / "pkg" / name).write_text(textwrap.dedent(source).lstrip("\n"))
    return f"pkg/{name}"


def _check(repo: Path, *args: str) -> Run:
    assert TOOL.is_file(), f"{TOOL} does not exist yet"
    shutil.copy(TOOL, repo / "tools" / "check_code_history.py")
    result = subprocess.run(
        [sys.executable, "tools/check_code_history.py", *args],
        cwd=repo,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return Run(result.returncode, result.stdout, result.stderr)


def _line(repo: Path, path: str, marker: str) -> int:
    """The 1-based number of the one line of `path` holding `marker`."""
    lines = (repo / path).read_text().splitlines()
    hits = [n for n, text in enumerate(lines, start=1) if marker in text]
    assert len(hits) == 1, f"{marker!r} is on lines {hits} of the fixture"
    return hits[0]


# -- the rules ------------------------------------------------------------------

# One instance of each rule, each on a line no other rule matches.
INSTANCES = {
    "increment": "Added in increment 12.",
    "label": "Kept as in 16c for now.",
    "ruling-id": "Follows R5.",
    "review": "Settled in the review round.",
    "person": "Asked for by Ola.",
    "hash": "Introduced in `abc1234`.",
    "pointer": "Derived in notes.md.",
}

DIRTY = (
    '"""A module with history in one function and in comments."""\n\n\n'
    "def f() -> None:\n"
    '    """Do nothing.\n\n'
    + "".join(f"    D {text}\n" for text in INSTANCES.values())
    + '    """\n\n\n'
    + "".join(f"# C {text}\n" for text in INSTANCES.values())
)


def test_each_rule_fires_in_a_function_docstring_and_in_a_comment(repo: Path) -> None:
    path = _write(repo, "dirty.py", DIRTY)

    run = _check(repo, path)

    expected = {
        (_line(repo, path, f"{where} {text}"), rule)
        for rule, text in INSTANCES.items()
        for where in ("D", "C")
    }
    assert run.located == expected, run.out + run.err
    assert len(run.reports) == 14, run.out


def test_class_and_async_function_docstrings_are_read(repo: Path) -> None:
    path = _write(
        repo,
        "kinds.py",
        '''
        class Tile:
            """A tile, as in increment 12."""


        async def fetch() -> None:
            """Fetch a tile, as in increment 13."""
        ''',
    )

    run = _check(repo, path)

    assert run.located == {
        (_line(repo, path, "increment 12"), "increment"),
        (_line(repo, path, "increment 13"), "increment"),
    }, run.out + run.err


def test_string_literals_in_code_are_not_prose(repo: Path) -> None:
    path = _write(
        repo,
        "code.py",
        '''
        """Clean."""

        X = "R5 16c"
        print("increment 12")
        x = 1
        y = f"{x} Ola"
        Z = 2
        """A bare string after an assignment, naming increment 12, is no docstring."""
        ''',
    )

    run = _check(repo, path)

    assert run.code == 0, run.out + run.err
    assert run.reports == []


def test_help_strings_are_prose_reported_where_the_value_starts(repo: Path) -> None:
    path = _write(
        repo,
        "cli.py",
        '''
        """Clean."""

        import typer

        app = typer.Typer()
        x = "v"
        option = typer.Option(help="increment 12")


        @app.command(
            help="Serve." f"{x}"
            " (increment 37)"
        )
        def serve() -> None:
            pass
        ''',
    )

    run = _check(repo, path)

    assert run.located == {
        (_line(repo, path, "option = typer.Option"), "increment"),
        (_line(repo, path, 'help="Serve."'), "increment"),
    }, run.out + run.err


@pytest.mark.parametrize("keyword", ["short_help", "epilog"])
def test_short_help_and_epilog_are_read_like_help(repo: Path, keyword: str) -> None:
    path = _write(
        repo,
        "cli.py",
        f'''
        """Clean."""

        import typer

        app = typer.Typer()


        @app.command({keyword}="Draw a map (increment 29).")
        def draw() -> None:
            pass
        ''',
    )

    run = _check(repo, path)

    assert run.located == {(_line(repo, path, "@app.command"), "increment")}, run.out + run.err


# -- what the rules leave alone -------------------------------------------------


def test_scientific_notation_format_specs_paths_and_d8_are_not_history(repo: Path) -> None:
    # The `.md` path is in the one place a pointer is allowed (the module
    # docstring, heading form), so that only the `label` rule is under test
    # for `15c` inside it; the pointer rule is tested below.
    path = _write(
        repo,
        "exclusions.py",
        '''
        """Tolerances.

        `docs/increments/15c-geographic-dem.md`, *3 The geographic grid*
        """


        def agree() -> None:
            """Agree to within 1e-6 of a metre, or 1e-4 in degrees; D8 flow directions."""


        # Agree to within 1e-6 of a metre, or 1e-4 in degrees.
        # Printed as f"{x:.0f}" metres, or f"{x:.2f}" with decimals.
        # D8 flow directions point downhill.
        ''',
    )

    run = _check(repo, path)

    assert run.code == 0, run.out + run.err
    assert run.reports == []


@pytest.mark.parametrize(
    "comment",
    ["# Kept as in 16c.", "# Kept as before 15c-2.", "# Kept as in 16c. Then more."],
)
def test_a_label_that_ends_a_sentence_is_a_label(repo: Path, comment: str) -> None:
    # A full stop after a label ends the sentence; it does not hide the label.
    path = _write(repo, "full_stop.py", f"X = 1\n{comment}\n")

    run = _check(repo, path)

    assert run.code == 1, run.out + run.err
    assert run.located == {(2, "label")}, run.out + run.err


def test_a_label_that_ends_a_docstring_sentence_is_a_label(repo: Path) -> None:
    path = _write(
        repo,
        "full_stop_doc.py",
        '''
        def f() -> None:
            """Do nothing.

            Kept as in 16c.
            """
        ''',
    )

    run = _check(repo, path)

    assert run.located == {(_line(repo, path, "16c."), "label")}, run.out + run.err


@pytest.mark.parametrize("comment", ["# Read 16c.txt first.", "# Scaled by 3f.5 here."])
def test_a_full_stop_before_a_word_character_still_hides_a_label(repo: Path, comment: str) -> None:
    # A file suffix (`16c.txt`) or a decimal (`3f.5`) is not a sentence end.
    path = _write(repo, "suffix.py", f"X = 1\n{comment}\n")

    run = _check(repo, path)

    assert run.code == 0, run.out + run.err
    assert run.reports == []


def test_a_label_inside_a_path_is_not_a_label_but_the_path_is_a_pointer(repo: Path) -> None:
    path = _write(
        repo,
        "path_comment.py",
        """
        # See `docs/increments/15c-geographic-dem.md`.
        X = 1
        """,
    )

    run = _check(repo, path)

    assert run.located == {(1, "pointer")}, run.out + run.err


def test_tool_directives_are_not_history(repo: Path) -> None:
    path = _write(
        repo,
        "directives.py",
        """
        # ruff: noqa: RUF001
        # mypy: disable-error-code="attr-defined"
        import os  # noqa: E501, C901, B008, F401
        X = os.sep  # type: ignore[attr-defined]
        Y = 1  # fmt: skip
        Z = 2  # pragma: no cover
        W = 3  # NOQA: D205
        """,
    )

    run = _check(repo, path)

    assert run.code == 0, run.out + run.err
    assert run.reports == []


def test_a_link_to_an_outside_bug_report_with_its_reason_passes(repo: Path) -> None:
    # Section 4 allows a link to an outside bug report where the code works
    # around that bug, with the reason beside it. None of the seven rules
    # matches such a link, issue numbers, pull requests and commits included.
    path = _write(
        repo,
        "workaround.py",
        """
        def peak(values: list[float]) -> float:
            # NumPy warns on an all-NaN slice and returns NaN
            # (https://github.com/numpy/numpy/issues/21004), so the default is
            # passed in rather than taken from the result.
            # PROJ returns inf for a point outside the area of use
            # (https://github.com/OSGeo/PROJ/issues/2345#issuecomment-1234567890,
            # fixed upstream by https://github.com/OSGeo/PROJ/pull/3512 and
            # https://github.com/OSGeo/PROJ/commit/3f2a9c1d; see numpy/numpy#21004
            # and https://github.com/OSGeo/PROJ/blob/master/NEWS.md), so inf is
            # caught here.
            return max(values)
        """,
    )

    run = _check(repo, path)

    assert run.code == 0, run.out + run.err
    assert run.reports == []


# -- pointers -------------------------------------------------------------------

ALLOWED = "`docs/increments/32-landcover-simplify.md`, *5.3 The band in plan, not over the DEM*"


def test_one_heading_form_pointer_in_the_module_docstring_is_allowed(repo: Path) -> None:
    path = _write(repo, "one.py", f'"""A module.\n\n{ALLOWED}\n"""\n')

    run = _check(repo, path)

    assert run.code == 0, run.out + run.err
    assert run.reports == []


def test_every_other_pointer_is_reported(repo: Path) -> None:
    source = (
        '"""A module with pointers.\n\n'
        f"{ALLOWED}\n"
        "`docs/increments/15-dem-mosaic.md`, *1 A second heading*\n"
        "See `x.md` section 3.\n"
        '"""\n\n'
        "import typer\n\n\n"
        "def f() -> None:\n"
        '    """See `docs/y.md` for the derivation."""\n\n\n'
        "# Derived in `docs/z.md`.\n"
        'option = typer.Option(help="See `docs/h.md`.")\n'
    )
    path = _write(repo, "pointers.py", source)

    run = _check(repo, path)

    markers = ("15-dem-mosaic.md", "`x.md`", "docs/y.md", "docs/z.md", "docs/h.md")
    assert run.located == {(_line(repo, path, m), "pointer") for m in markers}, run.out + run.err


# -- exit status and output -----------------------------------------------------


def test_a_clean_file_exits_0_with_no_output(repo: Path) -> None:
    path = _write(repo, "clean.py", '"""A clean module."""\n\n# A clean comment.\nX = 1\n')

    run = _check(repo, path)

    assert (run.code, run.out, run.err) == (0, "", "")


def test_a_dirty_file_exits_1_with_one_report_line_per_match(repo: Path) -> None:
    path = _write(repo, "dirty.py", DIRTY)

    run = _check(repo, path)

    assert run.code == 1, run.err
    lines = run.out.splitlines()
    assert len(lines) == 14, run.out
    for text in lines:
        match = REPORT.match(text)
        assert match, f"not of the form path:line: rule: text: {text!r}"
        assert Path(match["path"]).name == "dirty.py"
        assert match["rule"] in RULES
        assert match["text"].strip()


# -- scope ----------------------------------------------------------------------


def test_an_empty_scope_checks_nothing(repo: Path) -> None:
    _write(repo, "dirty.py", DIRTY)

    run = _check(repo)

    assert run.code == 0, run.out + run.err
    assert run.reports == []


def test_a_listed_directory_is_checked(repo: Path) -> None:
    _write(repo, "dirty.py", DIRTY)
    _pyproject(repo, ["pkg"])

    run = _check(repo)

    assert run.code == 1, run.out + run.err
    assert len(run.reports) == 14, run.out
    assert {Path(p).name for p, _, _ in run.reports} == {"dirty.py"}


def test_a_path_on_the_command_line_overrides_the_scope(repo: Path) -> None:
    _write(repo, "dirty.py", DIRTY)

    run = _check(repo, "pkg")

    assert run.code == 1, run.out + run.err
    assert len(run.reports) == 14, run.out


def test_a_listed_path_that_does_not_exist_is_an_error_naming_it(repo: Path) -> None:
    _pyproject(repo, ["no_such_dir"])

    run = _check(repo)

    assert run.code != 0
    assert "no_such_dir" in run.out + run.err
    assert run.reports == []


# -- same-code mode -------------------------------------------------------------

MODULE = '''\
"""A small module."""

# ruff: noqa: RUF001
from __future__ import annotations

LIMIT = 10  # NOQA: E501


class Missing(Exception):  # noqa: N818 a name kept for callers
    """Raised when a tile is missing."""


def load(path: str) -> bytes:
    """Read a file."""
    # Open read-only.
    with open(path, "rb") as handle:
        return handle.read()


def size(data: bytes) -> int:
    """Count bytes."""
    return len(data)  # type: ignore[arg-type]
'''


def _replace(old: str, new: str) -> Callable[[str], str]:
    def edit(text: str) -> str:
        assert old in text, f"{old!r} is not in the fixture"
        return text.replace(old, new, 1)

    return edit


def _prepend(line: str) -> Callable[[str], str]:
    return lambda text: f"{line}\n{text}"


def _committed(repo: Path, source: str) -> str:
    """Commit `source` as `pkg/mod.py`; return the commit."""
    (repo / "pkg" / "mod.py").write_text(source)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "base")
    return _git(repo, "rev-parse", "HEAD")


SAME = {
    "nothing changed": lambda text: text,
    "docstring rewritten": _replace('"""Read a file."""', '"""The bytes of ``path``."""'),
    "docstring reflowed": _replace(
        '"""Read a file."""', '"""Read\n    a file, all of it.\n\n    Opened read-only.\n    """'
    ),
    "module docstring rewritten": _replace('"""A small module."""', '"""Tiles and sizes."""'),
    "comment rewritten": _replace("# Open read-only.", "# The file is only ever read."),
    "noqa reason reworded": _replace("a name kept for callers", "callers import this name"),
}

DIFFERS = {
    "string literal changed": _replace('"rb"', '"r+b"'),
    "noqa dropped": _replace("  # noqa: N818 a name kept for callers", ""),
    "upper-case NOQA dropped": _replace("  # NOQA: E501", ""),
    "type: ignore dropped": _replace("  # type: ignore[arg-type]", ""),
    "ruff noqa added as the first line": _prepend("# ruff: noqa"),
    "ruff noqa RUF001 dropped": _replace("# ruff: noqa: RUF001\n", ""),
    "ruff code added": _replace("RUF001", "RUF001, RUF002"),
    "mypy directive added": _prepend("# mypy: ignore-errors"),
    "pyright directive added": _prepend("# pyright: basic"),
    "isort directive added": _prepend("# isort: skip_file"),
    "only docstring replaced by ...": _replace('"""Raised when a tile is missing."""', "..."),
}


@pytest.mark.parametrize("edit", SAME.values(), ids=SAME.keys())
def test_same_code_accepts_a_prose_only_change(repo: Path, edit: Callable[[str], str]) -> None:
    base = _committed(repo, MODULE)
    (repo / "pkg" / "mod.py").write_text(edit(MODULE))

    run = _check(repo, "--same-code", base)

    assert run.code == 0, run.out + run.err


@pytest.mark.parametrize("edit", DIFFERS.values(), ids=DIFFERS.keys())
def test_same_code_reports_a_code_or_directive_change_naming_the_file(
    repo: Path, edit: Callable[[str], str]
) -> None:
    base = _committed(repo, MODULE)
    (repo / "pkg" / "mod.py").write_text(edit(MODULE))

    run = _check(repo, "--same-code", base)

    assert run.code == 1, run.out + run.err
    assert "pkg/mod.py" in run.out + run.err


def test_same_code_reports_a_new_file_naming_it(repo: Path) -> None:
    base = _committed(repo, MODULE)
    (repo / "pkg" / "new.py").write_text('"""New."""\n\nX = 1\n')
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "a new file")

    run = _check(repo, "--same-code", base)

    assert run.code == 1, run.out + run.err
    assert "pkg/new.py" in run.out + run.err


BLOCKS = '''\
"""Blocks beyond a function body."""


def walk(items: list[int]) -> int:
    total = 0
    try:
        total = sum(items)
    except ValueError:
        total = -1
    else:
        total += 1
    finally:
        items.clear()
    for item in items:
        total += item
    else:
        total += 2
    return total
'''


def test_same_code_ignores_bare_strings_under_except_else_finally_and_for_else(
    repo: Path,
) -> None:
    base = _committed(repo, BLOCKS)
    edited = BLOCKS
    for block, text in [
        ("except ValueError:", "Bad input."),
        ("else:\n        total += 1", "Fine."),
        ("finally:", "Always."),
        ("else:\n        total += 2", "Done."),
    ]:
        header, _, rest = block.partition("\n")
        edited = _replace(block, f'{header}\n        "{text}"\n{rest}'.rstrip("\n"))(edited)
    (repo / "pkg" / "mod.py").write_text(edited)

    run = _check(repo, "--same-code", base)

    assert run.code == 0, run.out + run.err


def test_same_code_still_sees_a_real_change_under_those_blocks(repo: Path) -> None:
    base = _committed(repo, BLOCKS)
    (repo / "pkg" / "mod.py").write_text(_replace("total += 2", "total += 3")(BLOCKS))

    run = _check(repo, "--same-code", base)

    assert run.code == 1, run.out + run.err
    assert "pkg/mod.py" in run.out + run.err
