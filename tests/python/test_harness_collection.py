"""The collection hook that lets the harness job run without the package (h17 §4b).

The spec is `docs/increments/h17-ci-test-time.md` §4b and §9 A. pytest imports
every test file before `-m` deselects, so a job that runs `-m harness` without
`tin_engine` installed would fail at collection on every product test file.
`pytest_ignore_collect` in `tests/python/conftest.py` closes that: when the mark
expression is exactly `harness`, a `test_*.py` file whose text lacks
`mark.harness` is not collected (so never imported); otherwise pytest decides.

Each case runs pytest in a subprocess over a fixture directory, with the suite's
conftest loaded as a plugin (`-p conftest`, `tests/python` on the path) and
`addopts` cleared, as T5 in `test_ci_changes.py` does. The directory, not a
file, is passed: pytest never applies `pytest_ignore_collect` to a path named on
the command line, and CI passes none (it collects `testpaths`).

The unmarked file writes a sentinel before its failing import, so "never
imported" is read from the sentinel's absence, not inferred from a green run.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from harness_fixtures import REAL, clean_env

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

CONFTEST = REAL / "tests" / "python" / "conftest.py"
MISSING = "h17_no_such_module"

MARKED = """\
import pytest

pytestmark = pytest.mark.harness


def test_marked() -> None:
    assert True
"""

MARKED_IMPORTING_MISSING = f"""\
import pytest

import {MISSING}

pytestmark = pytest.mark.harness


def test_marked() -> None:
    assert True
"""

UNMARKED_IMPORTING_MISSING = f"""\
from pathlib import Path

(Path(__file__).parent / "imported.sentinel").write_text("imported")

import {MISSING}


def test_unmarked() -> None:
    assert True
"""

Run = Callable[[dict[str, str], list[str]], subprocess.CompletedProcess[str]]


@pytest.fixture
def run_pytest(tmp_path: Path) -> Run:
    """Write `files` into a fresh directory and run pytest over it with `args`."""

    def run(files: dict[str, str], args: list[str]) -> subprocess.CompletedProcess[str]:
        assert CONFTEST.exists(), "tests/python/conftest.py does not exist"
        suite = tmp_path / "suite"
        suite.mkdir()
        for name, text in files.items():
            (suite / name).write_text(text)
        env = {**clean_env(), "PYTHONPATH": str(CONFTEST.parent)}
        command = [sys.executable, "-m", "pytest", "-q", "-p", "conftest", "-p", "no:cacheprovider"]
        command += [
            "-o",
            "addopts=",
            "-o",
            "markers=harness: h17 fixture",
            "--rootdir",
            str(tmp_path),
        ]
        return subprocess.run(
            [*command, *args, str(suite)],
            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120,
        )  # fmt: skip

    return run


def sentinel(tmp_path: Path) -> Path:
    return tmp_path / "suite" / "imported.sentinel"


def test_under_m_harness_an_unmarked_file_is_never_imported(
    run_pytest: Run, tmp_path: Path
) -> None:
    done = run_pytest(
        {"test_marked.py": MARKED, "test_unmarked.py": UNMARKED_IMPORTING_MISSING},
        ["-m", "harness"],
    )
    output = done.stdout + done.stderr
    assert done.returncode == 0, output
    assert "1 passed" in output, output
    assert MISSING not in output, output
    assert not sentinel(tmp_path).exists(), "the unmarked file was imported under -m harness"


def test_under_m_harness_a_marked_file_that_fails_to_import_still_errors(
    run_pytest: Run,
) -> None:
    """The loud mode stays: a harness test that needs a missing module fails collection."""
    done = run_pytest({"test_marked.py": MARKED_IMPORTING_MISSING}, ["-m", "harness"])
    output = done.stdout + done.stderr
    assert done.returncode == pytest.ExitCode.INTERRUPTED, output
    assert f"ModuleNotFoundError: No module named '{MISSING}'" in output, output


@pytest.mark.parametrize("args", [["-m", "not harness"], []], ids=["m-not-harness", "no-m"])
def test_under_any_other_expression_an_unmarked_file_is_imported(
    run_pytest: Run, tmp_path: Path, args: list[str]
) -> None:
    """Only the exact expression `harness` ignores files; the product legs collect as before."""
    done = run_pytest(
        {"test_marked.py": MARKED, "test_unmarked.py": UNMARKED_IMPORTING_MISSING}, args
    )
    output = done.stdout + done.stderr
    assert sentinel(tmp_path).exists(), output
    assert done.returncode == pytest.ExitCode.INTERRUPTED, output
    assert f"ModuleNotFoundError: No module named '{MISSING}'" in output, output
