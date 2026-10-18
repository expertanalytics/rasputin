"""The test session never reads the user's config or home (increment 68).

The spec is `docs/increments/68-pytest-home.md` §4 and §5. The autouse fixtures
that give each test an empty `XDG_CONFIG_HOME` and `XDG_CACHE_HOME` are
function-scoped, so a module-, class- or session-scoped fixture is set up before
them and sees the environment pytest was started with: the user's `HOME`, their
sources file, and their `RASPUTIN_` settings. A `pytest_configure` hook in
`tests/python/conftest.py` moves `HOME` and the four XDG directories into one
empty temporary directory and removes the package's three user settings, before
collection and before any fixture.

Two tests here are product tests (they call `tin_engine`, imported inside the
test, so the harness job, which runs without the package, can import this
file); two are harness tests that run an inner pytest in a subprocess, with the
suite's conftest loaded as a plugin (`-p conftest`, `tests/python` on the path)
and `addopts` cleared, as `test_harness_collection.py` does.
"""

from __future__ import annotations

import os
import pwd
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from harness_fixtures import REAL, clean_env

CONFTEST = REAL / "tests" / "python" / "conftest.py"

#: The package's user settings, which the session hook removes (§4 item 1).
USER_SETTINGS = ("RASPUTIN_DATA", "RASPUTIN_THREADS", "RASPUTIN_MCP_OUT_ROOT")

#: Three bytes that `tomllib` refuses: an unreadable sources file.
BAD_TOML = "[[["

#: The account's home directory, from the password database, not from `HOME`.
ACCOUNT_HOME = Path(pwd.getpwuid(os.getuid()).pw_dir).resolve()


def is_under(path: Path, root: Path) -> bool:
    return path.resolve().is_relative_to(root.resolve())


# --- product tests ----------------------------------------------------------


@pytest.fixture(scope="module")
def module_view() -> dict[str, Any]:
    """What a module-scoped fixture sees: set up before every function-scoped
    autouse fixture, so only the session hook stands between it and the
    environment pytest was started with."""
    from tin_engine.sources import user_source_candidates

    return {
        "HOME": os.environ.get("HOME"),
        "candidates": user_source_candidates(os.environ),
        "settings": {name: os.environ[name] for name in USER_SETTINGS if name in os.environ},
    }


def test_a_module_fixture_sees_no_user_sources_file(module_view: dict[str, Any]) -> None:
    home = module_view["HOME"]
    assert home is not None, "HOME is unset in a module fixture; the hook sets it"
    assert not is_under(Path(home), ACCOUNT_HOME), (
        f"a module fixture sees HOME={home}, inside the account's home {ACCOUNT_HOME}"
    )
    candidates: tuple[Path, ...] = module_view["candidates"]
    assert candidates, "no candidate at all: the XDG configuration directory went missing"
    for candidate in candidates:
        assert not candidate.exists(), f"a module fixture would read {candidate}"
        assert not is_under(candidate, ACCOUNT_HOME), (
            f"sources candidate {candidate} is under the account's home {ACCOUNT_HOME}"
        )
    assert module_view["settings"] == {}, (
        f"user settings reach a module fixture: {module_view['settings']}"
    )


def test_a_monkeypatched_home_still_wins(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A test that sets its own HOME gets its own places, not the session's."""
    from tin_engine.sources import user_source_candidates

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    assert user_source_candidates(os.environ) == (
        tmp_path / ".config" / "rasputin" / "sources.toml",
    )


# --- harness tests: an inner pytest run ---------------------------------------

Run = Callable[[str, dict[str, str]], subprocess.CompletedProcess[str]]


@pytest.fixture
def run_inner(tmp_path: Path) -> Run:
    """Write `text` as the one test file of a fresh directory and run pytest
    over it, with the suite's conftest as a plugin and `overrides` laid over
    the outer environment."""

    def run(text: str, overrides: dict[str, str]) -> subprocess.CompletedProcess[str]:
        assert CONFTEST.exists(), "tests/python/conftest.py does not exist"
        suite = tmp_path / "suite"
        suite.mkdir()
        (suite / "test_inner.py").write_text(text)
        env = {**clean_env(), **overrides, "PYTHONPATH": str(CONFTEST.parent)}
        command = [sys.executable, "-m", "pytest", "-q", "-p", "conftest", "-p", "no:cacheprovider"]
        command += ["-o", "addopts=", "--rootdir", str(tmp_path), str(suite)]
        return subprocess.run(
            command, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120
        )

    return run


#: The inner test of `test_a_planted_user_config_...`: its module fixture
#: parses every sources file it can find, as the CLI does, so a planted
#: unreadable file is an error, the incident's shape.
INNER_READS_SOURCES = """\
import os
import tomllib
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def sources():
    places = [
        Path(os.environ.get("HOME", "/nonexistent")) / ".config" / "rasputin" / "sources.toml",
        Path(os.environ.get("XDG_CONFIG_HOME", "/nonexistent")) / "rasputin" / "sources.toml",
        Path(os.environ.get("RASPUTIN_DATA", "/nonexistent")) / "sources.toml",
    ]
    return {str(p): tomllib.loads(p.read_text()) for p in places if p.is_file()}


def test_reads_no_user_sources(sources):
    assert sources == {}
"""


@pytest.mark.harness
def test_a_planted_user_config_does_not_reach_a_module_fixture(
    run_inner: Run, tmp_path: Path
) -> None:
    home = tmp_path / "planted-home"
    config = tmp_path / "planted-config"
    data = tmp_path / "planted-data"
    for folder in (home / ".config" / "rasputin", config / "rasputin", data):
        folder.mkdir(parents=True)
        (folder / "sources.toml").write_text(BAD_TOML)
    done = run_inner(
        INNER_READS_SOURCES,
        {"HOME": str(home), "XDG_CONFIG_HOME": str(config), "RASPUTIN_DATA": str(data)},
    )
    output = done.stdout + done.stderr
    assert done.returncode == pytest.ExitCode.OK, output
    assert "1 passed" in output, output


#: The inner test of `test_the_install_trap_switch_...`: CI's trap step reads
#: the switch in a `skipif` at collection time, so it is read at import too.
INNER_READS_TRAP_SWITCH = """\
import os

AT_IMPORT = os.environ.get("RASPUTIN_INSTALL_TRAP")


def test_switch_is_set():
    assert AT_IMPORT == "1"
    assert os.environ.get("RASPUTIN_INSTALL_TRAP") == "1"
"""


@pytest.mark.harness
def test_the_install_trap_switch_reaches_the_tests(run_inner: Run) -> None:
    """The hook removes the user settings only, not every `RASPUTIN_` name:
    CI's "Hardening install trap" step would otherwise skip and stay green."""
    done = run_inner(INNER_READS_TRAP_SWITCH, {"RASPUTIN_INSTALL_TRAP": "1"})
    output = done.stdout + done.stderr
    assert done.returncode == pytest.ExitCode.OK, output
    assert "1 passed" in output, output
