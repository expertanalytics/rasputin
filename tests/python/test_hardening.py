"""Release hardening as the user sees it: increment 24, T4 to T7.

`docs/increments/24-release-hardening.md` sections 3 and 12. The bounds-checked
build is the default; which build is installed is read from the loaded
extension, ``tin_engine._core.hardening``, and repeated on line two of
``rasputin version`` and as one line of the ``--stats`` report. The strings are
as Ola ruled them (Q2): ``bounds checks: on (libc++ fast)`` and
``bounds checks: off``.

The ``stats.render`` half of T7 is in ``test_stats.py``, beside the golden
report; T8 (``tools/bench.py``) is in ``test_bench.py``.

T6, the stale build-directory trap, takes two real non-editable installs with
``--no-build-isolation`` into one throwaway virtual environment. Not for every
leg and not for a local ``pytest``, so it runs only with
``RASPUTIN_INSTALL_TRAP=1``, and **one CI leg** sets that: CI rather than a
manual check, because it is what keeps the ``pyproject.toml`` default and
scikit-build-core's behaviour checked after a release of it changes. The
command, locally or in that leg::

    RASPUTIN_INSTALL_TRAP=1 pytest --no-cov tests/python/test_hardening.py -k trap

Not invariant-critical, so no mutation round.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tomllib
import venv
from pathlib import Path

import pytest
from typer.testing import CliRunner

import tin_engine._core as core
from harness_fixtures import Tool
from tin_engine import cli, installed_version
from tin_engine.cli import app

REPO = Path(__file__).resolve().parents[2]

#: h11's classifier, ``tools/ci_changes.py``: the trap's copy leaves out the
#: files its ``is_prose`` accepts (``copy_source_tree``).
ci = Tool("ci_changes")

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

#: What a default (RASPUTIN_HARDENING=ON) build reports, keyed on the platform's
#: standard library: libc++ on macOS, libstdc++ on Linux.
HARDENED = {"darwin": "libc++ fast", "linux": "libstdc++ assertions"}


@pytest.fixture
def platform_mode() -> str:
    if sys.platform not in HARDENED:
        pytest.skip(f"no hardened mode is named for {sys.platform}")
    return HARDENED[sys.platform]


@pytest.fixture
def unchecked(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make the loaded extension look like an OFF build.

    `cli.py` may read the mode as ``_core.hardening`` at call time or import the
    name ``hardening`` from ``_core`` into its own namespace; both names are
    patched, so the test pins the behaviour and not which of the two it reads.
    """
    monkeypatch.setattr(core, "hardening", "none", raising=False)
    monkeypatch.setattr(cli, "hardening", "none", raising=False)


def version_lines() -> list[str]:
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0, result.output
    return result.stdout.splitlines()


def stats_report(tmp_path: Path) -> str:
    """``rasputin mesh catchment --flat --stats -``: the report after the path line."""
    out = tmp_path / "c.vtk"
    result = runner.invoke(app, ["mesh", "catchment", "--flat", "--out", str(out), "--stats", "-"])
    assert result.exit_code == 0, result.output
    _path_line, report = result.stdout.split("\n", 1)
    return report


def bounds_lines_before_first_section(report: str) -> tuple[list[str], list[str]]:
    """Every ``bounds checks:`` line, and those of them above the first ``## ``."""
    lines = report.splitlines()
    first_section = next(i for i, line in enumerate(lines) if line.startswith("## "))
    every = [line for line in lines if line.startswith("bounds checks:")]
    above = [line for line in lines[:first_section] if line.startswith("bounds checks:")]
    return every, above


class TestCoreAttribute:
    """T4: the installed, default build is the hardened one."""

    def test_the_default_install_reports_the_platforms_hardened_mode(
        self, platform_mode: str
    ) -> None:
        assert getattr(core, "hardening", None) == platform_mode

    def test_the_attribute_is_a_str(self) -> None:
        assert isinstance(getattr(core, "hardening", None), str)


class TestVersion:
    """T5: ``rasputin version`` prints the version, then the bounds-check line."""

    def test_two_lines_the_first_unchanged(self) -> None:
        lines = version_lines()
        assert len(lines) == 2, lines
        assert lines[0] == installed_version()

    def test_the_second_line_names_the_loaded_mode(self) -> None:
        lines = version_lines()
        assert len(lines) == 2, lines
        assert lines[1] == f"bounds checks: on ({getattr(core, 'hardening', None)})"

    def test_the_default_install_says_on(self, platform_mode: str) -> None:
        assert version_lines()[1:] == [f"bounds checks: on ({platform_mode})"]

    @pytest.mark.usefixtures("unchecked")
    def test_an_unchecked_build_says_off(self) -> None:
        lines = version_lines()
        assert lines == [installed_version(), "bounds checks: off"]


class TestStatsLine:
    """T7, the CLI half: the --stats report names the build it was made with."""

    def test_exactly_one_line_before_the_first_section_equal_to_version(
        self, tmp_path: Path
    ) -> None:
        every, above = bounds_lines_before_first_section(stats_report(tmp_path))
        assert len(every) == 1, every
        assert above == every
        assert every == version_lines()[1:]

    @pytest.mark.usefixtures("unchecked")
    def test_an_unchecked_build_says_off(self, tmp_path: Path) -> None:
        every, above = bounds_lines_before_first_section(stats_report(tmp_path))
        assert every == above == ["bounds checks: off"]


# ------------------------------------------------------------------------ T6
#
# Where the trap lives. scikit-build-core deletes a build directory's
# CMakeCache.txt when the scikit_build_core package it runs from is not the one
# recorded there (scikit_build_core/cmake.py, "New isolated environment ...,
# clearing cache"). Under pip's default build isolation every install gets a
# fresh environment, so the cache never survives and an OFF install cannot leak:
# a T6 run that way passes with or without pyproject.toml's
# [tool.scikit-build.cmake.define] block and so tests nothing. The cache
# survives only when both installs build with the same scikit-build-core:
# ``--no-build-isolation`` from one environment, which is the everyday
# developer's rebuild. So both installs go into one venv that holds the build
# requirements, and a sentinel entry planted in the cache after the first
# install proves, after the second, that the cache really was reused.

SENTINEL = "RASPUTIN_T6_SENTINEL:STRING=kept"


def copy_source_tree(destination: Path) -> Path:
    """The tracked and untracked-but-not-ignored files of this checkout, so the
    copy has no build directory of its own until the first install makes one.

    Prose is left out: every path ``tools/ci_changes.py``'s ``is_prose`` accepts
    (h11, ``docs/increments/h11-ci-path-filter.md`` section 6.1, point 2).
    ``shutil.copy2`` opens each file in this process, so copying prose would be
    a prose read, which ``conftest.py``'s session-end check (h11's T5) fails on;
    and the build must not need prose, because a prose-only pull request runs
    no build. The ``NOT_PROSE`` files, ``README.md`` among them, are still
    copied. The installs run in a ``pip`` subprocess that the hook cannot see,
    so this trap is where a build that starts to need a prose file shows up: as
    a failed install, whose message names ``NOT_PROSE``."""
    listed = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-co", "--exclude-standard", "-z"],
        capture_output=True,
        check=True,
    ).stdout.decode()
    for name in filter(None, listed.split("\0")):
        if ci.is_prose(name):
            continue
        source = REPO / name
        if source.is_file():
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    return destination


def build_environment(env_dir: Path) -> Path:
    """A venv holding ``[build-system].requires``, read from this checkout's
    ``pyproject.toml``; its python. uv when present (it serves the wheels from
    its cache, offline after the first fetch), pip otherwise."""
    venv.EnvBuilder(with_pip=True).create(env_dir)
    python = env_dir / "bin" / "python"
    requires = tomllib.loads((REPO / "pyproject.toml").read_text())["build-system"]["requires"]
    uv = shutil.which("uv")
    argv = (
        [uv, "pip", "install", "--quiet", "--python", str(python), *requires]
        if uv
        else [str(python), "-m", "pip", "install", "--quiet", *requires]
    )
    done = subprocess.run(argv, capture_output=True, text=True, check=False, timeout=600)
    assert done.returncode == 0, done.stderr[-4000:]
    return python


#: The one sentence h11 section 6.1 adds to a failed install's message.
NOT_COPIED = (
    "\nIf the build now needs a file the source copy leaves out (prose, as"
    " tools/ci_changes.py decides), that file belongs in NOT_PROSE."
)


def install(python: Path, source: Path, *settings: str) -> Path:
    """A non-editable, non-isolated ``pip install`` of ``source``; the installed
    extension's path. Dependencies are not installed: the extension is loaded by
    path below, which needs none of them."""
    argv = [str(python), "-m", "pip", "install", "--no-build-isolation", "--no-deps"]
    argv += ["--quiet", str(source), *settings]
    done = subprocess.run(argv, capture_output=True, text=True, check=False, timeout=1800)
    assert done.returncode == 0, done.stderr[-4000:] + NOT_COPIED
    (built,) = python.parent.parent.glob("lib/python*/site-packages/tin_engine/_core*.so")
    return built


def loaded_mode(extension: Path) -> str:
    """``hardening`` of the extension at ``extension``, imported by path in a
    fresh interpreter, so neither this session's ``_core`` nor the editable
    install's finder is involved."""
    probe = (
        "import importlib.util, sys\n"
        f"spec = importlib.util.spec_from_file_location('tin_engine._core', {str(extension)!r})\n"
        "module = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(module)\n"
        "print(getattr(module, 'hardening', '<no attribute>'))\n"
    )
    done = subprocess.run(
        [sys.executable, "-I", "-c", probe], capture_output=True, text=True, check=False
    )
    assert done.returncode == 0, done.stderr
    return done.stdout.strip()


@pytest.mark.skipif(
    os.environ.get("RASPUTIN_INSTALL_TRAP") != "1",
    reason="two real installs: set RASPUTIN_INSTALL_TRAP=1 (one CI leg does)",
)
def test_trap_an_off_install_does_not_leak_into_the_next_default_install(
    tmp_path: Path, platform_mode: str
) -> None:
    """T6: section 3's stale build-directory trap, and the override precedence.

    One source tree, so one ``build/{wheel_tag}``; one environment, so its
    CMakeCache.txt survives between installs. A:
    ``-C cmake.define.RASPUTIN_HARDENING=OFF``, which must override the ON in
    ``pyproject.toml``. B: no setting, which must not pick up A's cached OFF.
    Fails if pyproject.toml stops passing ON explicitly.
    """
    source = copy_source_tree(tmp_path / "src")
    python = build_environment(tmp_path / "venv")

    off = install(python, source, "--config-settings", "cmake.define.RASPUTIN_HARDENING=OFF")
    assert loaded_mode(off) == "none", "the command-line OFF did not override pyproject's ON"
    (cache,) = (source / "build").glob("*/CMakeCache.txt")
    assert "RASPUTIN_HARDENING:BOOL=OFF" in cache.read_text()
    cache.write_text(cache.read_text() + SENTINEL + "\n")

    default = install(python, source)
    assert SENTINEL in cache.read_text(), (
        "scikit-build-core cleared the cache between the installs: the trap was not exercised"
    )
    assert loaded_mode(default) == platform_mode, "A's cached OFF leaked into the default install"
