"""`tools/rebuild_core.py`: rebuild the worktree's `_core` and install it where its venv imports it.

docs/increments/h36-rule-cut.md §3 and §6, tests 1-4. Each test runs a copy of
the tool in a temporary tree: `build-pyext/` stands for the configured build
directory, a fake `cmake` on `PATH` records its arguments and plays the build,
and `.venv` is a real virtual environment (no pip) whose site-packages holds a
stand-in `tin_engine` package with an installed `_core<suffix>` file. That
file is not a loadable extension, so the tool must find it without importing
it.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import sysconfig
from dataclasses import dataclass
from pathlib import Path

import pytest

pytestmark = pytest.mark.harness

REAL = Path(__file__).resolve().parents[2]
TOOL = "tools/rebuild_core.py"
SUFFIX = sysconfig.get_config_var("EXT_SUFFIX")
#: A suffix no interpreter running this suite has, for the wrong-build case.
OTHER_SUFFIX = ".cpython-39-other-os.so"
INSTALLED = b"the installed, stale extension\n"
BUILT = b"the freshly built extension\n"

FAKE_CMAKE = """\
#!/bin/sh
# Record the arguments, write the configured artefact, exit as configured.
printf '%s\\n' "$@" > "{log}"
if [ -f "{artefact_name}" ]; then cp "{artefact_bytes}" "{build}/$(cat "{artefact_name}")"; fi
exit "$(cat "{status}")"
"""


@dataclass
class Tree:
    root: Path
    bin: Path
    destination: Path

    @property
    def build(self) -> Path:
        return self.root / "build-pyext"

    @property
    def cmake_log(self) -> Path:
        return self.bin / "cmake.args"

    def cmake_plays(self, status: int, writes: str | None) -> None:
        """The fake build exits with `status`; it writes `build-pyext/<writes>` when given."""
        (self.bin / "status").write_text(f"{status}\n")
        if writes is not None:
            (self.bin / "artefact-name").write_text(writes)

    def run(self) -> subprocess.CompletedProcess[str]:
        env = {k: v for k, v in os.environ.items()
               if k not in ("VIRTUAL_ENV", "PYTHONPATH", "PYTHONHOME")}  # fmt: skip
        env["PATH"] = f"{self.bin}{os.pathsep}{env.get('PATH', '')}"
        return subprocess.run(
            [sys.executable, str(self.root / TOOL)],
            capture_output=True,
            text=True,
            cwd=self.root,
            env=env,
            timeout=120,
            check=False,
        )


def _lines(result: subprocess.CompletedProcess[str]) -> list[str]:
    return [line for line in (result.stdout + result.stderr).splitlines() if line.strip()]


@pytest.fixture
def tree(tmp_path: Path) -> Tree:
    root = tmp_path / "wt"
    (root / "tools").mkdir(parents=True)
    source = REAL / TOOL
    assert source.exists(), f"{TOOL} is missing from the checkout"
    shutil.copy2(source, root / TOOL)

    build = root / "build-pyext"
    build.mkdir()
    (build / "CMakeCache.txt").write_text("CMAKE_PROJECT_NAME:STATIC=rasputin\n")

    subprocess.run([sys.executable, "-m", "venv", "--without-pip", str(root / ".venv")],
                   check=True, capture_output=True, timeout=120)  # fmt: skip
    purelib = subprocess.run(
        [str(root / ".venv" / "bin" / "python"), "-c",
         "import sysconfig; print(sysconfig.get_paths()['purelib'])"],
        check=True, capture_output=True, text=True, timeout=60,
    ).stdout.strip()  # fmt: skip
    package = Path(purelib) / "tin_engine"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    destination = package / f"_core{SUFFIX}"
    destination.write_bytes(INSTALLED)

    fake = tmp_path / "bin"
    fake.mkdir()
    (fake / "artefact-bytes").write_bytes(BUILT)
    cmake = fake / "cmake"
    cmake.write_text(FAKE_CMAKE.format(
        log=fake / "cmake.args", artefact_name=fake / "artefact-name",
        artefact_bytes=fake / "artefact-bytes", build=build, status=fake / "status",
    ))  # fmt: skip
    cmake.chmod(0o755)
    made = Tree(root=root, bin=fake, destination=destination)
    made.cmake_plays(0, f"_core{SUFFIX}")
    return made


# --- 1. not configured ---------------------------------------------------------------------


def test_without_a_configured_build_it_refuses_naming_new_worktree_and_runs_no_build(
    tree: Tree,
) -> None:
    (tree.build / "CMakeCache.txt").unlink()
    result = tree.run()
    assert result.returncode == 2, result.stdout + result.stderr
    [line] = _lines(result)
    assert line.startswith("rebuild_core:")
    assert "tools/new_worktree.py" in line
    assert not tree.cmake_log.exists(), "cmake ran without a configured build directory"
    assert tree.destination.read_bytes() == INSTALLED


# --- 2. the build fails --------------------------------------------------------------------


@pytest.mark.parametrize("status", [1, 3])
def test_a_failed_build_returns_its_status_and_copies_nothing(tree: Tree, status: int) -> None:
    # The fake build still leaves an artefact, so a tool that copied anyway is caught.
    tree.cmake_plays(status, f"_core{SUFFIX}")
    result = tree.run()
    assert result.returncode == status, result.stdout + result.stderr
    [line] = _lines(result)
    assert line.startswith("rebuild_core:")
    assert str(status) in line
    assert tree.destination.read_bytes() == INSTALLED


# --- 3. the build passes -------------------------------------------------------------------


def test_a_passing_build_runs_the_core_target_and_installs_the_built_file(tree: Tree) -> None:
    result = tree.run()
    assert result.returncode == 0, result.stdout + result.stderr

    args = tree.cmake_log.read_text().splitlines()
    assert args[0] == "--build"
    assert Path(args[1]).resolve() == tree.build.resolve()
    assert args[2:] == ["-j", "--target", "_core"]

    assert tree.destination.read_bytes() == BUILT
    assert tree.destination.read_bytes() == (tree.build / f"_core{SUFFIX}").read_bytes()

    [line] = _lines(result)
    digest = hashlib.sha256(BUILT).hexdigest()[:12]
    assert digest in line.split()
    assert str(tree.destination) in line or str(tree.destination.resolve()) in line


# --- 4. the built file is not for this venv ------------------------------------------------


def test_a_build_without_the_venvs_suffix_exits_2_and_copies_nothing(tree: Tree) -> None:
    tree.cmake_plays(0, f"_core{OTHER_SUFFIX}")
    result = tree.run()
    assert result.returncode == 2, result.stdout + result.stderr
    assert (tree.build / f"_core{OTHER_SUFFIX}").exists()
    assert not (tree.build / f"_core{SUFFIX}").exists()
    assert tree.destination.read_bytes() == INSTALLED
    [line] = _lines(result)
    assert line.startswith("rebuild_core:")


def test_a_venv_without_an_installed_core_exits_2_and_installs_nothing(tree: Tree) -> None:
    # §3: "exit 2 if either is missing"; here the destination is the missing one.
    tree.destination.unlink()
    result = tree.run()
    assert result.returncode == 2, result.stdout + result.stderr
    assert not tree.destination.exists()
    [line] = _lines(result)
    assert line.startswith("rebuild_core:")
