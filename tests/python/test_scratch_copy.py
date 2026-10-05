"""h16 T2: `tools/scratch_copy.py <rev> <dir>`.

Spec: `docs/increments/h16-harness-fixes.md` §2, T2. Each test runs the
script from a copy committed in a temporary repository, from that
repository, as `python3 tools/scratch_copy.py <rev> <dir>`. The repository
has an ignored `.venv` holding a stand-in `_core` and a `bin/python` that
runs this interpreter with a planted import hook: a meta-path finder named
like scikit-build-core's editable finder, which sends `tin_engine` to a decoy
package outside the copy. The printed command must drop that finder, or the
generated test it runs fails; `test_the_planted_finder_wins_without_the_drop`
shows the probe can fail.
"""

from __future__ import annotations

import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from harness_fixtures import REAL, clean_env, git

SCRIPT = "tools/scratch_copy.py"
CORE = "_core.cpython-312-darwin.so"
SITE = ".venv/lib/python3.12/site-packages/tin_engine"

FINDER = '''\
"""A stand-in for scikit-build-core's editable finder: tin_engine -> the decoy."""
import importlib.abc
import importlib.util
import sys

DECOY = {decoy!r}


class _EditableFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name == "tin_engine":
            return importlib.util.spec_from_file_location(
                name, DECOY + "/__init__.py", submodule_search_locations=[DECOY]
            )
        return None


sys.meta_path.insert(0, _EditableFinder())
'''

PROBE = """\
from pathlib import Path

import tin_engine


def test_tin_engine_comes_from_the_copy():
    where = Path(tin_engine.__file__).resolve()
    assert where.is_relative_to(Path({copy!r}).resolve()), where
"""


def env() -> dict[str, str]:
    out = clean_env()
    out.pop("PYTHONPATH", None)
    return out


def write(root: Path, relative: str, text: str) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """Two commits: `base` and `cxx` (a change under include/), then `py`
    (a change under src_python/ only). The `.venv` is ignored, as in a worktree."""
    base = tmp_path.resolve()
    root = base / "repo"
    root.mkdir()
    git(root, "init", "-q", "-b", "master")
    git(root, "config", "user.email", "t@example.invalid")
    git(root, "config", "user.name", "T")
    source = REAL / SCRIPT
    if source.exists():
        (root / "tools").mkdir()
        shutil.copy2(source, root / SCRIPT)
    write(root, ".gitignore", ".venv\n")
    write(root, "src_python/tin_engine/__init__.py", "")
    write(root, "src_python/tin_engine/mod.py", "VALUE = 'base'\n")
    write(root, "include/x.hpp", "int a;\n")
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "base")
    write(root, "include/x.hpp", "int a;\nint b;\n")
    git(root, "commit", "-q", "-am", "cxx")
    write(root, "src_python/tin_engine/mod.py", "VALUE = 'py'\n")
    git(root, "commit", "-q", "-am", "py")

    decoy = base / "decoy" / "tin_engine"
    write(decoy, "__init__.py", "")
    plant = base / "plant"
    write(plant, "sitecustomize.py", FINDER.format(decoy=str(decoy)))
    python = write(
        root, ".venv/bin/python",
        f'#!/bin/sh\nPYTHONPATH="{plant}" exec "{sys.executable}" "$@"\n',
    )  # fmt: skip
    python.chmod(python.stat().st_mode | stat.S_IXUSR)
    (root / SITE).mkdir(parents=True)
    (root / SITE / CORE).write_bytes(b"\x7fELF stand-in")
    return root


def run_copy(repo: Path, *args: str) -> subprocess.CompletedProcess[str]:
    script = repo / SCRIPT
    if not script.exists():
        pytest.fail(f"{SCRIPT} is missing from the checkout")
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        cwd=repo,
        env=env(),
        timeout=120,
        check=False,
    )


def copied(repo: Path, rev: str, target: Path) -> subprocess.CompletedProcess[str]:
    result = run_copy(repo, rev, str(target))
    assert result.returncode == 0, f"exit {result.returncode}; stderr: {result.stderr}"
    return result


def refused(result: subprocess.CompletedProcess[str], says: str) -> None:
    assert result.returncode == 2, f"exit {result.returncode}; stderr: {result.stderr}"
    lines = [line for line in result.stderr.splitlines() if line.strip()]
    assert len(lines) == 1, f"not one stderr line: {result.stderr!r}"
    assert says in lines[0], f"the refusal does not say {says!r}: {lines[0]!r}"


# ---------------------------------------------------------------- 1. refusals


def test_a_directory_that_is_not_empty_is_refused(repo: Path, tmp_path: Path) -> None:
    target = tmp_path / "full"
    write(target, "keep.txt", "mine\n")
    refused(run_copy(repo, "HEAD", str(target)), says="not empty")
    assert (target / "keep.txt").read_text() == "mine\n"
    assert sorted(p.name for p in target.iterdir()) == ["keep.txt"]


def test_a_directory_inside_the_repository_is_refused(repo: Path) -> None:
    target = repo / "scratch"
    refused(run_copy(repo, "HEAD", str(target)), says="inside")
    assert not (target / "include").exists()


def test_an_empty_directory_is_accepted(repo: Path, tmp_path: Path) -> None:
    target = tmp_path / "empty"
    target.mkdir()
    copied(repo, "HEAD", target)
    assert (target / "include/x.hpp").is_file()


def test_a_missing_directory_is_made(repo: Path, tmp_path: Path) -> None:
    target = tmp_path / "new"
    copied(repo, "HEAD", target)
    assert (target / "include/x.hpp").is_file()


# ---------------------------------------------------------------- 2. the copy


def test_the_copy_is_not_a_git_work_tree(repo: Path, tmp_path: Path) -> None:
    target = tmp_path / "copy"
    copied(repo, "HEAD", target)
    assert not (target / ".git").exists()
    probe = subprocess.run(
        ["git", "-C", str(target), "rev-parse", "--is-inside-work-tree"],
        capture_output=True, text=True, env=env(), check=False,
    )  # fmt: skip
    assert probe.returncode != 0, f"{target} is inside a work tree: {probe.stdout}"


def test_the_copy_holds_the_revision_not_the_working_tree(repo: Path, tmp_path: Path) -> None:
    write(repo, "src_python/tin_engine/mod.py", "VALUE = 'uncommitted'\n")
    write(repo, "untracked.txt", "u\n")
    target = tmp_path / "copy"
    copied(repo, "HEAD~2", target)
    assert (target / "src_python/tin_engine/mod.py").read_text() == "VALUE = 'base'\n"
    assert (target / "include/x.hpp").read_text() == "int a;\n"
    assert not (target / "untracked.txt").exists()


def test_replace_refs_do_not_change_the_copy(repo: Path, tmp_path: Path) -> None:
    planted = tmp_path / "planted"
    planted.write_text("VALUE = 'planted'\n")
    fake = git(repo, "hash-object", "-w", str(planted)).strip()
    real = git(repo, "rev-parse", "HEAD:src_python/tin_engine/mod.py").strip()
    git(repo, "replace", real, fake)
    assert "planted" in git(repo, "show", "HEAD:src_python/tin_engine/mod.py")
    target = tmp_path / "copy"
    copied(repo, "HEAD", target)
    assert (target / "src_python/tin_engine/mod.py").read_text() == "VALUE = 'py'\n"


# ---------------------------------------------------------------- 3. _core and warnings


def test_the_built_core_is_copied_into_the_package(repo: Path, tmp_path: Path) -> None:
    target = tmp_path / "copy"
    result = copied(repo, "HEAD", target)
    assert (target / "src_python/tin_engine" / CORE).read_bytes() == b"\x7fELF stand-in"
    assert result.stderr.strip() == "", f"a clean copy warned: {result.stderr!r}"


def test_no_built_core_warns_and_still_copies(repo: Path, tmp_path: Path) -> None:
    (repo / SITE / CORE).unlink()
    target = tmp_path / "copy"
    result = copied(repo, "HEAD", target)
    assert "_core" in result.stderr, f"no warning about the missing _core: {result.stderr!r}"
    assert (target / "src_python/tin_engine/mod.py").is_file()
    assert len(result.stdout.strip().splitlines()) == 1


def test_a_revision_before_a_cxx_change_warns_the_core_is_from_head(
    repo: Path, tmp_path: Path
) -> None:
    result = copied(repo, "HEAD~2", tmp_path / "copy")
    assert "HEAD" in result.stderr, f"no warning that _core was built from HEAD: {result.stderr!r}"
    assert "_core" in result.stderr


def test_a_revision_with_the_same_cxx_does_not_warn(repo: Path, tmp_path: Path) -> None:
    """HEAD~1 differs from HEAD under src_python/ only."""
    result = copied(repo, "HEAD~1", tmp_path / "copy")
    assert result.stderr.strip() == "", result.stderr


# ---------------------------------------------------------------- 4. the printed command


def runnable(printed: str, test_file: Path) -> str:
    """The printed command with `tests/python/` replaced by `test_file`, or, if
    the command names no test path, with `test_file` appended."""
    line = printed.strip()
    for tail in (" tests/python/", " tests/python"):
        if line.endswith(tail):
            return line.removesuffix(tail) + f" {test_file}"
    return f"{line} {test_file}"


def test_the_printed_command_is_one_line_naming_the_copy_and_the_venv(
    repo: Path, tmp_path: Path
) -> None:
    target = tmp_path / "copy"
    lines = copied(repo, "HEAD", target).stdout.strip().splitlines()
    assert len(lines) == 1, lines
    assert f"cd {target}" in lines[0]
    assert str(repo / ".venv/bin/python") in lines[0]


def test_the_planted_finder_wins_without_the_drop(repo: Path, tmp_path: Path) -> None:
    """The control: with the copy merely first on sys.path, the decoy is imported."""
    target = tmp_path / "copy"
    copied(repo, "HEAD", target)
    program = (
        f"import sys; sys.path.insert(0, {str(target / 'src_python')!r}); "
        "import tin_engine; print(tin_engine.__file__)"
    )
    result = subprocess.run(
        [str(repo / ".venv/bin/python"), "-c", program],
        capture_output=True, text=True, env=env(), check=True,
    )  # fmt: skip
    assert "decoy" in result.stdout, result.stdout


def test_the_printed_command_runs_pytest_against_the_copy(repo: Path, tmp_path: Path) -> None:
    target = tmp_path / "copy"
    printed = copied(repo, "HEAD", target).stdout
    probe = write(tmp_path / "gen", "test_where.py", PROBE.format(copy=str(target)))
    result = subprocess.run(
        ["/bin/sh", "-c", runnable(printed, probe)],
        capture_output=True, text=True, cwd=tmp_path, env=env(), timeout=120, check=False,
    )  # fmt: skip
    assert result.returncode == 0, f"exit {result.returncode}\n{result.stdout}\n{result.stderr}"
    assert "1 passed" in result.stdout, result.stdout
