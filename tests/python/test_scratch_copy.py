"""h16 T2: `tools/scratch_copy.py <rev> <dir>`.

Spec: `docs/increments/h16-harness-fixes.md` §2, T2. Each test runs the
script from a copy committed in a temporary repository, from that
repository, as `python3 tools/scratch_copy.py <rev> <dir>`. The repository
has an ignored `.venv`, a real virtual environment of this interpreter,
holding a stand-in `_core` and, as an editable install does, a `.pth` file
that installs a meta-path finder named like scikit-build-core's editable
finder at every interpreter start, sending `tin_engine` to a decoy package
outside the copy. A second `.pth` line puts this interpreter's own
site-packages (pytest) on the path. The printed command must drop that
finder, in its own process and in any Python process a test starts (review
round 1, blocking 2), or the generated test it runs fails;
`test_the_planted_finder_wins_without_the_drop` shows the probe can fail.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from harness_fixtures import REAL, clean_env, git

SCRIPT = "tools/scratch_copy.py"
CORE = "_core.cpython-312-darwin.so"

FINDER = '''\
"""A stand-in for scikit-build-core's editable finder: tin_engine -> the decoy.

Imported by a `.pth` line, so it is installed at every interpreter start."""
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

#: A generated test whose child Python process must see the copy too: three
#: suites start such children (review round 1, blocking 2), and in a mutant run
#: a child importing the original code would report a false survivor.
CHILD_PROBE = """\
import subprocess
import sys
from pathlib import Path


def test_a_child_python_imports_tin_engine_from_the_copy():
    out = subprocess.run(
        [sys.executable, "-c", "import tin_engine; print(tin_engine.__file__)"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    where = Path(out).resolve()
    assert where.is_relative_to(Path({copy!r}).resolve()), where
"""


def env() -> dict[str, str]:
    out = clean_env()
    out.pop("PYTHONPATH", None)
    return out


def site_packages(root: Path) -> Path:
    """The `.venv`'s site-packages directory, whatever the Python version."""
    (found,) = (root / ".venv").glob("lib/python3.*/site-packages")
    return found


def write(root: Path, relative: str, text: str) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """Two commits: `base` and `cxx` (a change under include/), then `py`
    (a change under src_python/ only). The `.venv` is ignored, as in a
    worktree, and so is `.claude/worktrees/`, as in the main checkout."""
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
    write(root, ".gitignore", ".venv\n.claude/worktrees/\n")
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
    subprocess.run(
        [sys.executable, "-m", "venv", "--without-pip", str(root / ".venv")],
        capture_output=True, env=env(), check=True,
    )  # fmt: skip
    site = site_packages(root)
    write(site, "_plant_finder.py", FINDER.format(decoy=str(decoy)))
    write(site, "_editable_tin_engine.pth", "import _plant_finder\n")
    write(site, "_host_packages.pth", str(Path(pytest.__file__).resolve().parent.parent) + "\n")
    (site / "tin_engine").mkdir()
    (site / "tin_engine" / CORE).write_bytes(b"\x7fELF stand-in")
    return root


def run_copy(
    repo: Path, *args: str, environ: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    """The copy committed in `repo` (a checkout or a worktree), by path, from `repo`."""
    script = repo / SCRIPT
    if not script.exists():
        pytest.fail(f"{SCRIPT} is missing from the checkout")
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        cwd=repo,
        env=environ if environ is not None else env(),
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


@pytest.mark.parametrize("where", ["scratch", ".claude/scratch"])
def test_from_a_worktree_a_directory_inside_the_main_checkout_is_refused(
    repo: Path, where: str
) -> None:
    """Review round 1: run from a worktree under `.claude/worktrees/`, a target
    outside the worktree but inside the main checkout that holds it passed."""
    worktree = repo / ".claude/worktrees/w"
    git(repo, "worktree", "add", "-q", "--detach", str(worktree), "HEAD")
    target = repo / where
    refused(run_copy(worktree, "HEAD", str(target)), says="inside")
    assert not target.exists() or not any(target.iterdir())


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


#: A `tar` that extracts part of the archive and then fails, as a full disk
#: would, with two lines of its own on stderr.
FAILING_TAR = """\
#!/bin/sh
while [ $# -gt 0 ]; do
    if [ "$1" = "-C" ]; then shift; dir="$1"; fi
    shift
done
mkdir -p "$dir/src_python" && echo half > "$dir/src_python/partial.py"
cat > /dev/null
echo "tar: first error" >&2
echo "tar: second error" >&2
exit 1
"""


@pytest.mark.parametrize("existing", [False, True], ids=["missing", "empty"])
def test_a_failed_tar_is_one_stderr_line_and_leaves_nothing_behind(
    repo: Path, tmp_path: Path, existing: bool
) -> None:
    """Review round 1: tar's own lines and a traceback-like line reached stderr,
    and the target was left partly filled. Now: exit 2, one `scratch_copy:`
    line, and the target missing or empty, as it was before."""
    fake = write(tmp_path / "bin", "tar", FAILING_TAR)
    fake.chmod(0o755)
    target = tmp_path / "copy"
    if existing:
        target.mkdir()
    environ = {**env(), "PATH": f"{fake.parent}:{env()['PATH']}"}
    result = run_copy(repo, "HEAD", str(target), environ=environ)
    refused(result, says="tar")
    assert result.stderr.lstrip().startswith("scratch_copy:"), result.stderr
    assert not target.exists() or not any(target.iterdir()), sorted(
        str(p.relative_to(target)) for p in target.rglob("*")
    )


# ---------------------------------------------------------------- 3. _core and warnings


def test_the_built_core_is_copied_into_the_package(repo: Path, tmp_path: Path) -> None:
    target = tmp_path / "copy"
    result = copied(repo, "HEAD", target)
    assert (target / "src_python/tin_engine" / CORE).read_bytes() == b"\x7fELF stand-in"
    assert result.stderr.strip() == "", f"a clean copy warned: {result.stderr!r}"


def test_no_built_core_warns_and_still_copies(repo: Path, tmp_path: Path) -> None:
    (site_packages(repo) / "tin_engine" / CORE).unlink()
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


def test_a_python_process_started_by_a_test_imports_the_copy(repo: Path, tmp_path: Path) -> None:
    """Review round 1, blocking 2: the drop in the printed command's own process
    does not reach a child, which installs the editable finder again at start."""
    target = tmp_path / "copy"
    printed = copied(repo, "HEAD", target).stdout
    probe = write(tmp_path / "gen", "test_child.py", CHILD_PROBE.format(copy=str(target)))
    result = subprocess.run(
        ["/bin/sh", "-c", runnable(printed, probe)],
        capture_output=True, text=True, cwd=tmp_path, env=env(), timeout=120, check=False,
    )  # fmt: skip
    assert result.returncode == 0, f"exit {result.returncode}\n{result.stdout}\n{result.stderr}"
    assert "1 passed" in result.stdout, result.stdout
