"""Shared fixture for the h18 suites: a temporary repository and a recording `run`.

The spec is `docs/increments/h18-worktree-and-merge-tools.md`, §3.6 and §4.9.
Both tools have `main(argv, run)`, and every subprocess goes through `run(argv,
cwd, *, capture=True)`. The tests pass a `Recorder` as `run`: it runs `git`
for real, against a temporary repository whose `origin` is a bare repository
on disk (no network), and answers `uv`, `cmake` and the venv's Python from a
table the test gives. Any other program is a test failure naming the argv.

The repository is made under a symbolic link (`<tmp>/link -> <tmp>/real`), and
every path is handed to the tools through the link, unresolved. The tables
answer with resolved paths, as Python reports a module's file, so a tool that
compares paths without `Path.resolve()` fails on every platform, not only on
macOS (where `/var` is `/private/var`).

Each tool is loaded from its copy inside the temporary main checkout, never
from this checkout, so "this repository" (the tool's own checkout, or the
current directory, whichever the tool uses) is the temporary one and no test can
make a worktree or a merge in the real repository. The copy is taken from this
checkout's `tools/` when the repository is made; while a tool does not exist,
loading it fails the test, naming the file.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType

import pytest

from harness_fixtures import COPIED, GUARD_PUSH, REAL

#: The two tools of h18, copied into the temporary repository beside COPIED
#: (which holds tools/brief.py, for `brief.WRITES`, and the governance guard,
#: for `governed()`, with what they import).
TOOLS = ("tools/new_worktree.py", "tools/merge_master.py")

#: The trailer the tests pass; it matches §4.6's pattern.
TRAILER = "Co-Authored-By: Claude Test <noreply@example.invalid>"

#: No user or system git configuration (no signing, no hooks path, no rerere),
#: an identity for commits, and no editor or prompt that could wait for input.
ENV = {
    **{k: v for k, v in os.environ.items() if not k.startswith("GIT_")},
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_AUTHOR_NAME": "T",
    "GIT_AUTHOR_EMAIL": "t@example.invalid",
    "GIT_COMMITTER_NAME": "T",
    "GIT_COMMITTER_EMAIL": "t@example.invalid",
    "GIT_EDITOR": "false",
    "GIT_TERMINAL_PROMPT": "0",
}

GITIGNORE = ".venv/\nbuild-pyext/\n.claude/worktrees/\n__pycache__/\n"

#: The files of the root commit, besides the copied tools.
ROOT_FILES = {
    ".gitignore": GITIGNORE,
    "CLAUDE.md": "# rules\n\nline one\n",
    "notes.txt": "notes\n",
    "a.py": "x = 1\ny = 2\n",
    "b.py": "b = 1\n",
    "gone.py": "g = 1\n",
    "include/core.h": "#pragma once\n",
}


def git(repo: Path, *args: str) -> str:
    """Run git in `repo` with the suite's environment; its stdout. Fails on error."""
    done = subprocess.run(
        ["git", "-C", str(repo), *args],
        capture_output=True, text=True, env=ENV, stdin=subprocess.DEVNULL, check=False,
    )  # fmt: skip
    assert done.returncode == 0, f"git {' '.join(args)}: {done.stderr}"
    return done.stdout


def write(repo: Path, changes: dict[str, str | None]) -> None:
    """Write each file; `None` removes it from the index and the tree."""
    for name, text in changes.items():
        if text is None:
            git(repo, "rm", "-q", name)
        else:
            (repo / name).parent.mkdir(parents=True, exist_ok=True)
            (repo / name).write_text(text)
            git(repo, "add", name)


def commit(repo: Path, changes: dict[str, str | None], message: str) -> str:
    """One commit of `changes` in `repo`; its full hash."""
    write(repo, changes)
    git(repo, "commit", "-q", "-m", message)
    return head(repo)


def head(repo: Path) -> str:
    return git(repo, "rev-parse", "HEAD").strip()


def git_dir(repo: Path) -> Path:
    """The checkout's own git dir (a linked worktree's is under the main `.git`)."""
    return Path(git(repo, "rev-parse", "--absolute-git-dir").strip())


@dataclass
class Repos:
    """`seed` holds master's history; `origin` is bare; `main` is a clone of it."""

    seed: Path
    origin: Path
    main: Path

    def publish(self) -> None:
        """Make origin's master the seed's (a fetch into the bare repository, not a push)."""
        git(self.origin, "fetch", "-q", str(self.seed), "+master:master")

    def pr(self, number: int, changes: dict[str, str | None]) -> None:
        """A first-parent merge on master, `Merge pull request #<n> from x/pr-<n>`."""
        branch = f"pr-{number}"
        git(self.seed, "checkout", "-q", "-b", branch, "master")
        commit(self.seed, changes, f"change for #{number}")
        git(self.seed, "checkout", "-q", "master")
        git(self.seed, "merge", "-q", "--no-ff", "-m",
            f"Merge pull request #{number} from x/{branch}", branch)  # fmt: skip

    def master_commit(self, changes: dict[str, str | None], message: str) -> None:
        """A plain commit on master, with no PR subject."""
        commit(self.seed, changes, message)

    def worktree(self, name: str) -> Path:
        return self.main / ".claude" / "worktrees" / name


def make_repos(tmp: Path) -> Repos:
    """The three repositories, made through `<tmp>/link`, a link to `<tmp>/real`."""
    (tmp / "real").mkdir()
    (tmp / "link").symlink_to(tmp / "real", target_is_directory=True)
    base = tmp / "link"
    seed = base / "seed"
    seed.mkdir()
    git(seed, "init", "-q", "-b", "master")
    for relative in (*COPIED, *TOOLS):
        source = REAL / relative
        if source.exists():
            (seed / relative).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, seed / relative)
    for name, text in ROOT_FILES.items():
        (seed / name).parent.mkdir(parents=True, exist_ok=True)
        (seed / name).write_text(text)
    git(seed, "add", "-A")
    git(seed, "commit", "-q", "-m", "root")
    origin, main = base / "origin.git", base / "main"
    subprocess.run(["git", "clone", "-q", "--bare", str(seed), str(origin)],
                   check=True, capture_output=True, env=ENV)  # fmt: skip
    subprocess.run(["git", "clone", "-q", str(origin), str(main)],
                   check=True, capture_output=True, env=ENV)  # fmt: skip
    return Repos(seed, origin, main)


def stub_venv(wt: Path, python: str = "3.12") -> Path:
    """`<wt>/.venv/bin/python` (never run: the recorder answers it) and its package dir."""
    bin_dir = wt / ".venv" / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    py = bin_dir / "python"
    py.write_text("#!/bin/sh\necho 'h18: the stub venv python was run for real' >&2\nexit 99\n")
    py.chmod(0o755)
    (wt / ".venv" / "lib" / f"python{python}" / "site-packages" / "tin_engine").mkdir(
        parents=True, exist_ok=True
    )
    return py


def load(checkout: Path, relative: str) -> ModuleType:
    """Import `<checkout>/<relative>` under a name of its own, or fail naming the file."""
    path = checkout / relative
    if not path.exists():
        pytest.fail(f"{relative} does not exist in the temporary checkout")
    name = f"h18_{path.stem}_{abs(hash(str(path)))}"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def isolated_imports(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Iterable[None]:
    """Undo what loading a tool copy does to `sys.path` and `sys.modules`.

    A tool may put its own `tools/` on `sys.path` and import `brief` or a guard
    from there; those modules come from this test's temporary repository and
    must not be found by a later test.
    """
    monkeypatch.setattr(sys, "path", list(sys.path))
    before = set(sys.modules)
    yield
    real_tmp = str(tmp_path.resolve())
    for name in set(sys.modules) - before:
        origin = getattr(sys.modules[name], "__file__", None) or ""
        if name.startswith("h18_") or str(Path(origin).resolve()).startswith(real_tmp):
            del sys.modules[name]


# --- the recorder ---------------------------------------------------------------------


@dataclass
class Answer:
    returncode: int = 0
    stdout: str = ""
    stderr: str = ""


@dataclass
class Call:
    argv: list[str]
    cwd: Path | None
    capture: bool


def program(argv: list[str]) -> str:
    return Path(argv[0]).name


def is_venv_python(argv: list[str]) -> bool:
    path = Path(argv[0])
    return (
        path.name.startswith("python")
        and path.parent.name == "bin"
        and (path.parent.parent.name == ".venv")
    )


def git_sub(argv: list[str]) -> str | None:
    """The git subcommand of an argv, past `-C <dir>`, `-c <k=v>` and other options."""
    if program(argv) != "git":
        return None
    i = 1
    while i < len(argv) and argv[i].startswith("-"):
        i += 2 if argv[i] in ("-C", "-c", "--git-dir", "--work-tree", "--namespace") else 1
    return argv[i] if i < len(argv) else None


Table = Callable[[list[str], Path | None], Answer]


@dataclass
class Recorder:
    """The `run` the tools are given: git for real, the rest from `table`.

    `git_override` may answer a git call instead of running it (a failing
    fetch). Output of a call made with `capture=False` is written to this
    process's stdout and stderr, as a passed-through subprocess would write it.
    """

    table: Table
    git_override: Callable[[list[str]], Answer | None] | None = None
    calls: list[Call] = field(default_factory=list)

    def __call__(
        self,
        argv: Iterable[str | os.PathLike[str]],
        cwd: str | os.PathLike[str] | None = None,
        *,
        capture: bool = True,
    ) -> subprocess.CompletedProcess[str]:
        words = [os.fspath(a) for a in argv]
        here = None if cwd is None else Path(cwd)
        self.calls.append(Call(words, here, capture))
        if program(words) == "git":
            answer = self.git_override(words) if self.git_override else None
            if answer is None:
                done = subprocess.run(words, cwd=here, capture_output=True, text=True, env=ENV,
                                      stdin=subprocess.DEVNULL, check=False)  # fmt: skip
                answer = Answer(done.returncode, done.stdout, done.stderr)
        else:
            answer = self.table(words, here)
        if not capture:
            sys.stdout.write(answer.stdout)
            sys.stderr.write(answer.stderr)
            return subprocess.CompletedProcess(words, answer.returncode, "", "")
        return subprocess.CompletedProcess(words, answer.returncode, answer.stdout, answer.stderr)

    def argvs(self) -> list[list[str]]:
        return [c.argv for c in self.calls]

    def git_subs(self) -> list[str]:
        return [s for s in (git_sub(c.argv) for c in self.calls) if s is not None]


def failing_fetch(argv: list[str]) -> Answer | None:
    """A `git_override` under which every fetch fails as a network failure would."""
    if git_sub(argv) == "fetch":
        return Answer(128, "", "fatal: unable to access 'origin': h18 test network down\n")
    return None


@dataclass
class Outcome:
    code: int
    out: str
    err: str
    recorder: Recorder

    @property
    def text(self) -> str:
        return self.out + self.err

    def err_lines(self) -> list[str]:
        return [line for line in self.err.splitlines() if line.strip()]


def invoke(
    tool: ModuleType, args: list[str], recorder: Recorder, capsys: pytest.CaptureFixture[str]
) -> Outcome:
    """`tool.main(args, run=recorder)`; its exit code (returned or raised) and output."""
    capsys.readouterr()
    try:
        code = tool.main(args, run=recorder)
    except SystemExit as stop:
        code = stop.code if isinstance(stop.code, int) else 1
    out, err = capsys.readouterr()
    return Outcome(int(code), out, err, recorder)


# --- the push guard, loaded by path as harness_fixtures does --------------------------


def guard_push() -> ModuleType:
    path = REAL / GUARD_PUSH
    spec = importlib.util.spec_from_file_location("h18_guard_push", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def push_guard_findings(argvs: Iterable[list[str]]) -> list[tuple[list[str], object]]:
    """Every argv the push guard would ask about, with its reasons."""
    guard = guard_push()
    found: list[tuple[list[str], object]] = []
    for argv in argvs:
        reasons, why = guard.publishes(argv), guard.segment_why(argv)
        if reasons != [] or why is not None:
            found.append((argv, (reasons, why)))
    return found
