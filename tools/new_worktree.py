#!/usr/bin/env python3
"""A worktree with a `.venv` of its own and a configured `build-pyext` (h18 T4).

    python3 tools/new_worktree.py <name> [--from <rev>] [--python <version>]
    python3 tools/new_worktree.py --existing <path> [--python <version>]

The first form makes `<main>/.claude/worktrees/<name>` on branch
`worktree-<name>` at `--from` (default `origin/master`, fetched first). Both
forms then make `<wt>/.venv` with uv, install this checkout into it editable
(`.[dev,codecs]` plus pybind11), configure `<wt>/build-pyext` for the
`_core` rebuild of `tools/rebuild_core.py`, and check that the
venv imports `tin_engine` from this worktree. A private (rasputin-lsm)
worktree also gets the core its pyproject pins, with the core's `codecs` extra
unless it declares imagecodecs itself; cmake is pointed at that core, and the
check is of `rasputin_lsm` from `<wt>/src`. Each step is skipped when its
product is already there, so `--existing` is safe to run again; it never
deletes anything.

Spec: docs/increments/h18-worktree-and-merge-tools.md §3. A refusal is exit 2
with its reason as the last line on stderr. Standard library only.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tomllib
from collections.abc import Callable, Sequence
from pathlib import Path

Run = Callable[..., subprocess.CompletedProcess[str]]
HERE = Path(__file__).resolve().parent
NAME = re.compile(r"[a-z0-9][a-z0-9._-]*")
#: Per layout: the package, its source directory, its extension, the rebuild line.
PUBLIC = ("tin_engine", "src_python", "_core", "python3 {wt}/tools/rebuild_core.py")
PRIVATE = ("rasputin_lsm", "src", "_lsm", "cmake --build {wt}/build-pyext -j --target _lsm")


def check_code(package: str, ext: str) -> str:
    return (f"import {package}, {package}.cli, {package}.{ext} as c; print({package}.__file__); "
            f"print({package}.cli.__file__); print(c.__file__)")  # fmt: skip


def private_install(wt: Path) -> list[str] | None:
    """What a rasputin-lsm checkout installs besides pybind11; None for any other."""
    path = wt / "pyproject.toml"
    project = tomllib.loads(path.read_text()).get("project", {}) if path.is_file() else {}
    if project.get("name") != "rasputin-lsm":
        return None
    extras, deps = project.get("optional-dependencies", {}), project.get("dependencies", [])
    names = [name for name in ("dev", "codecs") if name in extras]
    pin = next((d for d in deps if d.replace(" ", "").startswith("rasputin@")), None)
    if pin is None:
        raise Refusal(f"{path} pins no `rasputin @ <url>` core in its dependencies")
    if not any("imagecodecs" in d for d in [*deps, *(d for k in names for d in extras[k])]):
        pin = "rasputin[codecs] @" + pin.split("@", 1)[1]
    return ["-e", f".[{','.join(names)}]", pin]


class Refusal(Exception):  # noqa: N818  (a refusal is not an error of the tool)
    """Exit 2; the message is the last line on stderr, after the tool's prefix."""


def run(
    argv: Sequence[str], cwd: str | Path | None = None, *, capture: bool = True
) -> subprocess.CompletedProcess[str]:
    """The one route to a subprocess; a missing program is exit 127, not a traceback."""
    try:
        return subprocess.run(list(argv), cwd=cwd, capture_output=capture, text=True, check=False)
    except FileNotFoundError:
        return subprocess.CompletedProcess(list(argv), 127, "", f"{argv[0]}: not found\n")


def last_line(done: subprocess.CompletedProcess[str]) -> str:
    lines = [line for line in (done.stderr or done.stdout or "").splitlines() if line.strip()]
    return lines[-1].strip() if lines else f"exit {done.returncode}"


def common_dir(where: str | Path, run: Run) -> Path | None:
    done = run(["git", "-C", str(where), "rev-parse", "--path-format=absolute",
                "--git-common-dir"])  # fmt: skip
    return Path(done.stdout.strip()).resolve() if done.returncode == 0 else None


def check_checkout(path: Path, run: Run) -> None:
    """Refuse unless `path` is the top of a checkout of the repository holding this file."""
    top = run(["git", "-C", str(path), "rev-parse", "--show-toplevel"])
    ours, theirs = common_dir(HERE, run), common_dir(path, run)
    on_top = top.returncode == 0 and Path(top.stdout.strip()).resolve() == path.resolve()
    if not on_top or theirs is None or theirs != ours:
        raise Refusal(f"{path} is not the top of a checkout of this repository")


def kept(wt: Path) -> str:
    return f"The worktree is kept; finish with: python3 tools/new_worktree.py --existing {wt}"


def step(run: Run, name: str, argv: list[str], wt: Path, cwd: Path | None = None) -> str:
    """Run one setup step, its output passed through unless it is a probe; its stdout."""
    done = run(argv, cwd, capture=argv[-1].startswith("--"))
    if done.returncode != 0:
        sys.stderr.write(done.stderr or "")
        raise Refusal(f"{name} failed (exit {done.returncode}); its output is above. {kept(wt)}")
    return done.stdout or ""


def check(
    run: Run, wt: Path, py: str, layout: tuple[str, ...] = PUBLIC
) -> tuple[str, str, str] | str:
    """The step-7 import check: the three paths, or why it does not pass."""
    package, src, ext = layout[:3]
    done = run([py, "-c", check_code(package, ext)], wt)
    lines = [line.strip() for line in done.stdout.splitlines() if line.strip()]
    if done.returncode != 0 or len(lines) != 3:
        return f"the import check failed (exit {done.returncode}: {last_line(done)}). {kept(wt)}"
    root = wt.resolve()
    where = (root / src / package, root / src / package, root / ".venv")
    for module, path, under in zip((package, f"{package}.cli", f"{package}.{ext}"), lines,
                                   where, strict=True):  # fmt: skip
        if not Path(path).resolve().is_relative_to(under):
            return (f"{wt}/.venv imports {module} from {path}, not from {wt}; "
                    "this venv runs another checkout's code")  # fmt: skip
    return lines[0], lines[1], lines[2]


def set_up(run: Run, wt: Path, version: str) -> None:
    """Steps 4 to 7 of §3.2, each skipped when its product exists."""
    py, private = str(wt / ".venv" / "bin" / "python"), private_install(wt)
    layout = PUBLIC if private is None else PRIVATE
    if not (wt / ".venv").exists():
        step(run, "uv venv", ["uv", "venv", "--python", version, str(wt / ".venv")], wt)
    found, configured = check(run, wt, py, layout), (wt / "build-pyext" / "CMakeCache.txt").exists()
    if isinstance(found, str) or not configured:  # the configure needs pybind11 in this venv
        packages = ["-e", ".[dev,codecs]"] if private is None else private
        step(run, "uv pip install", ["uv", "pip", "install", "--python", py, *packages,
                                     "pybind11"], wt, wt)  # fmt: skip
    if not configured:
        cmakedir = step(run, "pybind11 --cmakedir", [py, "-m", "pybind11", "--cmakedir"], wt, wt)
        # The installed core's CMake config, which only scikit-build-core finds unaided.
        core = sorted((wt / ".venv" / "lib").glob("python*/site-packages/tin_engine"))
        root = [] if private is None else [f"-Drasputin_ROOT={core[-1] if core else ''}"]
        step(run, "cmake configure", ["cmake", "-S", str(wt), "-B", str(wt / "build-pyext"), *root,
             "-DCMAKE_BUILD_TYPE=Release", "-DRASPUTIN_BUILD_PYTHON=ON",
             "-DRASPUTIN_BUILD_TESTS=OFF", "-DRASPUTIN_HARDENING=ON", f"-DPYTHON_EXECUTABLE={py}",
             f"-DPython_EXECUTABLE={py}", f"-Dpybind11_DIR={cmakedir.strip()}"], wt)  # fmt: skip
    if isinstance(found, str):
        found = check(run, wt, py, layout)
        if isinstance(found, str):
            raise Refusal(found)
    for tool in ("ruff", "mypy"):
        step(run, f"{tool} --version", [py, "-m", tool, "--version"], wt, wt)
    branch = run(["git", "-C", str(wt), "rev-parse", "--abbrev-ref", "HEAD"]).stdout.strip()
    sha = run(["git", "-C", str(wt), "rev-parse", "--short", "HEAD"]).stdout.strip()
    python = run([py, "--version"], wt).stdout.strip()
    print(f"worktree  {wt}  (branch {branch} at {sha})\npython    {py}  ({python})\n"
          f"code      {layout[0]} and {layout[0]}.cli from {wt}/{layout[1]}\n"
          f"{layout[2]:<9} {found[2]}\nrebuild   {layout[3].format(wt=wt)}"
          "   (configured, not built)")  # fmt: skip


def probe(run: Run) -> None:
    for name in ("uv", "cmake"):
        done = run([name, "--version"])
        if done.returncode != 0:
            raise Refusal(f"{name} is not installed or does not run: {last_line(done)}")


def make(run: Run, name: str, rev: str) -> Path:
    """Steps 1 to 3 of §3.2: the refusals, the fetch and `git worktree add`."""
    if not NAME.fullmatch(name):
        raise Refusal(f"name \"{name}\": use lower-case letters, digits, '.', '_' and '-' only")
    common = common_dir(HERE, run)
    if common is None:
        raise Refusal(f"{HERE} is not in a checkout of this repository")
    main, wt = common.parent, common.parent / ".claude" / "worktrees" / name
    if wt.exists() or wt.is_symlink():
        raise Refusal(f"{wt} already exists; to give it a venv run: "
                      f"python3 tools/new_worktree.py --existing {wt}")  # fmt: skip
    branch = f"worktree-{name}"
    if run(["git", "-C", str(main), "rev-parse", "--verify", "--quiet",
            f"refs/heads/{branch}"]).returncode == 0:  # fmt: skip
        raise Refusal(f"branch {branch} already exists; pick another name")
    probe(run)
    if rev == "origin/master":
        fetch = run(["git", "-C", str(main), "fetch", "origin", "master"])
        if fetch.returncode != 0:
            raise Refusal(f"could not fetch origin/master ({last_line(fetch)}); "
                          "blocked on network: stop and hand back")  # fmt: skip
    add = run(["git", "-C", str(main), "worktree", "add", "-b", branch, str(wt), rev],
              capture=False)  # fmt: skip
    if add.returncode != 0:
        raise Refusal(f"git worktree add failed (exit {add.returncode}); its output is above. "
                      "No worktree was made")  # fmt: skip
    return wt


def main(argv: list[str] | None = None, run: Run = run) -> int:
    parser = argparse.ArgumentParser(prog="new_worktree", description=__doc__.split("\n")[0])
    parser.add_argument("name", nargs="?")
    parser.add_argument("--existing", type=Path)
    parser.add_argument("--from", dest="rev", default="origin/master")
    parser.add_argument("--python", default="3.14")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    if (args.name is None) == (args.existing is None):
        parser.error("give a <name> or --existing <path>, not both")
    try:
        if args.existing is not None:
            check_checkout(args.existing, run)
            probe(run)
            wt = args.existing
        else:
            wt = make(run, args.name, args.rev)
        set_up(run, wt, args.python)
    except Refusal as exc:
        print(f"new_worktree: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
