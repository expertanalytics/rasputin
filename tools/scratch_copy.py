#!/usr/bin/env python3
"""A copy of one revision to run the suites against, outside any work tree (h16 T2).

    python3 tools/scratch_copy.py <rev> <dir>

Extracts `git archive <rev>` into `<dir>` (missing or empty, outside this
repository and outside the main checkout that holds it), copies the built
extensions (`*.so`) of each package in the copy's source roots (`src_python`,
`src`) from this worktree's `.venv` into the copy's package, and prints one
line: the command that runs `pytest` on the copy with the worktree's
interpreter. To run some tests only, replace the `tests/python/` it ends with
by their paths.

The command sets `PYTHONPATH` to a generated `sitecustomize.py` in
`<dir>/.scratch_copy/` and to the copy's source roots. The sitecustomize drops
the editable finder at every interpreter start, after the `.pth` file installed
it, and the worktree's source roots a `.pth` file put on `sys.path`, so the
packages are imported from the copy in pytest's process and in any
Python process a test starts. It cannot reach a child whose environment
replaces `PYTHONPATH`, such as one built from scratch
(`tests/python/test_io_geotiff.py:1275`): that child imports the worktree's
code, and a mutant run there can report a false survivor.

Spec: docs/increments/h16-harness-fixes.md §2, T2. A refusal or a git error
is exit 2 with one `scratch_copy:` line on stderr.
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
#: As in count_loc.py: no user or system config, no replace refs.
GIT_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
           "GIT_NO_REPLACE_OBJECTS": "1"}  # fmt: skip
#: The C++ an extension is built from; the private `src/` holds Python only.
CXX = ("*.cpp", "*.hpp", "*.h", "*CMakeLists.txt", ":(exclude)tests/")
SOURCE_ROOTS = ("src_python", "src")
#: Run at every interpreter start by way of `PYTHONPATH`, after a `DROP = [...]` line naming
#: the worktree's source roots: drop the editable finders and those `sys.path` entries only
#: (dropping everything under the worktree would drop the venv's site-packages too).
SITECUSTOMIZE = (
    "import os, sys\n"
    'sys.meta_path[:] = [f for f in sys.meta_path if "editable" not in repr(f).lower()]\n'
    "sys.path[:] = [p for p in sys.path if not p or os.path.realpath(p) not in DROP]\n"
)
PROGRAM = "import sys, pytest; sys.exit(pytest.main(sys.argv[1:]))"


def _main_checkout() -> Path:
    """The work tree of the repository's common git dir: the main checkout."""
    common = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--path-format=absolute",
                             "--git-common-dir"], capture_output=True, text=True, env=GIT_ENV,
                            check=False)  # fmt: skip
    return Path(common.stdout.strip()).resolve().parent if common.returncode == 0 else ROOT


def _extract(archive: bytes, target: Path) -> None:
    """`tar -x` into `target`; on failure, leave it missing or empty, as it was."""
    existed = target.exists()
    target.mkdir(parents=True, exist_ok=True)
    tar = subprocess.run(["tar", "-x", "-C", str(target)], input=archive, capture_output=True,
                         check=False)  # fmt: skip
    if tar.returncode == 0:
        return
    if existed:
        for child in target.iterdir():
            shutil.rmtree(child) if child.is_dir() and not child.is_symlink() else child.unlink()
    else:
        shutil.rmtree(target)
    first = (tar.stderr.decode(errors="replace").strip().splitlines() or ["failed"])[0]
    raise RuntimeError(f"tar -x: {first}")


def copy(rev: str, target: Path) -> str:
    """Make the copy; return the command line to print. `RuntimeError` refuses."""
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise RuntimeError(f"{target} exists and is not empty")
    for repo in (ROOT, _main_checkout()):
        if target.is_relative_to(repo):
            raise RuntimeError(f"{target} is inside the repository {repo}")
    archive = subprocess.run(["git", "-C", str(ROOT), "archive", "--end-of-options", rev],
                             capture_output=True, env=GIT_ENV, check=False)  # fmt: skip
    if archive.returncode != 0:
        first = (archive.stderr.decode(errors="replace").strip().splitlines() or ["failed"])[0]
        raise RuntimeError(f"git archive {rev}: {first}")
    _extract(archive.stdout, target)
    roots = [target / name for name in SOURCE_ROOTS if (target / name).is_dir()]
    packages = [p for r in roots for p in sorted(r.iterdir()) if (p / "__init__.py").is_file()]
    built = [(so, package / so.relative_to(installed)) for package in packages
             for installed in sorted(ROOT.glob(f".venv/lib/python3.*/site-packages/{package.name}"))
             for so in sorted(installed.rglob("*.so"))]  # fmt: skip
    if not built:
        names = ", ".join(p.name for p in packages) or "no package"
        print(f"scratch_copy: no built extension of {names} in {ROOT}/.venv; suites needing it "
              "will fail", file=sys.stderr)  # fmt: skip
    for so, place in built:
        place.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(so, place)
    same = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", rev, "HEAD", "--", *CXX],
                          capture_output=True, env=GIT_ENV, check=False)  # fmt: skip
    if built and same.returncode != 0:
        print(f"scratch_copy: the copied extensions were built from HEAD, whose C++ differs from "
              f"{rev}", file=sys.stderr)  # fmt: skip
    site = target / ".scratch_copy"
    site.mkdir()
    drop = [str((ROOT / name).resolve()) for name in SOURCE_ROOTS]
    (site / "sitecustomize.py").write_text(f"DROP = {drop!r}\n{SITECUSTOMIZE}")
    path = os.pathsep.join(map(str, [site, *roots]))
    python = ROOT / ".venv" / "bin" / "python"
    quoted = (shlex.quote(str(target)), shlex.quote(path), shlex.quote(str(python)),
              shlex.quote(PROGRAM))  # fmt: skip
    return "cd {} && PYTHONPATH={} {} -c {} tests/python/".format(*quoted)


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("scratch_copy: usage: scratch_copy.py <rev> <dir>", file=sys.stderr)
        return 2
    try:
        print(copy(args[0], Path(args[1]).resolve()))
    except RuntimeError as exc:
        print(f"scratch_copy: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
