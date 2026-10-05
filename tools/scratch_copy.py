#!/usr/bin/env python3
"""A copy of one revision to run the suites against, outside any work tree (h16 T2).

    python3 tools/scratch_copy.py <rev> <dir>

Extracts `git archive <rev>` into `<dir>` (missing or empty, outside this
repository), copies the built `_core` from this worktree's `.venv` into the
copy's package, and prints one line: the command that runs `pytest` on the
copy with the worktree's interpreter, the editable finder dropped so that
`tin_engine` is imported from the copy. Append test paths to it, or keep the
`tests/python/` it ends with.

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
CXX = ("include", "src", "bindings", "CMakeLists.txt")
#: Run by the worktree's interpreter: drop the editable finder, put the copy first.
PROGRAM = (
    "import sys; sys.meta_path[:] = "
    '[f for f in sys.meta_path if "editable" not in repr(f).lower()]; '
    'sys.path.insert(0, "src_python"); import pytest; sys.exit(pytest.main(sys.argv[1:]))'
)


def copy(rev: str, target: Path) -> str:
    """Make the copy; return the command line to print. `RuntimeError` refuses."""
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise RuntimeError(f"{target} exists and is not empty")
    if target == ROOT or target.is_relative_to(ROOT):
        raise RuntimeError(f"{target} is inside the repository {ROOT}")
    archive = subprocess.run(["git", "-C", str(ROOT), "archive", "--end-of-options", rev],
                             capture_output=True, env=GIT_ENV, check=False)  # fmt: skip
    if archive.returncode != 0:
        first = (archive.stderr.decode(errors="replace").strip().splitlines() or ["failed"])[0]
        raise RuntimeError(f"git archive {rev}: {first}")
    target.mkdir(parents=True, exist_ok=True)
    subprocess.run(["tar", "-x", "-C", str(target)], input=archive.stdout, check=True)
    cores = sorted(ROOT.glob(".venv/lib/python3.*/site-packages/tin_engine/_core*.so"))
    package = target / "src_python" / "tin_engine"
    if not cores:
        print(f"scratch_copy: no built _core in {ROOT}/.venv; suites needing it will fail",
              file=sys.stderr)  # fmt: skip
    for core in cores:
        package.mkdir(parents=True, exist_ok=True)
        shutil.copy2(core, package / core.name)
    same = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", rev, "HEAD", "--", *CXX],
                          capture_output=True, env=GIT_ENV, check=False)  # fmt: skip
    if cores and same.returncode != 0:
        print(f"scratch_copy: the copied _core was built from HEAD, whose C++ differs from {rev}",
              file=sys.stderr)  # fmt: skip
    python = ROOT / ".venv" / "bin" / "python"
    quoted = shlex.quote(str(target)), shlex.quote(str(python)), shlex.quote(PROGRAM)
    return "cd {} && {} -c {} tests/python/".format(*quoted)


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("scratch_copy: usage: scratch_copy.py <rev> <dir>", file=sys.stderr)
        return 2
    try:
        print(copy(args[0], Path(args[1]).resolve()))
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"scratch_copy: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
