#!/usr/bin/env python3
"""Rebuild this worktree's `_core` and install it where its venv imports it.

    python3 tools/rebuild_core.py

Runs `cmake --build <root>/build-pyext -j --target _core`, then copies
`build-pyext/_core<suffix>` over `tin_engine/_core<suffix>` in the
site-packages of `<root>/.venv`, where `<suffix>` is that interpreter's
extension suffix. The installed file is located from the venv's own
`sysconfig` paths, never by importing it: the stale file may not load.
Prints one line, the destination and the first 12 hex digits of its SHA-256.

Exit 2 with one `rebuild_core:` line when `build-pyext` is not configured
(`tools/new_worktree.py` configures it), or the built or installed file is
missing; a failed build returns cmake's own exit status, and nothing is
copied. Standard library only.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
#: Asked of the venv's interpreter: its extension suffix, then its site-packages.
QUERY = (
    "import sysconfig; p = sysconfig.get_paths(); "
    "print(sysconfig.get_config_var('EXT_SUFFIX')); print(p['platlib']); print(p['purelib'])"
)


def fail(message: str, code: int = 2) -> int:
    print(f"rebuild_core: {message}", file=sys.stderr)
    return code


def main() -> int:
    build = ROOT / "build-pyext"
    if not (build / "CMakeCache.txt").is_file():
        return fail(f"{build} is not configured; python3 tools/new_worktree.py --existing {ROOT}")
    status = subprocess.run(
        ["cmake", "--build", str(build), "-j", "--target", "_core"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    ).returncode  # fmt: skip
    if status != 0:
        return fail(f"cmake --build exited {status}; nothing copied", status)
    python = ROOT / ".venv" / "bin" / "python"
    try:
        asked = subprocess.run([str(python), "-I", "-c", QUERY], capture_output=True,
                               text=True, check=True, timeout=60).stdout.split("\n")  # fmt: skip
    except (OSError, subprocess.SubprocessError) as error:
        return fail(f"cannot ask {python} for its site-packages: {error}")
    suffix, *sites = asked[:3]
    built = build / f"_core{suffix}"
    if not built.is_file():
        return fail(f"no {built.name} in {build} for the venv's suffix {suffix}; nothing copied")
    found = [Path(s) / "tin_engine" / built.name for s in sites]
    installed = [p for p in found if p.is_file()]
    if not installed:
        return fail(f"the venv has no installed {found[0]}; nothing copied")
    shutil.copyfile(built, installed[0])
    digest = hashlib.sha256(installed[0].read_bytes()).hexdigest()[:12]
    print(f"installed {installed[0]} {digest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
