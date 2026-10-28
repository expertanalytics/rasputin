#!/usr/bin/env python3
"""Refuse a branch that changes a full-lane path without an increment file (h34 §3).

    python3 tools/check_lane.py <base> <head>

Lists `git diff --name-only --no-renames <base>...<head>`; a path starting
with a `FULL_LANE` entry needs a `docs/increments/<name>.md` (one level,
`README.md` excluded) in the same diff. Exit 1 names each such path; exit 2
on a git failure. Spec: docs/increments/h34-lean-process.md §3. Standard
library only.
"""

from __future__ import annotations

import re
import subprocess
import sys

#: The meshing core (Ola's words), then elevation sampling, hydrology,
#: catchments and map projection (the main session's defaults, his call open).
FULL_LANE = (
    "include/terrain/cdt/", "include/terrain/predicates/", "include/terrain/refinement/",
    "include/terrain/noding/", "include/terrain/core/", "include/terrain/mesh/", "src/",
    "include/terrain/raster/", "src_python/tin_engine/elevation.py",
    "src_python/tin_engine/raster.py", "src_python/tin_engine/dem_input.py",
    "src_python/tin_engine/mosaic.py", "src_python/tin_engine/target_grid.py",
    "include/terrain/hydrology/", "src_python/tin_engine/catchment.py",
    "src_python/tin_engine/catchment_core.py", "src_python/tin_engine/network.py",
    "src_python/tin_engine/sub_catchments.py", "src_python/tin_engine/crs.py",
    "src_python/tin_engine/geodesy.py",
)  # fmt: skip
INCREMENT = re.compile(r"docs/increments/(?!README\.md$)[^/]+\.md")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("check_lane: usage: check_lane.py <base> <head>", file=sys.stderr)
        return 2
    base, head = argv
    done = subprocess.run(["git", "diff", "--name-only", "--no-renames", f"{base}...{head}"],
                          capture_output=True, text=True, check=False)  # fmt: skip
    if done.returncode != 0:
        print(f"check_lane: git diff {base}...{head}: {done.stderr.strip()}", file=sys.stderr)
        return 2
    paths = done.stdout.splitlines()
    full = [p for p in paths if p.startswith(FULL_LANE)]
    if not full or any(INCREMENT.fullmatch(p) for p in paths):
        return 0
    print("check_lane: full-lane paths changed without a docs/increments/<name>.md file:")
    for path in full:
        print(f"  {path}")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
