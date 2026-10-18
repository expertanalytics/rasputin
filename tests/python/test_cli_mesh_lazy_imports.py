"""LZ1: an uncut run does not import the cut path.

`docs/increments/23-basin-scale.md`, "23c-3: each piece's memory its own",
D8 and LZ1. In a fresh interpreter (a subprocess, so no other test's import
counts), `rasputin mesh` of the square scene with `--domain` and
`--tolerance 1` exits 0; that is one cell at the default `--piece-budget`,
so the run goes through the one-cell test, and after it `sys.modules` holds
none of the cut path's modules. A bare `import tin_engine.cli` would not
do: the run calls the partition.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import piece_spies as ps

CUT_PATH = ("tin_engine.basin_run", "tin_engine.pieces", "tin_engine.io.mesh_index")

#: The console script's own entry point, then what is loaded.
CHILD = """
import json, sys
from tin_engine.plugins import main
sys.argv = ["rasputin", *sys.argv[1:]]
try:
    main()
    code = 0
except SystemExit as exc:
    code = exc.code if isinstance(exc.code, int) else 1
loaded = sorted(m for m in sys.modules if m.startswith("tin_engine"))
print("LOADED " + json.dumps({"code": code, "modules": loaded}))
"""


def test_lz1_an_uncut_domain_tolerance_run_loads_none_of_the_cut_path(tmp_path: Path) -> None:
    args = ps.square_scene(tmp_path)
    out = tmp_path / "x.vtk"
    done = subprocess.run(
        [sys.executable, "-c", CHILD, "mesh", *args, "--out", str(out)],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
        cwd=tmp_path,
    )
    line = next((s for s in done.stdout.splitlines() if s.startswith("LOADED ")), None)
    assert line is not None, done.stdout + done.stderr
    found = json.loads(line.removeprefix("LOADED "))
    assert found["code"] == 0, done.stderr
    assert out.exists() and not (tmp_path / "x.vtk.pieces").exists(), "not an uncut run"
    assert "tin_engine.cli" in found["modules"], "the child did not run the CLI"
    assert [m for m in CUT_PATH if m in found["modules"]] == []
