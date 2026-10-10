"""Increment 34 design probe: the measured figures of the test fixtures that
section 9 specifies (docs/increments/34-probes/README.md, item 4).

V1 is the valley of section 9: a flat floor, a 40-degree wall, a plateau, with
bumps so that refinement has work everywhere. V1 is 65 x 65 nodes at 10 m (the
C++ suites and the binding tests); V2 is the same valley on the CLI fixtures'
micro_tiff grid, 33 rows x 41 cols at 10 m (x) by 5 m (y).

For each: the share of nodes per steepness class band (Horn, half-degree
classes rounded up, section 3), and the triangle counts of rasputin master's
uniform meshes at F = 10 and N = 2, from the outline of the node rectangle,
and the Python simulation's counts for the rules (greedy_sim.py).
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TREE = HERE.parents[2]
sys.path.insert(0, str(TREE / "tests" / "python"))
from greedy_sim import run, tol_of  # noqa: E402
from slope_stats import horn  # noqa: E402

from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff  # noqa: E402

PY, RASPUTIN, PKG = sys.argv[1], sys.argv[2], sys.argv[3]


def valley(rows: int, cols: int, dx: float, dy: float) -> np.ndarray:
    """z in metres: floor (x < 200 m), a 40-degree wall to x = 440 m, a plateau;
    plus 3 m bumps, sin(2 pi y / 80 m) sin(2 pi x / 130 m)."""
    r, c = np.indices((rows, cols), dtype=np.float64)
    x, y = c * dx, r * dy
    base = np.tan(np.radians(40.0)) * np.clip(x - 200.0, 0.0, 240.0)
    return (base + 3.0 * np.sin(2 * np.pi * y / 80.0) * np.sin(2 * np.pi * x / 130.0)).astype(
        np.float32
    )


def counts(z: np.ndarray, dx: float, dy: float, tol: str) -> str:
    rows, cols = z.shape
    with tempfile.TemporaryDirectory() as tmp:
        tif = Path(tmp) / "valley.tif"
        tif.write_bytes(micro_tiff(z, scale=(dx, dy, 0.0)).getvalue())
        x1, y1 = TIE_X + (cols - 1) * dx, TIE_Y - (rows - 1) * dy
        ring = [(TIE_X, y1), (x1, y1), (x1, TIE_Y), (TIE_X, TIE_Y), (TIE_X, y1)]
        dom = Path(tmp) / "box.geojson"
        dom.write_text(
            json.dumps(
                {
                    "type": "Polygon",
                    "coordinates": [ring],
                    "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::25833"}},
                }
            )
        )
        out = subprocess.run(
            [
                PY,
                RASPUTIN,
                "mesh",
                "--dem",
                str(tif),
                "--domain",
                str(dom),
                "--tolerance",
                tol,
                "--out",
                str(Path(tmp) / "m.vtk"),
            ],
            env={"PYTHONPATH": PKG},
            capture_output=True,
            text=True,
        )
        return next(
            (
                line.split(".")[0]
                for line in out.stdout.splitlines() + out.stderr.splitlines()
                if "triangles." in line
            ),
            out.stderr[-300:],
        )


def main() -> None:
    for name, (rows, cols, dx, dy) in {
        "V1": (65, 65, 10.0, 10.0),
        "V2": (33, 41, 10.0, 5.0),
    }.items():
        z = valley(rows, cols, dx, dy)
        s = horn(z, dx, dy)
        cls = np.ceil(2.0 * s) / 2.0
        print(f"{name}: {rows} x {cols} nodes, {dx:g} x {dy:g} m, z {z.min():.1f}..{z.max():.1f} m")
        print(
            "  classes: "
            + " ".join(f">={a}: {100 * np.mean(cls >= a):.1f}%" for a in (15, 25, 30, 35, 40))
        )
        print(
            f"  nodes held to N=2 by step 30: {np.sum(tol_of(cls, 2.0, 10.0, 30.0, 30.0) == 2.0)} "
            f"of {z.size}; by ramp 25..35, below F: "
            f"{np.sum(tol_of(cls, 2.0, 10.0, 25.0, 35.0) < 10.0)}"
        )
        for tol in ("10", "2"):
            print(f"  rasputin uniform {tol}: {counts(z, dx, dy, tol)}")
        if dx == dy:
            for rule in ("uniform-F", "uniform-N", "node", "tri+halo"):
                ntri, rounds, over, _ = run(z.astype(np.float64), s, rule, 2.0, 10.0, 30.0, 30.0)
                print(
                    f"  simulation {rule:9s} {ntri} triangles, {rounds} rounds, "
                    f"nodes over {over:.3f}%"
                )


if __name__ == "__main__":
    main()
