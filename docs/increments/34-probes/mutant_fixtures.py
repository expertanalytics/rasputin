"""Increment 34 design probe: the layouts that let tests 4 (b) and 6 kill M6 and
M3, with measured figures (README item 4b; design review round 1, B3).

M6 (the foot's epsilon from the triangle part's tolerance F instead of the
node's N). On V1 with a domain whose south-east edge crosses the wall
diagonally, off-node, rasputin master meshes at --tolerance 2 with the
constraint feet on (the default). Every output vertex that is a DEM node lies
at least eps(2) from the constraint, or it would have gone in as its foot;
the ones closer than eps(10) to the edge would have been footed under M6, so
the mesh would differ. eps(tol) = clamp(tol / G, dx / 100, dx / 2), G the
largest bilinear slope bound over the cells around the node
(refine.hpp, foot_epsilon).

The line of tests 5 and 7: how many nodes each driver holds tighter.

M3 (a check point's class from its nearest node instead of its cell's largest
corner). Check points as test_core_refine_points.scattered(65, 4, 7) places
them, in its draw order: 4 per cell at dyadic offsets, z the bilinear surface
plus noise of up to 4 m. A point is exposed when its cell's corners straddle
the step at 30 degrees and its nearest corner is below it: M3 holds it to
F = 10 where the rule holds it to N = 2.

Usage: python mutant_fixtures.py PYTHON RASPUTIN_PY PKG
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "tests" / "python"))
from fixture_figures import valley  # noqa: E402
from slope_stats import horn  # noqa: E402

from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff  # noqa: E402

PY, RASPUTIN, PKG = sys.argv[1], sys.argv[2], sys.argv[3]
DX = 10.0
# The node rectangle, its south-east corner cut by an edge across the wall.
A, B = (640.0, -300.0), (150.0, -640.0)  # metres from (TIE_X, TIE_Y), off-node
RING = [(0.0, -640.0), (B[0], B[1]), (A[0], A[1]), (640.0, 0.0), (0.0, 0.0)]


def foot_eps(z, r, c, tol):
    g = 0.0
    for rr in (r - 1, r):
        for cc in (c - 1, c):
            if 0 <= rr < z.shape[0] - 1 and 0 <= cc < z.shape[1] - 1:
                z00, z01, z10, z11 = z[rr, cc], z[rr, cc + 1], z[rr + 1, cc], z[rr + 1, cc + 1]
                gx = max(abs(z01 - z00), abs(z11 - z10)) / DX
                gy = max(abs(z10 - z00), abs(z11 - z01)) / DX
                g = max(g, float(np.hypot(gx, gy)))
    cap = DX / 2
    return min(max(tol / g, DX / 100), cap) if g > 0 else cap


def seg_dist(p, a, b):
    p, a, b = np.asarray(p), np.asarray(a), np.asarray(b)
    t = np.clip(np.dot(p - a, b - a) / np.dot(b - a, b - a), 0.0, 1.0)
    return float(np.hypot(*(p - (a + t * (b - a)))))


def m6(z):
    with tempfile.TemporaryDirectory() as tmp:
        tif = Path(tmp) / "v1.tif"
        tif.write_bytes(micro_tiff(z.astype(np.float32), scale=(DX, DX, 0.0)).getvalue())
        ring = [[TIE_X + x, TIE_Y + y] for x, y in RING]
        dom = Path(tmp) / "cut.geojson"
        crs = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::25833"}}
        dom.write_text(json.dumps({"type": "Polygon", "coordinates": [ring], "crs": crs}))
        ply, rec = Path(tmp) / "m.ply", Path(tmp) / "m.json"
        subprocess.run(
            [PY, RASPUTIN, "mesh", "--dem", str(tif), "--domain", str(dom), "--tolerance", "2",
             "--out", str(ply), "--ascii", "--record", str(rec)],
            env={"PYTHONPATH": PKG}, check=True, capture_output=True,
        )  # fmt: skip
        lines = ply.read_text().splitlines()
        n = int(next(x for x in lines if x.startswith("element vertex")).split()[-1])
        start = lines.index("end_header") + 1
        xy = np.array([[float(v) for v in x.split()[:2]] for x in lines[start : start + n]])
        record = json.loads(rec.read_text())
    feet = record.get("points_snapped_to_lines")
    on_wall = exposed = 0
    for x, y in xy - [TIE_X, TIE_Y]:
        c, r = x / DX, -y / DX
        if c != round(c) or r != round(r):
            continue  # not a DEM node: a foot or an outline vertex off the lattice
        r, c = round(r), round(c)
        d = seg_dist((x, y), A, B)
        if 200.0 <= x <= 440.0 and d < DX / 2:
            on_wall += 1
            exposed += foot_eps(z, r, c, 2.0) <= d < foot_eps(z, r, c, 10.0)
    print(
        f"M6: --tolerance 2 on V1 cut by the edge {A} to {B}: {len(xy)} vertices, "
        f"{feet} feet; DEM-node vertices on the wall within {DX / 2:g} m of the edge: "
        f"{on_wall}, of them between eps(2) and eps(10) from it (footed under M6): {exposed}"
    )


def m3(z, s):
    """The check points in tests/python/test_core_refine_points.py's
    scattered(65, 4, 7) draw order: all column offsets, then all row offsets,
    then all noise, over the cells in row-major order, 4 per cell."""
    cls = np.ceil(2.0 * s)
    step = 60  # 30 degrees in half degrees
    n, per_cell = z.shape[0], 4
    rng = np.random.default_rng(7)
    cells = (n - 1) * (n - 1) * per_cell
    base_r, base_c = np.divmod(np.repeat(np.arange((n - 1) * (n - 1)), per_cell), n - 1)
    col = base_c + rng.integers(1, 1024, cells) / 1024.0
    row = base_r + rng.integers(1, 1024, cells) / 1024.0
    noise = rng.integers(-256, 257, cells) / 64.0
    lo = np.minimum.reduce([cls[:-1, :-1], cls[:-1, 1:], cls[1:, :-1], cls[1:, 1:]])
    hi = np.maximum.reduce([cls[:-1, :-1], cls[:-1, 1:], cls[1:, :-1], cls[1:, 1:]])
    mixed_cell = (hi >= step) & (lo < step)
    in_mixed = mixed_cell[base_r, base_c]
    near = cls[np.round(row).astype(int), np.round(col).astype(int)] < step
    exposed = in_mixed & near
    print(
        f"M3: V1 cells whose corners straddle 30 deg: {int(mixed_cell.sum())}; check points "
        f"(scattered(65, 4, 7)) in them nearest a corner below 30 deg: {int(exposed.sum())}; "
        f"of those with noise above N = 2 m: {int((exposed & (np.abs(noise) > 2.0)).sum())}"
    )


def line(s):
    """V1's line for tests 5 and 7: (100, -50) to (600, -600) m from the
    north-west node, N 1 m, ramp 0 to 200 m, with the slope's N 2 from 30 deg."""
    rows, cols = s.shape
    a, b = (100.0, -50.0), (600.0, -600.0)
    d = np.array([[seg_dist((c * DX, -r * DX), a, b) for c in range(cols)] for r in range(rows)])
    t_line = np.where(d >= 200.0, 10.0, 1.0 + 9.0 * d / 200.0)
    t_slope = np.where(np.ceil(2.0 * s) >= 60, 2.0, 10.0)
    print(
        f"line: V1 nodes within 200 m of it {np.sum(d < 200.0)} of {d.size}; held tighter by "
        f"the line than by the slope {np.sum(t_line < t_slope)}, by the slope than by the "
        f"line {np.sum(t_slope < t_line)}"
    )


def main():
    z = valley(65, 65, DX, DX).astype(np.float64)
    s = horn(z, DX, DX)
    m6(z)
    m3(z, s)
    line(s)


if __name__ == "__main__":
    main()
