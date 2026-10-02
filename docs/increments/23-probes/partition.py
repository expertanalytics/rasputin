"""Increment 23's arithmetic for the piece partition of Ola's 2026-10-01 ruling.

    python docs/increments/23-probes/partition.py DATA

DATA is `../rasputin_data/sao_francisco_piece` (see its README.md). For the
Velhas piece (its 30 m EPSG:31983 grid as `basin-piece/prep_dem.py` built it)
and for the basin (BHO level 2, used here only as the basin's extent), on a
30 m lattice in EPSG:31983: the window in nodes, then the partition rule of
`23-basin-scale.md` ("Choosing the cuts") at the default `--memory-budget`,
for several tolerances and `--pieces` requests, with
the cells that meet the domain, the share of the domain's area in the largest
cell against the mean, and the seam length inside the domain. A measurement
script, not production code; nothing imports it.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import pyproj
import shapely
from shapely.geometry import shape
from shapely.ops import transform

H_M = 30.0
BUDGET = 16 << 30  # bytes; the default --memory-budget
# Bytes per window node at a tolerance (m), from the basin-piece sweep: 17 B per
# node plus 310 B per triangle times the Velhas piece's triangles per domain node;
# at 0 every node is a vertex (two triangles per node).
BYTES = ((0.0, 637), (1.0, 267), (2.0, 155), (5.0, 65), (10.0, 37), (20.0, 25), (50.0, 19))


def bytes_per_node(tol: float) -> int:
    """Linear between the table's points, rounded up; the last value beyond it."""
    for (t0, b0), (t1, b1) in itertools.pairwise(BYTES):
        if tol <= t1:
            return math.ceil(b0 + (b1 - b0) * (tol - t0) / (t1 - t0))
    return BYTES[-1][1]


def partition(cols: int, rows: int, tol: float, pieces: int) -> tuple[int, int, int, int]:
    """(nx, ny, dx, dy): the rule as designed, in whole nodes."""
    b = bytes_per_node(tol)
    p = max(pieces, -(-cols * rows * b // BUDGET))
    if p <= 1:
        return 1, 1, cols, rows
    nx = max(1, round(math.sqrt(p * cols / rows)))
    ny = max(1, round(p / nx))
    dx, dy = math.ceil(cols / nx), math.ceil(rows / ny)
    while dx * dy * b > BUDGET:
        if dx >= dy:
            nx += 1
        else:
            ny += 1
        dx, dy = math.ceil(cols / nx), math.ceil(rows / ny)
    return math.ceil(cols / dx), math.ceil(rows / dy), dx, dy


def report(name: str, dom: shapely.Geometry, x0: float, y1: float, cols: int, rows: int) -> None:
    print(
        f"\n{name}: window {cols} x {rows} = {cols * rows / 1e6:.1f} M nodes, "
        f"domain {dom.area / 1e6:,.0f} km2, "
        f"cover {dom.area / ((cols - 1) * (rows - 1) * H_M * H_M):.0%}"
    )
    shapely.prepare(dom)
    for tol, pieces in ((50, 1), (10, 1), (5, 1), (1, 1), (0.5, 1), (0, 1), (1, 16), (1, 64)):
        nx, ny, dx, dy = partition(cols, rows, tol, pieces)
        areas, seam = [], 0.0
        for j in range(ny):
            for i in range(nx):
                box = shapely.box(
                    x0 + i * dx * H_M,
                    y1 - (j + 1) * dy * H_M,
                    x0 + (i + 1) * dx * H_M,
                    y1 - j * dy * H_M,
                )
                if dom.intersects(box):
                    areas.append(dom.intersection(box).area)
        lines = [
            shapely.LineString([(x0 + i * dx * H_M, y1), (x0 + i * dx * H_M, y1 - rows * H_M)])
            for i in range(1, nx)
        ]
        lines += [
            shapely.LineString([(x0, y1 - j * dy * H_M), (x0 + cols * H_M, y1 - j * dy * H_M)])
            for j in range(1, ny)
        ]
        seam = sum(dom.intersection(ln).length for ln in lines)
        mean = sum(areas) / len(areas)
        print(
            f"  {tol:>4} m, {bytes_per_node(tol):>3} B/node, --pieces {pieces:>2}: "
            f"{nx:>3} x {ny:>3} cells of {dx} x {dy} nodes "
            f"({dx * dy / 1e6:.2f} M), {len(areas):>4} meet the domain, "
            f"largest/mean area {max(areas) / mean:.2f}, seams {seam / 1e3:,.0f} km"
        )


def main(data: Path) -> None:
    to_utm = pyproj.Transformer.from_crs("EPSG:4674", "EPSG:31983", always_xy=True).transform
    piece = transform(
        to_utm,
        shape(
            json.loads((data / "bho2017_5k_76949_outline_epsg4674.geojson").read_text())[
                "features"
            ][0]["geometry"]
        ),
    )
    meta = json.loads((data / "derived" / "bho2017_5k_76949_glo30_epsg31983_30m.json").read_text())
    rows, cols = meta["shape"]  # numpy order
    x0, y1 = meta["first_node"]
    report("Velhas piece", piece, x0, y1, cols, rows)
    basin = transform(
        to_utm,
        shapely.union_all(
            [
                shape(f["geometry"])
                for f in json.loads((data / "bho2017_level2_76_raw.geojson").read_text())[
                    "features"
                ]
            ]
        ),
    )
    bx0, by0, bx1, by1 = basin.bounds
    x0, y1 = math.floor(bx0 / H_M) * H_M - H_M, math.ceil(by1 / H_M) * H_M + H_M
    cols = round((math.ceil(bx1 / H_M) * H_M + H_M - x0) / H_M) + 1
    rows = round((y1 - (math.floor(by0 / H_M) * H_M - H_M)) / H_M) + 1
    report("Basin", basin, x0, y1, cols, rows)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
