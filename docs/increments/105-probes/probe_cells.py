"""Probe for docs/increments/105-cell-table.md: the cell table Kultsjön would get.

Reads a coarse mesh written by `rasputin mesh --triangles N --out x.gpkg`
(its `triangles` layer), cuts every triangle against the CORINE 2018 polygons
in the clc2018 fetch cache (EPSG:3035, moved to the mesh's EPSG:3006 with
pyproj), and sums the cut areas into Shyft's five land-type fractions. The
reservoir is the CORINE water polygon that holds the SMHI dam register's
"Kultsjön" dam point (the delivery's reference-dams.csv), standing in for a
reservoir polygon file. Scratch, not product code: the product route is the
design's section 4.

    .venv/bin/python docs/increments/105-probes/probe_cells.py MESH.gpkg
"""

from __future__ import annotations

import sqlite3
import sys
import time
from collections import defaultdict

import numpy as np
import shapely
from pyproj import Transformer
from shapely import STRtree

from tin_engine.io.geopackage import decode_geometry

CLC = "/Users/skavhaug/.cache/rasputin/clc2018/features.gpkg"
FOREST, LAKE, GLACIER = {311, 312, 313}, {511, 512, 521, 522, 523}, {335}
DAM_XY = (520615.8, 7204484.4)  # SMHI dam "Kultsjön", EPSG:3006 (reference-dams.csv)


def triangles(path: str) -> tuple[list[shapely.Polygon], np.ndarray]:
    conn = sqlite3.connect(path)
    rows = conn.execute("SELECT fid, geom FROM triangles ORDER BY fid").fetchall()
    polys = [decode_geometry(blob)[1] for _, blob in rows]
    xyz = np.array([shapely.get_coordinates(p, include_z=True)[:3] for p in polys])
    return polys, xyz


def corine(bounds: tuple[float, float, float, float]) -> tuple[list[shapely.Geometry], list[int]]:
    to_3035 = Transformer.from_crs("EPSG:3006", "EPSG:3035", always_xy=True)
    to_3006 = Transformer.from_crs("EPSG:3035", "EPSG:3006", always_xy=True)
    x0, y0, x1, y1 = bounds
    xs, ys = to_3035.transform([x0, x1, x0, x1], [y0, y0, y1, y1])
    conn = sqlite3.connect(CLC)
    rows = conn.execute(
        "SELECT f.geom, f.Code_18 FROM features f JOIN rtree_features_geom r ON f.OBJECTID = r.id "
        "WHERE r.maxx >= ? AND r.minx <= ? AND r.maxy >= ? AND r.miny <= ?",
        (min(xs), max(xs), min(ys), max(ys)),
    ).fetchall()
    geoms = [
        shapely.transform(
            decode_geometry(b)[1], lambda c: np.column_stack(to_3006.transform(c[:, 0], c[:, 1]))
        )
        for b, _ in rows
    ]
    return geoms, [int(c) for _, c in rows]


def main(mesh: str) -> None:
    t0 = time.perf_counter()
    polys, xyz = triangles(mesh)
    geoms, codes = corine(shapely.total_bounds(polys))
    reservoir = {
        i
        for i, g in enumerate(geoms)
        if codes[i] in LAKE and g.distance(shapely.Point(DAM_XY)) < 200
    }
    tree = STRtree(geoms)
    n = len(polys)
    frac = np.zeros((n, 5))  # glacier, lake, reservoir, forest, unspecified
    by_code: dict[int, float] = defaultdict(float)
    area = shapely.area(polys)
    for t, tri in enumerate(polys):
        for i in tree.query(tri, predicate="intersects"):
            a = shapely.intersection(tri, geoms[i]).area
            if a <= 0:
                continue
            c = codes[i]
            by_code[c] += a
            k = (
                0
                if c in GLACIER
                else 2
                if i in reservoir
                else 1
                if c in LAKE
                else 3
                if c in FOREST
                else 4
            )
            frac[t, k] += a / area[t]
    frac[:, 4] += np.maximum(0.0, 1.0 - frac.sum(axis=1))  # uncovered slivers: unspecified
    e1, e2 = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
    normal = np.cross(e1, e2)
    normal[normal[:, 2] < 0] *= -1
    slope = np.degrees(np.arctan2(np.hypot(normal[:, 0], normal[:, 1]), normal[:, 2]))
    aspect = np.degrees(np.arctan2(normal[:, 0], normal[:, 1])) % 360.0
    mid = xyz.mean(axis=1)
    seconds = time.perf_counter() - t0
    print(
        f"{mesh}: {n} triangles, {len(geoms)} CORINE polygons in the box, "
        f"{len(reservoir)} reservoir polygon(s), {seconds:.1f} s"
    )
    print("cell  x_m      y_m       z_m    area_km2 slope aspect glacier lake  resv  forest unspec")
    for t in [*range(8), n - 1]:
        print(
            f"{t + 1:4d} {mid[t, 0]:9.0f} {mid[t, 1]:9.0f} {mid[t, 2]:7.1f} {area[t] / 1e6:8.3f} "
            f"{slope[t]:5.1f} {aspect[t]:6.1f} " + " ".join(f"{f:5.3f}" for f in frac[t])
        )
    total = area.sum()
    print(f"totals over {total / 1e6:.1f} km2 (area-weighted):")
    for k, name in enumerate(("glacier", "lake", "reservoir", "forest", "unspecified")):
        km2 = (frac[:, k] * area).sum() / 1e6
        print(f"  {name:12s} {km2:8.2f} km2  {100 * km2 * 1e6 / total:5.2f} %")
    pure = (frac.max(axis=1) >= 0.999).sum()
    print(
        f"triangles of one category (>= 0.999): {pure} of {n}; row sums in "
        f"[{frac.sum(axis=1).min():.6f}, {frac.sum(axis=1).max():.6f}]"
    )
    print(
        "by CORINE code, km2:",
        {c: round(a / 1e6, 2) for c, a in sorted(by_code.items(), key=lambda x: -x[1])},
    )
    cut = sum(by_code.values()) / 1e6
    print(f"area of triangles {total / 1e6:.3f} km2; area of CORINE cut {cut:.3f} km2")


if __name__ == "__main__":
    main(sys.argv[1])
