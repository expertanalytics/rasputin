"""Probe for docs/increments/105-cell-table.md: what the two new GeoPackage
tables cost to write, beside the `triangles` layer they join to.

Writes a GeoPackage file in EPSG:3006 with N random triangles through the
product's `write_triangles`, then the design's `cells` layer (POINT Z, one row
per triangle, the columns of section 3, R-tree indexed) and the
`land_cover_fractions` attributes table (E rows, index on `triangle`), built
with the writer module's own blob helpers, and times each of the three. The
cells and fractions values are random: the cost is the SQLite insert and the
indexes, not the arithmetic. Scratch, not product code.

    .venv/bin/python docs/increments/105-probes/probe_write.py N E OUT.gpkg
"""

from __future__ import annotations

import os
import sqlite3
import sys
import time
from pathlib import Path

import numpy as np

from tin_engine.crs import srs_entry
from tin_engine.io.geopackage_write import (
    _records,
    _register_attributes,
    _split,
    build_rtree,
    create_geopackage,
    register_layer,
    write_triangles,
)

CELLS = (
    "CREATE TABLE cells (fid INTEGER PRIMARY KEY, geom POINT, catchment_id INTEGER, "
    "area_m2 REAL, z_m REAL, slope_deg REAL, aspect_deg REAL, glacier REAL, lake REAL, "
    "reservoir REAL, forest REAL, unspecified REAL)"
)
FRACTIONS = (
    "CREATE TABLE land_cover_fractions (fid INTEGER PRIMARY KEY, triangle INTEGER, "
    "code INTEGER, fraction REAL)"
)


def main(n: int, entries: int, out: str) -> None:
    rng = np.random.default_rng(105)
    x = rng.uniform(4.5e5, 5.0e5, 3 * n)
    y = rng.uniform(7.2e6, 7.25e6, 3 * n)
    vertices = np.column_stack([x, y, rng.uniform(300, 1500, 3 * n)])
    triangles = np.arange(3 * n).reshape(n, 3)
    srs = srs_entry("EPSG:3006")
    path = Path(out)
    path.unlink(missing_ok=True)
    conn = sqlite3.connect(path)
    create_geopackage(conn, srs)
    t0 = time.perf_counter()
    write_triangles(conn, vertices, triangles, None, None, srs.srs_id)
    conn.commit()
    t_tri = time.perf_counter() - t0

    t0 = time.perf_counter()
    mid = vertices[triangles].mean(axis=1)
    rec = _records(mid[:, None, :], srs.srs_id, 1001, [("xyz", "<f8", (3,))])
    rec["xyz"] = mid
    blobs, boxes = _split(rec)
    shares = rng.dirichlet(np.ones(5), n)
    columns = [
        np.ones(n, dtype=np.int64).tolist(),
        rng.uniform(1e3, 1e6, n).tolist(),
        mid[:, 2].tolist(),
        rng.uniform(0, 60, n).tolist(),
        rng.uniform(0, 360, n).tolist(),
        *(shares[:, k].tolist() for k in range(5)),
    ]
    conn.execute(CELLS)
    conn.execute("BEGIN")
    conn.executemany(
        "INSERT INTO cells VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        zip(range(1, n + 1), blobs, *columns, strict=True),
    )
    conn.commit()
    bounds = (boxes[:, 0].min(), boxes[:, 2].min(), boxes[:, 1].max(), boxes[:, 3].max())
    register_layer(conn, "cells", "POINT", srs.srs_id, bounds, z=True)
    build_rtree(conn, "cells", boxes)
    conn.commit()
    t_cells = time.perf_counter() - t0

    t0 = time.perf_counter()
    tri = np.sort(rng.integers(0, n, entries))
    code = rng.integers(1, 45, entries)
    frac = rng.uniform(0, 1, entries).astype(np.float32).astype(np.float64)
    conn.execute(FRACTIONS)
    conn.execute("BEGIN")
    conn.executemany(
        "INSERT INTO land_cover_fractions (triangle, code, fraction) VALUES (?,?,?)",
        zip(tri.tolist(), code.tolist(), frac.tolist(), strict=True),
    )
    conn.commit()
    conn.execute("CREATE INDEX land_cover_fractions_triangle ON land_cover_fractions (triangle)")
    _register_attributes(conn, "land_cover_fractions")
    conn.commit()
    t_frac = time.perf_counter() - t0
    conn.close()
    size = os.path.getsize(path)
    print(
        f"{n} triangles, {entries} fraction entries: triangles layer {t_tri:.2f} s; "
        f"cells layer {t_cells:.2f} s; land_cover_fractions table {t_frac:.2f} s; "
        f"file {size / 1e6:.1f} MB"
    )


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
