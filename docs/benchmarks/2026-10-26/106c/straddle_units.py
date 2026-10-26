"""106c acceptance: in a `mesh --catchments` GeoPackage, how far any cell
reaches outside the unit holding most of it, and whether each cell's
catchment_id is that unit's id.

The units are rebuilt by the CLI's own `_units` (select, repair, simplify at
the band the run record gives), so they are the polygons the mesh was
constrained by, in their own CRS (SVAR: EPSG:3006, the mesh's --out-crs).
Run from the worktree's venv, with the package of the measured head:
    .venv/bin/python docs/benchmarks/2026-10-26/106c/straddle_units.py \
        MESH.gpkg UNITS.gpkg DOMAIN.geojson BAND_M
"""

import sqlite3
import sys
from pathlib import Path

import numpy as np
import shapely

from tin_engine import cli
from tin_engine.feature_input import read_source
from tin_engine.io.domain_file import read_domain

mesh, units_path, domain_path, band = (
    Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), float(sys.argv[4]),
)  # fmt: skip
units = cli._units(units_path, None, units_path.stem, read_domain(domain_path), None, band)
ids = np.array([i for i, _ in units.units])
polys = np.array([g for _, g in units.units], dtype=object)

layers = [r[0] for r in sqlite3.connect(mesh).execute("select table_name from gpkg_contents")]
box = (-1e12, -1e12, 1e12, 1e12)
rows = read_source(mesh, "triangles", None, lambda _c: box).rows
tris = np.array([shapely.force_2d(g) for _, g, _ in rows], dtype=object)
cells = read_source(mesh, "cells", "catchment_id", lambda _c: box).rows  # a point per cell
cid = np.array([int(v) for _, _, v in cells])
points = np.array([shapely.force_2d(g) for _, g, _ in cells], dtype=object)
assert len(points) == len(tris) and shapely.covers(tris, points).all(), "cells not in row order"

t, u = shapely.STRtree(polys).query(tris, predicate="intersects")
share = shapely.area(shapely.intersection(tris[t], polys[u])) / shapely.area(tris[t])
main = np.full(len(tris), -1)
best = np.zeros(len(tris))
for k in np.argsort(share):  # the last write per triangle is its largest share
    main[t[k]], best[t[k]] = u[k], share[k]
outside = 1.0 - best
assert (main >= 0).all(), "a cell meets no unit"
wrong = int((ids[main] != cid).sum())
per = np.bincount(main, minlength=len(polys))
print(
    f"{mesh.name}: {len(tris)} cells, {len(polys)} units, band {band:.3f} m; "
    f"largest share of a cell outside the unit holding most of it {outside.max():.2e} "
    f"(cells over 1e-4: {int((outside > 1e-4).sum())}); "
    f"catchment_id not that unit's id: {wrong}; distinct ids {len(set(cid.tolist()))}; "
    f"cells per unit min {per.min()}, median {int(np.median(per))}, max {per.max()}; "
    f"layers {layers}"
)
