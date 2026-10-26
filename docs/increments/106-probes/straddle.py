"""Fix 3 probe: in a mesh .gpkg, the triangles that straddle two units of a
units .geojson (both EPSG:3006), and the cells each unit gets when every
triangle takes the unit that holds its centroid.
    .venv/bin/python docs/increments/106-probes/straddle.py MESH.gpkg UNITS.geojson [SHARE]
SHARE (default 1e-6) is the share of a triangle in a second unit that counts.
"""

import json
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import shape

from tin_engine.feature_input import read_source as read_rows

mesh, units_path = Path(sys.argv[1]), Path(sys.argv[2])
limit = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
rows = read_rows(mesh, "triangles", None, lambda _c: (-1e12, -1e12, 1e12, 1e12)).rows
tris = np.array([shapely.force_2d(g) for _, g, _ in rows], dtype=object)
feats = json.loads(units_path.read_text())["features"]
units = np.array([shape(f["geometry"]) for f in feats], dtype=object)
t, u = shapely.STRtree(units).query(tris, predicate="intersects")
share = shapely.area(shapely.intersection(tris[t], units[u])) / shapely.area(tris[t])
real = share > limit
per = np.bincount(t[real], minlength=len(tris))
straddle = int((per > 1).sum())
cent = shapely.centroid(tris)
ct, cu = shapely.STRtree(units).query(cent, predicate="within")
owner = np.full(len(tris), -1)
owner[ct] = cu
cells = np.bincount(owner[owner >= 0], minlength=len(units))
print(
    f"{mesh.name}: {len(tris)} triangles; straddling two units (share > {limit:g}), "
    "largest share outside its main unit "
    f"{np.sort(np.bincount(t, weights=share * (share < 0.5)))[-1]:.1e}: {straddle}; "
    f"centroid in no unit: {int((owner < 0).sum())}; cells per unit min {cells.min()}, "
    f"median {int(np.median(cells))}, max {cells.max()}"
)
