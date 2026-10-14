"""Increment 46 PR 6 acceptance, the one profile: horizon-search time on
Bygdin against the search distance, to see whether the pruning shortens the
rays (time flat in distance) or not (time growing with it).

    python search_reach.py <mesh.vtk> <dem path> <receivers> <out.json>

A fixed random subset of the mesh's receivers (seed 46), one thread, the
near raster read once for 20 km and searched to 2.5, 5, 10 and 20 km.
"""

import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np

from tin_engine.io.vtk_read import read_mesh
from tin_engine.radiation.frame import receivers_from_mesh
from tin_engine.radiation.request import RadiationRequest
from tin_engine.radiation.search import near_raster, search_horizons

mesh, dem, m, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
a = read_mesh(Path(mesh).read_bytes())
v = a.vertices
pick = np.sort(np.random.default_rng(46).choice(len(a.triangles), m, replace=False))
frame = receivers_from_mesh(v, a.triangles[pick], a.crs)
box = (*v[:, :2].min(axis=0), *v[:, :2].max(axis=0))
near = near_raster((Path(dem),), box, RadiationRequest())
rows = []
for d in (2500.0, 5000.0, 10000.0, 20000.0):
    req = RadiationRequest(d_near_m=d)
    c = time.perf_counter()
    g = search_horizons(frame, near, req, threads=1)
    s = time.perf_counter() - c
    rows.append({"d_near_m": d, "seconds": s, "ms_per_receiver": s / m * 1e3,
                 "us_per_ray": s / (m * req.n_directions) * 1e6,
                 "median_horizon_deg": float(np.median(g.horizon) / 100)})
    print(rows[-1], flush=True)
Path(out).write_text(json.dumps({"receivers": m, "rows": rows}, indent=2) + "\n")
