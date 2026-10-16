"""Q1 for 15f-3's acceptance (@perf): the mesh's constraint-line vertices by
where they sit on the DEM (or target) lattice: on a node, on one grid line
(a crossing), or off every grid line. Comparing the counts for the base mesh and
for 15f-3's gives how many crossings and how many midpoints the strip inserted.
Does not import tin_engine.
python q1_vertices.py SPACING MESH.vtk [MESH.vtk ...]  (lattice: x/h, -y/h; nodes on multiples of h)"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from strip_check import read_vtk  # noqa: E402

h = float(sys.argv[1])
for f in sys.argv[2:]:
    pts, edges = read_vtk(Path(f))
    v = np.unique(edges)
    c, r = pts[v, 0] / h, -pts[v, 1] / h
    oc = np.abs(c - np.round(c)) < 1e-6
    orow = np.abs(r - np.round(r)) < 1e-6
    print(f"{f}: constraint vertices {len(v)}, on nodes {int((oc & orow).sum())}, "
          f"on one grid line {int((oc ^ orow).sum())}, off the grid lines {int((~oc & ~orow).sum())}")
