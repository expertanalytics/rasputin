"""worst_angles.py A.vtk [...] -- the five smallest-angle triangles of each mesh:
angle, whether each corner is a DEM node (world coordinates on the 10 m grid of
the benchmark DEM), and whether a constraint edge bounds the triangle."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "tools"))
import bench  # noqa: E402

for path in sys.argv[1:]:
    m = bench.read_vtk_ascii(Path(path))
    xy, t = m.points[:, :2], m.triangles
    p = xy[t]
    ang = []
    for k in range(3):
        u, v = p[:, (k + 1) % 3] - p[:, k], p[:, (k + 2) % 3] - p[:, k]
        c = (u * v).sum(1) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)
        ang.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    a = np.min(ang, axis=0)
    cons = {frozenset(e) for e in m.edges.tolist()}
    x0, y0 = xy[:, 0].min(), xy[:, 1].max()
    print(Path(path).name)
    for i in np.argsort(a)[:5]:
        node = ["node" if abs(((xy[j, 0] - 0.0) / 10) % 1) < 1e-9 and abs((xy[j, 1] / 10) % 1) < 1e-9
                else "off" for j in t[i]]
        onc = any(frozenset((t[i][k], t[i][(k + 1) % 3])) in cons for k in range(3))
        print(f"  {a[i]:.4f} deg at ({xy[t[i]][:, 0].mean():.0f}, {xy[t[i]][:, 1].mean():.0f}) "
              f"corners {node} constraint edge {onc}")
