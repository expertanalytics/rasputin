"""Every constraint segment that only one piece has must lie on the domain's
outline (@perf). Midpoints of those segments against the outline projected
into EPSG:31983 (vertices projected, straight between, as the mesh's outline).
python once_on_outline.py PIECES_DIR OUTLINE.geojson"""
import sys
from pathlib import Path

import numpy as np
import pyproj
import shapely
import shapely.ops

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/docs/benchmarks/2026-10-04/15f-3-acceptance")
from strip_check import read_vtk  # noqa: E402  (no tin_engine)

pieces, outline = Path(sys.argv[1]), Path(sys.argv[2])
dom = shapely.from_geojson(outline.read_text())
if dom.geom_type == "GeometryCollection":
    dom = shapely.union_all(list(shapely.get_parts(dom)))
tr = pyproj.Transformer.from_crs("EPSG:4674", "EPSG:31983", always_xy=True)
ring = shapely.ops.transform(tr.transform, dom).boundary
count: dict = {}
for f in sorted(pieces.glob("*.vtk")):
    pts, lines = read_vtk(f)
    for a, b in lines:
        ka, kb = tuple(np.round(pts[a, :2], 6)), tuple(np.round(pts[b, :2], 6))
        k = (ka, kb) if ka <= kb else (kb, ka)
        count[k] = count.get(k, 0) + 1
once = [k for k, n in count.items() if n == 1]
mid = np.array([[(a[0] + b[0]) / 2, (a[1] + b[1]) / 2] for a, b in once])
d = shapely.distance(ring, shapely.points(mid))
print(f"segments in one piece only: {len(once)}; midpoint distance to the projected outline: "
      f"max {d.max():.6f} m; over 0.01 m: {int((d > 0.01).sum())}; over 1 m: {int((d > 1).sum())}")
print(f"segments in two pieces: {sum(1 for n in count.values() if n == 2)}; in more than two: "
      f"{sum(1 for n in count.values() if n > 2)}")
