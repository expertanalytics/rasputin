"""Master's CLI with prototype fixes patched in, to time them and compare outputs.

1. every buffer of the same polygon with the same arguments computed once (memo);
2. feature_input.source_region from the domain's convex hull grown by 100 m,
   instead of the convex hull of the domain grown by 100 m;
3. a mitre buffer of a long outline (5000 vertices or more, no holes) grown in
   pieces of 1000 edges, each overlapping the next by one edge, then united with
   the polygon.

Patch 3 replaces the method patch 1 wraps, so "123" runs as "23".
Layout as run.py. Usage: python run_patched.py <any of 1, 2, 3> -- <rasputin args...>
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in type(f).__module__]
sys.path.insert(0, str(HERE / "pkg"))

import shapely  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402
from shapely.geometry.base import BaseGeometry  # noqa: E402

from tin_engine import feature_input  # noqa: E402

which, args = sys.argv[1], sys.argv[3:]

if "1" in which:
    original = BaseGeometry.buffer
    memo: dict = {}

    def buffer(self, distance, *a, **k):
        key = (shapely.to_wkb(self), distance, a, tuple(sorted(k.items())))
        if key not in memo:
            memo[key] = original(self, distance, *a, **k)
        return memo[key]

    BaseGeometry.buffer = buffer

if "2" in which:

    def source_region(domain, dem_crs, source_crs):
        hull = shapely.convex_hull(domain.polygon)
        ring = shapely.segmentize(hull.buffer(feature_input.MARGIN).exterior, feature_input.DENSIFY)
        xy = shapely.get_coordinates(ring)
        if not feature_input.same_crs(source_crs, dem_crs):
            xy = feature_input.reprojector(dem_crs, source_crs)(xy)
        out = shapely.convex_hull(shapely.multipoints(xy))
        assert isinstance(out, Polygon)
        return out

    feature_input.source_region = source_region

if "3" in which:
    import numpy as np

    plain = BaseGeometry.buffer

    def piecewise(self, distance, *a, **k):
        mitre = k.get("join_style") == "mitre"
        if not (isinstance(self, Polygon) and not self.interiors and mitre):
            return plain(self, distance, *a, **k)
        xy = shapely.get_coordinates(self.exterior)
        n = len(xy) - 1
        if n < 5000:
            return plain(self, distance, *a, **k)
        size = 1000
        lines = []
        for s in range(0, n, size):
            stop = s + size + 2
            piece = xy[s:stop] if stop <= n + 1 else np.vstack([xy[s:], xy[1 : stop - n]])
            lines.append(shapely.linestrings(piece))
        grown = shapely.buffer(np.array(lines), distance, join_style="mitre", cap_style="flat")
        return shapely.union_all(np.concatenate([[self], grown]))

    BaseGeometry.buffer = piecewise

from tin_engine.cli import app  # noqa: E402

sys.argv = ["rasputin", *args]
app()
