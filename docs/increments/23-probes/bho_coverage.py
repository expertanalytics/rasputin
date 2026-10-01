"""Increment 23's probe of BHO 2017 5k as a source of natural cuts.

    python docs/increments/23-probes/bho_coverage.py DATA

DATA is `../rasputin_data/sao_francisco_piece` (see its README.md for where the
files came from). Prints, for the 1,163 elementary catchments of ottobasin
76949: how often each polygon edge occurs (an exact coverage has every interior
edge exactly twice, bit for bit), the area of the union against the sum of the
areas, and, per Pfafstetter level (the ottocode prefix length), the number of
units, their sizes in EPSG:31983 and the length of the shared boundaries. Last,
the same seam length for a grid of 2048-node blocks at 30 m. A measurement
script, not production code; nothing imports it.
"""

from __future__ import annotations

import collections
import itertools
import json
import math
import sys
from pathlib import Path

import pyproj
import shapely
from shapely.geometry import shape
from shapely.ops import transform


def main(data: Path) -> None:
    raw = json.loads((data / "bho2017_5k_76949_elementary.geojson").read_text())
    feats = raw["features"]
    edges: collections.Counter[tuple[tuple[float, float], ...]] = collections.Counter()
    for f in feats:
        g = shape(f["geometry"])
        for ring in (g.exterior, *g.interiors):
            c = list(ring.coords)
            for a, b in itertools.pairwise(c):
                edges[tuple(sorted((a, b)))] += 1
    print("polygons", len(feats), "edge multiplicity", dict(collections.Counter(edges.values())))
    geos = [shape(f["geometry"]) for f in feats]
    union = shapely.union_all(geos)
    total = sum(g.area for g in geos)
    print("sum of areas / union area - 1:", total / union.area - 1.0)

    to_utm = pyproj.Transformer.from_crs(4674, 31983, always_xy=True).transform
    utm = [transform(to_utm, g) for g in geos]
    codes = [str(f["properties"]["COBACIA"]) for f in feats]
    outline = shapely.union_all(utm)
    print(f"piece {outline.area / 1e6:.1f} km2, outline {outline.exterior.length / 1e3:.1f} km")
    for level in (6, 7, 8):
        groups: dict[str, list[shapely.Geometry]] = collections.defaultdict(list)
        for code, g in zip(codes, utm, strict=True):
            groups[code[:level]].append(g)
        units = [shapely.union_all(v) for v in groups.values()]
        km2 = sorted(u.area / 1e6 for u in units)
        mnodes = max((b[2] - b[0]) / 30 * (b[3] - b[1]) / 30 for b in (u.bounds for u in units))
        seam = (sum(u.boundary.length for u in units) - outline.boundary.length) / 2
        print(
            f"level {level}: {len(units)} units, km2 min/median/max "
            f"{km2[0]:.0f}/{km2[len(km2) // 2]:.0f}/{km2[-1]:.0f}, "
            f"largest bbox {mnodes / 1e6:.1f} M nodes at 30 m, seams {seam / 1e3:.0f} km"
        )
    side = 2048 * 30.0
    x0, y0, x1, y1 = outline.bounds
    seam = 0.0
    for i in range(math.floor(x0 / side) + 1, math.ceil(x1 / side)):
        line = shapely.LineString([(i * side, y0 - 1), (i * side, y1 + 1)])
        seam += shapely.intersection(line, outline).length
    for j in range(math.floor(y0 / side) + 1, math.ceil(y1 / side)):
        line = shapely.LineString([(x0 - 1, j * side), (x1 + 1, j * side)])
        seam += shapely.intersection(line, outline).length
    print(f"2048-node blocks at 30 m: seams {seam / 1e3:.0f} km")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
