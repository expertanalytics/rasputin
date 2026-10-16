"""How GEOS's overlay rewrites a polygon it does not cut: for every CORINE
polygon of Numedalslagen that the read region contains properly, compare
`intersection(g, region)` with g and with g re-oriented or normalised.
Run as .venv/bin/python docs/increments/35-probes/ring_order.py numedalslagen|lagan."""

from __future__ import annotations

import sys
from collections import Counter

import numpy as np
import shapely
from clean_probe import read_cover
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient

from tin_engine.feature_input import _polygonal


def main(case: str) -> None:
    _domain, _dem_crs, region, cover = read_cover(case)
    geoms = np.array([g for *_, g in cover], dtype=object)
    inside = shapely.contains_properly(region, geoms)
    seen: Counter[str] = Counter()
    shown = 0
    for g in geoms[inside]:
        cut = _polygonal(shapely.intersection(g, region))
        g = _polygonal(g)
        if shapely.to_wkb(cut) == shapely.to_wkb(g):
            seen["same"] += 1
            continue
        if g.geom_type != "Polygon" or cut.geom_type != "Polygon":
            seen[f"{g.geom_type}->{cut.geom_type}"] += 1
            continue
        def revrot(ring):
            xy = shapely.get_coordinates(ring)[:-1]
            return np.concatenate([xy[::-1], xy[-1:]])
        rr = g
        if g.exterior.is_ccw and all(not h.is_ccw for h in g.interiors):
            rr = Polygon(revrot(g.exterior), [revrot(h) for h in g.interiors])
        tries = (("revrot", rr), ("cw", orient(g, -1)), ("ccw", orient(g, 1)))
        for name, h in (*tries, ("normalize", shapely.normalize(g))):
            if shapely.to_wkb(cut) == shapely.to_wkb(h):
                seen[name] += 1
                break
        else:
            seen["other"] += 1
            if shown < 3:
                shown += 1
                a, b = shapely.get_coordinates(g.exterior), shapely.get_coordinates(cut.exterior)
                print(f"example: {len(a)} -> {len(b)} exterior points, {len(g.interiors)} -> "
                      f"{len(cut.interiors)} holes, ccw in {g.exterior.is_ccw} "
                      f"out {cut.exterior.is_ccw}, "
                      f"start in {a[0]} out {b[0]}, in starts at out index "
                      f"{np.flatnonzero((b == a[0]).all(axis=1))[:3]}, "
                      f"equals {shapely.equals(g, cut)}, "
                      f"holes ccw in {[h.is_ccw for h in g.interiors][:4]} "
                      f"out {[h.is_ccw for h in cut.interiors][:4]}")
    print(f"{case}: {int(inside.sum())} polygons inside the region: {dict(seen)}")


if __name__ == "__main__":
    main(sys.argv[1])
