"""Does the next stage wash out the overlay's ring order? C2 (no clip to the
read region for a polygon the region contains properly) feeds
`coverage_clean`; C6 (the clip to the domain part by part) feeds
`simplify_borders`. Each candidate's input differs from master's in ring order
only (ring_order.py); here the stage after it is run on both and compared as
WKB. Also C1': pyproj's transform of Lagan's 15.9 M vertices in chunks on
threads, one Transformer per thread, against the one-call result, bit for bit.
Run as .venv/bin/python docs/increments/35-probes/stage_check.py numedalslagen|lagan."""

from __future__ import annotations

import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import shapely
from clean_probe import ASK, CASES, read_cover
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Polygon

from tin_engine.border_simplify import simplify_borders
from tin_engine.feature_input import _outline_lines, _polygonal, snap_to_outline


def wkb(geoms):
    return [shapely.to_wkb(g) for g in geoms]


def main(case: str) -> None:
    domain, dem_crs, region, cover = read_cover(case)
    dom = domain.polygon
    shapely.prepare(dom)
    geoms = np.array([g for *_, g in cover], dtype=object)

    if CASES[case][3] is not None:  # C1': threads
        xy = shapely.get_coordinates(geoms)
        src, dst = CASES[case][3], dem_crs
        t = time.perf_counter()
        one = Transformer.from_crs(src, dst, always_xy=True)
        x, y = one.transform(xy[:, 0], xy[:, 1])
        t_one = time.perf_counter() - t
        for workers in (2, 4, 8):
            chunks = np.array_split(xy, workers * 4)
            def move(part, src=src, dst=dst):
                tr = Transformer.from_crs(src, dst, always_xy=True)
                return tr.transform(part[:, 0], part[:, 1])
            t = time.perf_counter()
            with ThreadPoolExecutor(workers) as pool:
                out = list(pool.map(move, chunks))
            t_thr = time.perf_counter() - t
            xs = np.concatenate([o[0] for o in out])
            ys = np.concatenate([o[1] for o in out])
            print(f"C1' {len(xy)} points: one call {t_one:.2f} s, {workers} threads "
                  f"{t_thr:.2f} s, same bits: {np.array_equal(xs, x) and np.array_equal(ys, y)}")

    # C2 then coverage_clean
    inside = shapely.contains_properly(region, geoms)
    master = [_polygonal(shapely.intersection(g, region)) for g in geoms]
    c2 = [
        _polygonal(g) if ok else _polygonal(shapely.intersection(g, region))
        for g, ok in zip(geoms, inside, strict=True)
    ]
    a = np.array([g for g in master if g is not None], dtype=object)
    b = np.array([g for g in c2 if g is not None], dtype=object)
    t = time.perf_counter()
    ca = shapely.coverage_clean(a, snapping_distance=1.0, gap_width=1.0, merge_strategy="min_area")
    t_a = time.perf_counter() - t
    cb = shapely.coverage_clean(b, snapping_distance=1.0, gap_width=1.0, merge_strategy="min_area")
    print(f"C2 after coverage_clean ({t_a:.2f} s): same WKB: {wkb(ca) == wkb(cb)}")

    # master's path on to the clip to the domain, then C6 beside it
    items = [c[:3] for c, g in zip(cover, master, strict=True) if g is not None]
    groups: dict[int, list[int]] = {}
    for i, (_, _, code) in enumerate(items):
        groups.setdefault(code, []).append(i)
    first = sorted(g[0] for g in groups.values())
    merged = [shapely.coverage_union_all(ca[groups[items[i][2]]]) for i in first]
    snapped = snap_to_outline(merged, dom, ASK["outline_snap_m"])
    polygons = snapped.polygons
    inside_dom = [
        Polygon() if (c := _polygonal(shapely.intersection(g, dom))) is None else c
        for g in polygons
    ]
    c6 = []
    for g in polygons:
        parts = shapely.get_parts(g)
        keep = shapely.contains_properly(dom, parts)
        cut = ~keep & shapely.intersects(dom, parts)
        rest = None
        if cut.any():
            rest = _polygonal(shapely.intersection(MultiPolygon(list(parts[cut])), dom))
        whole = [*parts[keep], *([] if rest is None else shapely.get_parts(rest))]
        c6.append(Polygon() if not whole else whole[0] if len(whole) == 1 else MultiPolygon(whole))
    sa = simplify_borders(inside_dom, ASK["tolerance_m"], ASK["repair_m"]).polygons
    sb = simplify_borders(c6, ASK["tolerance_m"], ASK["repair_m"]).polygons
    same_poly = wkb(sa) == wkb(sb)
    la = [wkb(_outline_lines(g, dom)) for g in sa]
    lb = [wkb(_outline_lines(g, dom)) for g in sb]
    print(f"C6 after simplify_borders: same polygons: {same_poly}; same lines: {la == lb}")


if __name__ == "__main__":
    main(sys.argv[1])
