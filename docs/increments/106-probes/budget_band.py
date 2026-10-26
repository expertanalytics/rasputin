"""Fix 3 probe, round 3: the border band taken from the triangle budget.
SMHI's SVAR units whose inner point lies in the Opal outline NAME, joined by
BARONR (SMHI's main sub-basins), repaired, then simplified with the frame at
band = FACTOR * h, h = sqrt(4 A / (sqrt(3) N)) the edge of an equilateral
triangle of area A / N (A the units' union area). Prints, per budget N, the
band, the border vertices (distinct points of the simplified sub-basins), the
border length and length / h (the triangle edges a border of that length
meets), and writes OUT/sub-N.geojson (the sub-basins, for --features) and
OUT/dom-N.geojson (their union, for --domain).
Run:  .venv/bin/python -I docs/increments/106-probes/budget_band.py OUT NAME FACTOR N [N ...]
"""

import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import mapping, shape

from tin_engine.border_simplify import simplify_borders
from tin_engine.io.network_file import read_source

SVAR = Path(
    "/Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/wfs-20261021/"
    "SVAR2022_delavrinningsomraden.gpkg"
)
out, name, factor = Path(sys.argv[1]), sys.argv[2], float(sys.argv[3])
budgets = [int(n) for n in sys.argv[4:]]
opal = shape(
    json.loads(
        (Path.home() / f"rasputin-output/opal-se2-2026-10-21/catchments/{name}.geojson").read_text()
    )["features"][0]["geometry"]
)  # EPSG:3006, as SVAR
rows = read_source(SVAR, None, "ARO_UUID", lambda _c: opal.buffer(5000).bounds, None, ("BARONR",))
groups: dict[int, list] = defaultdict(list)
for (_, g, _v), (b,) in zip(rows.rows, rows.extra, strict=True):
    if opal.contains(g.point_on_surface()):
        groups[int(b)].append(g)
units = shapely.coverage_clean(
    np.array([g for v in groups.values() for g in v], dtype=object),
    snapping_distance=1,
    gap_width=1,
)
keys = [k for k, v in groups.items() for _ in v]
subs = np.array(
    [shapely.coverage_union_all(units[[i for i, k in enumerate(keys) if k == b]]) for b in groups],
    dtype=object,
)
union = shapely.coverage_union_all(subs)
area = union.area
x0, y0, x1, y1 = union.bounds
a0 = shapely.area(subs)


def border_stats(polys: np.ndarray) -> tuple[int, float]:
    pts = np.concatenate([shapely.get_coordinates(p) for p in polys])
    length = (sum(shapely.length(shapely.boundary(polys))) + shapely.length(union.boundary)) / 2
    return len(np.unique(pts, axis=0)), float(length)


v0, l0 = border_stats(subs)
print(
    f"{name}: {len(units)} units in {len(subs)} main sub-basins (BARONR {sorted(groups)}), "
    f"{area / 1e6:.2f} km2, source border vertices {v0}, border length {l0 / 1e3:.1f} km"
)
for n in budgets:
    h = math.sqrt(4 * area / (math.sqrt(3) * n))
    band = factor * h
    frame = shapely.box(x0 - 2 * band - 1000, y0 - 2 * band - 1000, x1 + 2 * band + 1000,
                        y1 + 2 * band + 1000).difference(union)  # fmt: skip
    simple = np.array(simplify_borders([*subs, frame], band, 1.0).polygons[:-1], dtype=object)
    dom = shapely.coverage_union_all(simple)
    v, length = border_stats(simple)
    print(
        f"  N {n}: h {h:.0f} m, band {band:.0f} m, border vertices {v}, "
        f"length {length / 1e3:.1f} km, "
        f"length/h {length / h:.0f}, mean segment / h {length / v / h:.2f}, "
        f"largest area change {np.max(np.abs(shapely.area(simple) - a0) / a0):.1e}, "
        f"union {dom.geom_type} holes "
        f"{len(dom.interiors) if dom.geom_type == 'Polygon' else '-'} area change "
        f"{dom.area - area:+.3f} m2, coverage valid {shapely.coverage_is_valid(simple)}"
    )
    fc = {"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "EPSG:3006"}}}
    feats = [
        {"type": "Feature", "id": b, "properties": {"property": "wall"}, "geometry": mapping(g)}
        for b, g in zip(groups, simple, strict=True)
    ]
    (out / f"sub-{n}.geojson").write_text(json.dumps({**fc, "features": feats}))
    dom_fc = {**fc, "features": [{"type": "Feature", "properties": {}, "geometry": mapping(dom)}]}
    (out / f"dom-{n}.geojson").write_text(json.dumps(dom_fc))
    (out / f"src-{name}.geojson").write_text(json.dumps({**fc, "features": [
        {"type": "Feature", "id": b, "properties": {}, "geometry": mapping(g)}
        for b, g in zip(groups, subs, strict=True)]}))  # fmt: skip
