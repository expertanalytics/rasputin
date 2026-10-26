"""Fix 3 probe, round 2: how far the Opal outline (the --domain) and the
union of the SVAR units selected by inner point differ, the area a cell
could lie in with no unit under it.
Run:  .venv/bin/python docs/increments/106-probes/opal_gaps.py NAME [NAME ...]
"""

import json
import sys
from pathlib import Path

import pyproj
import shapely
from shapely.geometry import shape
from shapely.ops import transform

from tin_engine.io.network_file import read_source

SVAR = Path(
    "/Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/wfs-20261021/"
    "SVAR2022_delavrinningsomraden.gpkg"
)
for name in sys.argv[1:]:
    doc = json.loads(
        (Path.home() / f"rasputin-output/opal-se2-2026-10-21/catchments/{name}.geojson").read_text()
    )
    opal = shape(doc["features"][0]["geometry"])
    crs = doc.get("crs", {}).get("properties", {}).get("name", "EPSG:4326")
    to3006 = pyproj.Transformer.from_crs(crs, "EPSG:3006", always_xy=True).transform
    d = transform(to3006, opal)
    box = d.buffer(5000).bounds
    near = read_source(SVAR, None, "ARO_UUID", lambda _crs, box=box: box).rows
    units = [g for _, g, _ in near if g is not None and d.contains(g.point_on_surface())]
    units = list(shapely.coverage_clean(units, snapping_distance=1, gap_width=1))
    u = shapely.coverage_union_all(units)
    gap, over = d.difference(u), u.difference(d)
    widths = [2 * p.area / p.length for p in shapely.get_parts(gap) if p.length > 0]
    print(
        f"{name}: crs {crs}, outline vertices {shapely.get_num_coordinates(d)}, "
        f"units {len(units)}; "
        f"outline {d.area / 1e6:.2f} km2, units' union {u.area / 1e6:.2f} km2; "
        f"in outline not in a unit {gap.area / 1e6:.3f} km2 ({100 * gap.area / d.area:.2f} %), "
        f"{len(widths)} pieces, widest mean width {max(widths, default=0):.1f} m; "
        f"in a unit outside the outline {over.area / 1e6:.3f} km2"
    )
    # The outline rule (feature_input.snap_to_outline) at distance s moves a
    # border within s of the outline onto it: a gap piece whose every vertex
    # is within s of the outline is closed by it.
    ring = d.exterior
    far = [
        (p.area, max(ring.distance(shapely.Point(c)) for c in shapely.get_coordinates(p)))
        for p in shapely.get_parts(gap)
    ]
    for s in (5.0, 50.0, 99.0):
        left = [a for a, m in far if m > s]
        print(f"  outline rule at {s:g} m: {len(left)} pieces left, {sum(left) / 1e6:.4f} km2")
    # Round 2's rule: the units and the frame around the outline cleaned as
    # one coverage, every gap given to a unit (min_index, the frame last),
    # then clipped to the outline; then the land-cover border simplification.
    import time

    import numpy as np

    from tin_engine.border_simplify import simplify_borders

    t0 = time.perf_counter()
    frame = shapely.box(*d.buffer(1000).bounds).difference(d)
    for width in (100.0, 1000.0):
        cleaned = shapely.coverage_clean(
            np.array([*units, frame], dtype=object),
            snapping_distance=1,
            gap_width=width,
            merge_strategy="min_index",
        )[:-1]
        tiles = shapely.intersection(cleaned, d)
        union = shapely.coverage_union_all(tiles)
        print(
            f"  gaps to units at {width:g} m: {time.perf_counter() - t0:.1f} s, coverage valid "
            f"{shapely.coverage_is_valid(tiles)}, union vs outline "
            f"{union.symmetric_difference(d).area:.1f} m2, "
            f"vertices {int(shapely.get_num_coordinates(tiles).sum())}"
        )
    for band in (50.0, 150.0, 300.0):
        t1 = time.perf_counter()
        out = simplify_borders(list(tiles), band, 1.0).polygons
        a0, a1 = shapely.area(tiles), shapely.area(np.array(out, dtype=object))
        print(
            f"  borders at {band:g} m: {time.perf_counter() - t1:.1f} s, vertices "
            f"{int(shapely.get_num_coordinates(np.array(out, dtype=object)).sum())}, "
            f"largest unit area change {np.max(np.abs(a1 - a0) / a0):.2e}, "
            f"valid {shapely.coverage_is_valid(np.array(out, dtype=object))}"
        )
