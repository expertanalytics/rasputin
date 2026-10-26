"""Fix 3 probe, round 2: the selected SVAR units and the frame around their
union simplified as one coverage by rasputin's own border simplification
(border_simplify.simplify_borders, the land-cover --features-tolerance
step), so the union's outline is simplified with the inner borders and no
gap can open. Writes OUT/units2-<band>.geojson (the units, for --features)
and OUT/domain2-<band>.geojson (their union, for --domain).
Run:  .venv/bin/python -I docs/increments/106-probes/units_borders.py OUT NAME [BAND ...]
"""

import json
import sys
import time
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
out, name = Path(sys.argv[1]), sys.argv[2]
bands = [float(b) for b in sys.argv[3:]] or [50.0, 100.0, 200.0, 300.0]
opal = shape(
    json.loads(
        (Path.home() / f"rasputin-output/opal-se2-2026-10-21/catchments/{name}.geojson").read_text()
    )["features"][0]["geometry"]
)  # EPSG:3006, as SVAR
t0 = time.perf_counter()
near = read_source(SVAR, None, "ARO_UUID", lambda _crs: opal.buffer(5000).bounds).rows
units = [g for _, g, _ in near if g is not None and opal.contains(g.point_on_surface())]
units = shapely.coverage_clean(np.array(units, dtype=object), snapping_distance=1, gap_width=1)
union = shapely.coverage_union_all(units)
x0, y0, x1, y1 = union.bounds  # bounds arithmetic: a buffer of the outline takes minutes
frame = shapely.box(x0 - 1000, y0 - 1000, x1 + 1000, y1 + 1000).difference(union)
print(
    f"{name}: {len(units)} units, {int(shapely.get_num_coordinates(units).sum())} vertices, "
    f"union {union.geom_type} {union.area / 1e6:.2f} km2, outline vertices "
    f"{shapely.get_num_coordinates(union)}, read, clean, union and frame "
    f"{time.perf_counter() - t0:.1f} s"
)
a0 = shapely.area(units)
for band in bands:
    t1 = time.perf_counter()
    simple = np.array(simplify_borders([*units, frame], band, 1.0).polygons[:-1], dtype=object)
    took = time.perf_counter() - t1
    dom = shapely.coverage_union_all(simple)
    a1 = shapely.area(simple)
    print(
        f"  band {band:g} m: {took:.1f} s, vertices "
        f"{int(shapely.get_num_coordinates(simple).sum())}, "
        f"outline {shapely.get_num_coordinates(dom)} ({dom.geom_type}, holes "
        f"{len(dom.interiors) if dom.geom_type == 'Polygon' else '-'}), coverage valid "
        f"{shapely.coverage_is_valid(simple)}, largest unit area change "
        f"{np.max(np.abs(a1 - a0) / a0):.1e}, union area change {(dom.area - union.area):+.3f} m2"
    )
    fc = {"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "EPSG:3006"}}}
    (out / f"units2-{band:g}.geojson").write_text(
        json.dumps(
            {
                **fc,
                "features": [
                    {
                        "type": "Feature",
                        "id": k + 1,
                        "properties": {"property": "wall", "unit": k + 1},
                        "geometry": mapping(g),
                    }
                    for k, g in enumerate(simple)
                ],
            }
        )
    )
    (out / f"domain2-{band:g}.geojson").write_text(
        json.dumps(
            {**fc, "features": [{"type": "Feature", "properties": {}, "geometry": mapping(dom)}]}
        )
    )
