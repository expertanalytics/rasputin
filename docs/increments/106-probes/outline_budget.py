"""Fix 5 probe: an official outline simplified with the band the budget sets.
Each Opal outline as SMHI's network joins it (catchments-original/, EPSG:3006,
unsimplified) goes alone to 106c's `sub_catchments.simplify` (the outline and
its frame) at band = h / 8, h the edge of an equilateral triangle of area
A / N. Prints, per outline and budget N: vertices before and after, the band,
the area change, validity, holes, seconds; writes OUT/NAME-N.geojson (EPSG:3006)
for `mesh --domain`.
Run:  .venv/bin/python -I docs/increments/106-probes/outline_budget.py OUT N [N ...]
"""

import json
import sys
import time
from pathlib import Path

from shapely.geometry import mapping, shape

from tin_engine import sub_catchments

ROOT = Path.home() / "rasputin-output/opal-se2-2026-10-21/catchments-original"
out, budgets = Path(sys.argv[1]), [int(n) for n in sys.argv[2:]]
out.mkdir(parents=True, exist_ok=True)
for path in sorted(ROOT.glob("*.geojson")):
    given = shape(json.loads(path.read_text())["features"][0]["geometry"])
    n0, a0 = len(given.exterior.coords) - 1, given.area
    for n in budgets:
        band = sub_catchments.band_m(a0, n)
        t = time.perf_counter()
        (kept,) = sub_catchments.simplify([given], band)
        s = time.perf_counter() - t
        print(
            f"{path.stem:14} A {a0 / 1e6:9.2f} km2 N {n:6} band {band:7.1f} m  vertices "
            f"{n0:7} -> {len(kept.exterior.coords) - 1:5}  area change {kept.area - a0:+.3g} m2  "
            f"valid {kept.is_valid} holes {len(kept.interiors)} type {kept.geom_type} {s:.2f} s"
        )
        doc = {
            "type": "FeatureCollection",
            "crs": {"type": "name", "properties": {"name": "EPSG:3006"}},
            "features": [{"type": "Feature", "properties": {}, "geometry": mapping(kept)}],
        }
        (out / f"{path.stem}-{n}.geojson").write_text(json.dumps(doc))
