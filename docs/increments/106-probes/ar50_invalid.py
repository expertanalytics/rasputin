"""Fix 1 probe: AR50 Innlandet polygons meeting the Nedre Sjodalsvatn catchment
(NVE gauge 2.13.0): how many are invalid, whether the polygon cut throws on
them, and what shapely.make_valid changes. Run with the worktree's venv:
    .venv/bin/python docs/increments/106-probes/ar50_invalid.py
"""

from pathlib import Path

import numpy as np
import shapely

from tin_engine.feature_input import read_class_polygons, read_source
from tin_engine.legends import load as load_legend
from tin_engine.polygon_fractions import assign

AR50 = Path("/Users/skavhaug/projects/rasputin_data/AR50_UTM32_Innlandet_20261019/ar50.gpkg")
GAUGES = Path.home() / ".cache/rasputin/nve-gauge-catchments/features.gpkg"

rows = read_source(GAUGES, None, "stasjonnr", lambda crs: (-1e9, -1e9, 1e9, 1e9)).rows
(outline,) = [g for _, g, v in rows if str(v) == "2.13.0"]
print(f"outline 2.13.0: {outline.area / 1e6:.2f} km2 in EPSG:25833")
xy = shapely.get_coordinates(outline)
legend = load_legend("ar50")
found = read_class_polygons(AR50, None, legend, "EPSG:25833", xy)
bad = [(c, g) for c, g in found if not g.is_valid]
print(f"AR50 polygons read: {len(found)}, invalid: {len(bad)}")
for c, g in bad:
    fixed = shapely.make_valid(g)
    polys = shapely.get_parts(fixed)
    kept = [p for p in polys if p.geom_type == "Polygon"]
    print(
        f"  class {c}: {shapely.is_valid_reason(g)}; area {g.area:.1f} m2 -> "
        f"{sum(p.area for p in kept):.1f} m2 in {len(kept)} polygon(s), "
        f"{len(polys) - len(kept)} non-polygon part(s) dropped"
    )
# one triangle per invalid polygon, its box's lower-left half: does the cut throw?
for c, g in bad:
    x0, y0, x1, y1 = g.bounds
    v = np.array([[x0, y0, 0], [x1, y0, 0], [x0, y1, 0]], float)
    try:
        assign(v, [[0, 1, 2]], [(c, g)], known=frozenset({c}))
        print(f"  class {c}: cut ran")
    except shapely.errors.GEOSException as exc:
        print(f"  class {c}: cut threw {type(exc).__name__}: {str(exc)[:80]}")
