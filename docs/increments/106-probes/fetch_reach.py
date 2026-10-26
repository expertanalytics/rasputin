"""Fix 5 probe, round 2: how far past the box `fetch --domain D` covers does
mesh's outline reach, with `--simplify-domain` (D simplified at the budget's
band) and with `--catchments` (the simplified union of the SVAR units whose
inner point lies in D)? `fetch` covers D's bounds grown by its margin (about
186 m for glo30, `fetch/plan.py`, `source_box`). Per Opal outline (EPSG:3006,
the terrain-traced Kultsjön moved there from EPSG:4326) at N = 500: the
largest of the four amounts by which each outline's bounds pass D's, in m.
Run:  .venv/bin/python -I docs/increments/106-probes/fetch_reach.py
"""

import json
from pathlib import Path

import shapely
from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import transform

from tin_engine import sub_catchments
from tin_engine.io.network_file import read_source

SVAR = Path(
    "/Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/wfs-20261021/"
    "SVAR2022_delavrinningsomraden.gpkg"
)
OPAL = Path.home() / "rasputin-output/opal-se2-2026-10-21"
N = 500
to_3006 = Transformer.from_crs("EPSG:4326", "EPSG:3006", always_xy=True).transform


def reach(outline: shapely.Geometry, d: shapely.Geometry) -> float:
    (a0, b0, a1, b1), (x0, y0, x1, y1) = outline.bounds, d.bounds
    return max(x0 - a0, y0 - b0, a1 - x1, b1 - y1)


cases = [(p.stem, p, False) for p in sorted((OPAL / "catchments-original").glob("*.geojson"))]
cases.append(("kultsjon-terrain", OPAL / "catchments/kultsjon-lake-seeded-check.geojson", True))
for name, path, degrees in cases:
    d = shape(json.loads(path.read_text())["features"][0]["geometry"])
    d = transform(to_3006, d) if degrees else d
    (flag,) = sub_catchments.simplify([d], sub_catchments.band_m(d.area, N))
    box = d.buffer(20000).bounds
    rows = read_source(SVAR, None, "ARO_UUID", lambda _crs, box=box: box).rows
    near = [g for _, g, _ in rows if g is not None]
    units = sub_catchments.repair([near[k] for k in sub_catchments.select(near, d)])
    band = sub_catchments.band_m(shapely.coverage_union_all(units).area, N)
    union = shapely.coverage_union_all(sub_catchments.simplify(units, band))
    print(
        f"{name:17} --simplify-domain reach {reach(flag, d):7.1f} m   "
        f"--catchments {len(units):4} units reach {reach(union, d):7.1f} m"
    )
