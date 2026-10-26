"""Fix 3 probe: Kultsjön's 75 SVAR 2022 units simplified as one coverage.
Needs the joined outline from
    rasputin network-catchment --out OUT/kultsjon-joined.geojson --network <SVAR gpkg>
      --id-attribute ARO_UUID --down-attribute MAINDOWN
      --unit 11B6DEC8-0A79-40FE-AF91-7413AE9B4703
Run:  .venv/bin/python docs/increments/106-probes/units_simplify.py OUT [OPAL-NAME]
Writes OUT/units-<m>.geojson (the units, for --features) and
OUT/domain-<m>.geojson (their union, for --domain), m in metres.
"""

import json
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import mapping, shape

from tin_engine.io.network_file import read_source, unit_geometries

DATA = Path("/Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/wfs-20261021")
SVAR = DATA / "SVAR2022_delavrinningsomraden.gpkg"
NAME = sys.argv[2] if len(sys.argv) > 2 else "kultsjon"
OPAL = Path.home() / f"rasputin-output/opal-se2-2026-10-21/catchments/{NAME}.geojson"
out = Path(sys.argv[1])
joined = json.loads((out / "kultsjon-joined.geojson").read_text())["features"][0]
ids = joined["properties"]["units"]
units = unit_geometries(SVAR, None, "ARO_UUID", ids)
print(
    f"units: {len(units)}, vertices {sum(shapely.get_num_coordinates(g) for g in units.values())}"
)

# selection rule: a unit belongs when its point on surface is inside the domain
opal = shape(json.loads(OPAL.read_text())["features"][0]["geometry"])
box = opal.buffer(5000).bounds
near = read_source(SVAR, None, "ARO_UUID", lambda _crs: box).rows
picked = {str(v) for _, g, v in near if g is not None and opal.contains(g.point_on_surface())}
print(
    f"units near the Opal outline: {len(near)}; point on surface inside it: {len(picked)}; "
    f"same set as the join: {picked == set(ids)}"
)

order = sorted(units)
geoms = np.array([units[k] for k in order], dtype=object)
bad = shapely.coverage_invalid_edges(geoms)
print(
    f"as published: coverage valid {shapely.coverage_is_valid(geoms)}; units with invalid "
    f"edges {int((~shapely.is_empty(bad)).sum())}, {shapely.length(bad).sum():.1f} m of them"
)
area0 = shapely.area(geoms)
geoms = shapely.coverage_clean(geoms)
print(
    f"after coverage_clean: valid {shapely.coverage_is_valid(geoms)}, "
    f"area change {(shapely.area(geoms).sum() - area0.sum()):+.1f} m2"
)
area0 = shapely.area(geoms)
for tol in (0.0, 50.0, 100.0, 200.0, 300.0):
    simple = geoms if tol == 0 else shapely.coverage_simplify(geoms, tol)
    d = (shapely.area(simple) - area0) / area0
    union = shapely.coverage_union_all(simple)
    valid = shapely.coverage_is_valid(simple)
    change = np.abs(d)
    total = (shapely.area(simple).sum() - area0.sum()) / 1e6
    holes = len(union.interiors) if union.geom_type == "Polygon" else "-"
    print(
        f"tol {tol:5.0f} m: vertices {int(shapely.get_num_coordinates(simple).sum()):7d}; "
        f"unit area change max {change.max() * 100:.2f} %, median {np.median(change) * 100:.3f} %; "
        f"units over 5 % {int((change > 0.05).sum())}; total {total:+.3f} km2; "
        f"coverage valid {valid}; union {union.geom_type}, holes {holes}, "
        f"outline vertices {shapely.get_num_coordinates(union)}; "
        f"smallest unit {area0.min() / 1e6:.3f} km2"
    )
    if tol in (100.0, 200.0, 300.0):
        crs = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::3006"}}
        feats = [
            {
                "type": "Feature",
                "properties": {"property": "wall", "unit": k},
                "geometry": mapping(g),
            }
            for k, g in zip(order, simple, strict=True)
        ]
        (out / f"units-{tol:.0f}.geojson").write_text(
            json.dumps({"type": "FeatureCollection", "crs": crs, "features": feats})
        )
        (out / f"domain-{tol:.0f}.geojson").write_text(
            json.dumps(
                {
                    "type": "FeatureCollection",
                    "crs": crs,
                    "features": [{"type": "Feature", "properties": {}, "geometry": mapping(union)}],
                }
            )
        )
