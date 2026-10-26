"""Fix 4 probe: SMHI's dam register joined to SMHI's lakes, inside Kultsjön.
.venv/bin/python docs/increments/106-probes/se_reservoirs.py
"""

import json
from pathlib import Path

import shapely

D = Path("/Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/wfs-20261021")
OUT = Path.home() / "rasputin-output/opal-se2-2026-10-21/catchments-original/kultsjon.geojson"

dams = json.loads((D / "dammar_3006.geojson").read_text())["features"]
print("dam fields:", sorted(dams[0]["properties"]))
keys = [d["properties"].get("vy_eucd") for d in dams]
named = {k for k in keys if k}
print(
    f"dams {len(dams)}, naming a water body (vy_eucd): {sum(1 for k in keys if k)}, "
    f"distinct {len(named)}"
)
print("examples:", sorted(named)[:3])
lakes = json.loads((D / "vattenytor_2016_3006.geojson").read_text())["features"]
outline = shapely.geometry.shape(json.loads(OUT.read_text())["features"][0]["geometry"])
stripped = {k.removeprefix("SE") for k in named}
for field in ("SJOID", "VYID"):
    hit = [f for f in lakes if f["properties"].get(field) in stripped]
    print(f"lakes whose {field} a dam names (SE stripped): {len(hit)}")
hit = [f for f in lakes if f["properties"].get("SJOID") in stripped]
inside = []
for f in hit:
    g = shapely.geometry.shape(f["geometry"])
    a = g.intersection(outline).area
    if a > 0:
        inside.append((f["properties"].get("NAME"), f["properties"]["SJOID"], a / 1e6))
for row in inside:
    print("  ", row)
print(f"in Kultsjon: {len(inside)} lakes, {sum(a for *_, a in inside):.2f} km2")
# the same with the prefix kept: does anything match without the strip?
print(
    "matches with the prefix kept:", sum(1 for f in lakes if f["properties"].get("SJOID") in named)
)
