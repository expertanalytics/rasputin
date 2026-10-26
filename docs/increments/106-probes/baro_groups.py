"""Fix 3 probe: the units of a joined catchment grouped by an attribute
(SVAR's BARONR, the main sub-basin number): how many groups, and whether each
group's union is one polygon.
    .venv/bin/python docs/increments/106-probes/baro_groups.py JOINED.geojson [ATTRIBUTE]
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import shapely

from tin_engine.feature_input import read_source

SVAR = Path(
    "/Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/wfs-20261021/SVAR2022_delavrinningsomraden.gpkg"
)
attr = sys.argv[2] if len(sys.argv) > 2 else "BARONR"
joined = json.loads(Path(sys.argv[1]).read_text())["features"][0]
ids = set(joined["properties"]["units"])
box = shapely.geometry.shape(joined["geometry"]).bounds
rows = read_source(SVAR, None, "ARO_UUID", lambda _c: box, None, (attr,))
groups: dict[object, list] = defaultdict(list)
for (_, g, v), (a,) in zip(rows.rows, rows.extra, strict=True):
    if str(v) in ids:
        groups[a].append(g)
clean = {
    k: shapely.coverage_union_all(shapely.coverage_clean(np.array(v, dtype=object)))
    for k, v in groups.items()
}
parts = {k: len(shapely.get_parts(u)) for k, u in clean.items()}
print(
    f"{len(ids)} units in {len(groups)} {attr} groups; units per group "
    f"{sorted(len(v) for v in groups.values())}; groups in more than one part: "
    f"{ {k: p for k, p in parts.items() if p > 1} }; areas km2 "
    f"{sorted(round(u.area / 1e6, 1) for u in clean.values())}"
)
