"""Fetch NVE's lakes along each river station's reach (@perf; the
input lake_above.py needs). Usage: lake_query.py DATA HERE OUT

`lakes.geojson` from `fetch-stations` holds only the lakes within ±100 m of
each station point (`fetch/nve.py`, `LAKE_ENVELOPE_HALF`), so a lake on the
reach farther up is missing from it. This queries NVE's lake layer 5
(`Innsjodatabase2/MapServer/5`) once per station seeded on the river or
placed on a river line (93 stations), with the envelope that holds both the station point ± 1000 m (the envelope of
`@architect`'s probe, "NVE's lakes" in the increment file) and the whole
mapped reach above P (placed as the batch does, `gauge.place`, the
defaults), each lake kept once by `objectid`. Writes OUT (GeoJSON, NVE's
data, not committed) and prints, per station, how far the reach reaches
past the ±1000 m box (0 when inside)."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

from shapely.geometry import LineString, box
from shapely.ops import substring

from tin_engine.crs import reprojector
from tin_engine.fetch.http import RangeClient, query_url
from tin_engine.fetch.nve import FIELDS, LAYERS
from tin_engine.gauge import Gauge, place
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_stations
from tin_engine.sources import STATION_SOURCES

HALF = 1000.0
data, here, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
src = STATION_SOURCES["nve-hrd"]
stations, s_crs = read_stations(data / "stations.geojson")
segments, crs, _ = read_segments(data / "rivers.geojson")
move = reprojector(s_crs, crs)
rows = {r["station"]: r for r in csv.DictReader(open(here / "results.csv", encoding="utf-8"))}
client = RangeClient()
lakes: dict[int, dict] = {}
n = 0
for st in stations:
    r = rows[st.station]
    if r["seeded_by"] != "river" and r["lake"] != "False":
        continue
    ((x, y),) = move([(st.x, st.y)])
    p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
    assert p is not None
    up = substring(LineString(p.reach.line), 0.0, p.reach.at)
    station_box = box(x - HALF, y - HALF, x + HALF, y + HALF)
    outside = up.difference(station_box).length
    x0, y0, x1, y1 = station_box.union(up).bounds
    params = {
        "where": "1=1",
        "geometry": f"{x0},{y0},{x1},{y1}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": "25833",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": ",".join(FIELDS[5]),
        "outSR": "25833",
        "f": "geojson",
    }
    doc = json.loads(client.get_text(query_url(f"{src.service_url}/{LAYERS[5]}/query", params)))
    if doc.get("exceededTransferLimit") or (doc.get("properties") or {}).get("exceededTransferLimit"):
        raise SystemExit(f"{st.station}: the answer was truncated")
    for f in doc["features"]:
        lakes.setdefault(int(f["properties"]["objectid"]), f)
    n += 1
    print(f"{st.station}\treach above P {up.length:.0f} m\toutside the ±1000 m box {outside:.0f} m")
coll = {
    "type": "FeatureCollection",
    "crs": {"type": "name", "properties": {"name": crs}},
    "features": list(lakes.values()),
}
out.write_text(json.dumps(coll, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"{n} queries, {len(lakes)} lakes")
