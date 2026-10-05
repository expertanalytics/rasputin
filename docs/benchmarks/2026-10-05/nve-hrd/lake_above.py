"""River-seeded rows with a lake on the mapped reach above P (@perf; the
list question 10 asks for). Usage: lake_above.py DATA HERE LAKES

Places each station seeded on the river, or placed on a river line (the
population of `@architect`'s probe, "NVE's lakes"), as the batch does (`gauge.place`, the
defaults) and tests whether the reach from its upstream end (up to 1000 m
above P) to P meets one of NVE's lake polygons in LAKES: lake_query.py's
output (lakes along every reach), not DATA/lakes.geojson, which holds only
the lakes within 100 m of each station point. Prints station, class,
how it was placed and seeded, and each lake's number and name, marked "(not in lakes.geojson)" when the
fetch's own file lacks it."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from shapely import STRtree
from shapely.geometry import LineString
from shapely.ops import substring

from tin_engine.crs import reprojector
from tin_engine.gauge import Gauge, place
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_lakes, read_stations

data, here, lakes_path = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
stations, s_crs = read_stations(data / "stations.geojson")
segments, crs, _ = read_segments(data / "rivers.geojson")
lakes, _ = read_lakes(lakes_path)
near = {(lk.number, lk.polygon.wkb) for lk in read_lakes(data / "lakes.geojson")[0]}
tree = STRtree([lk.polygon for lk in lakes])
move = reprojector(s_crs, crs)
rows = {r["station"]: r for r in csv.DictReader(open(here / "results.csv", encoding="utf-8"))}
n = 0
for st in stations:
    r = rows[st.station]
    if r["seeded_by"] != "river" and r["lake"] != "False":
        continue
    ((x, y),) = move([(st.x, st.y)])
    p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
    assert p is not None
    hits = tree.query(substring(LineString(p.reach.line), 0.0, p.reach.at), predicate="intersects")
    if len(hits):
        n += 1
        names = "; ".join(
            f"{lakes[int(k)].number} {lakes[int(k)].name}"
            + ("" if (lakes[int(k)].number, lakes[int(k)].polygon.wkb) in near else " (not in lakes.geojson)")
            for k in hits
        )
        line = "lake line" if r["lake"] == "True" else "river line"
        print(f"{st.station}\t{r['name']}\t{r['class']}\t{line}, seeded by {r['seeded_by']}\t{names}")
print(f"{n} stations")
