"""River-seeded rows with a lake on the mapped reach above P (@perf; the
list question 10 asks for). Usage: lake_above.py DATA HERE

Places each river-seeded station as the batch does (`gauge.place`, the
defaults) and tests whether the reach from its upstream end (up to 1000 m
above P) to P meets one of NVE's lake polygons in lakes.geojson. Prints
station, class, and the lake's number and name."""

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

data, here = Path(sys.argv[1]), Path(sys.argv[2])
stations, s_crs = read_stations(data / "stations.geojson")
segments, crs, _ = read_segments(data / "rivers.geojson")
lakes, _ = read_lakes(data / "lakes.geojson")
tree = STRtree([lk.polygon for lk in lakes])
move = reprojector(s_crs, crs)
rows = {r["station"]: r for r in csv.DictReader(open(here / "results.csv", encoding="utf-8"))}
n = 0
for st in stations:
    r = rows[st.station]
    if r["seeded_by"] != "river":
        continue
    ((x, y),) = move([(st.x, st.y)])
    p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
    assert p is not None
    hits = tree.query(substring(LineString(p.reach.line), 0.0, p.reach.at), predicate="intersects")
    if len(hits):
        n += 1
        names = "; ".join(f"{lakes[int(k)].number} {lakes[int(k)].name}" for k in hits)
        print(f"{st.station}\t{r['name']}\t{r['class']}\t{names}")
print(f"{n} stations")
