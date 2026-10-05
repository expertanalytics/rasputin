"""Step 6's checks and increment 22's outline guarantees, per station (@perf).

Usage: checks.py DATA DEM BATCH_DIR OUT_CSV

Reads the batch's results.csv and catchment files in BATCH_DIR, and for every
station that is not refused:

* oracle: the placed node's count (`a0` / cell area) equals the catchment's
  node count before reduction (`nodes`); river rows only (a lake row has no
  placed node);
* monotone, drains: copied from the row (the sensitivity's own flags);
* lowered_off_chain: the burn re-run through `burn.burn_reach` on stage A's
  first window (the reach's bounds plus the corridor and 2 km, as
  `catchment._gauged` builds it), every node whose height changed, less the
  chain's nodes; must be 0;
* the reduced outline from the catchment file: `simple` (shapely is_valid and
  the exterior is_simple), `area_diff_m2` (reduced minus fine, from the
  file's properties), `seed_inside` (strictly inside the reduced outline:
  the placed node from that burn on a river row, the lake seed's point on a
  lake row);
* lake_above_p: on a river row, whether the mapped reach upstream of P
  (up to `reach_up`, 1000 m) meets a lake polygon (the list "NVE's lakes"
  counts 24).
No figure here is timed.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from shapely import STRtree
from shapely.geometry import LineString, Point, shape
from shapely.ops import substring

from tin_engine.burn import burn_reach
from tin_engine.catchment import WINDOW_MARGIN_M, _plan
from tin_engine.crs import reprojector
from tin_engine.dem_input import repository_for
from tin_engine.gauge import Gauge, lake_seed, place
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_lakes, read_stations
from tin_engine.mosaic import Bounds, assemble


def main(data: Path, dem: Path, batch: Path, out: Path) -> None:
    stations, s_crs = read_stations(data / "stations.geojson")
    segments, crs, _ = read_segments(data / "rivers.geojson")
    lakes, _ = read_lakes(data / "lakes.geojson")
    tree = STRtree([lk.polygon for lk in lakes])
    repository, _ = repository_for((dem,))
    footprints = repository.footprints()
    move = reprojector(s_crs, crs)
    rows = {r["station"]: r for r in csv.DictReader(open(batch / "results.csv"))}
    fields = [
        "station", "class", "seeded_by", "oracle", "placed_count", "nodes", "monotone",
        "drains", "lowered_off_chain", "simple", "fine_area_m2", "area_diff_m2", "seed_inside", "lake_above_p",
    ]  # fmt: skip
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for st in stations:
            r = rows.get(st.station)
            if r is None or r["class"] == "refused":
                continue
            ((x, y),) = move([(st.x, st.y)])
            gauge = Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river)
            placement = place(gauge, segments)
            rec: dict[str, object] = {"station": st.station, "class": r["class"]}
            rec["seeded_by"] = r["seeded_by"]
            seed_xy = None
            if r["seeded_by"] == "lake":
                seed = lake_seed(gauge, placement, lakes)
                assert seed is not None, st.station
                seed_xy = seed.point
            else:
                assert placement is not None, st.station
                reach = placement.reach
                line = np.asarray(reach.line)
                pad = reach.corridor + WINDOW_MARGIN_M
                (x0, y0), (x1, y1) = line.min(axis=0) - pad, line.max(axis=0) + pad
                plan = _plan(footprints, Bounds(x_min=x0, y_min=y0, x_max=x1, y_max=y1))
                tile = assemble(plan, repository.load).tile
                burnt, path = burn_reach(tile, reach)
                changed = {tuple(v) for v in np.argwhere(burnt.array != tile.array).tolist()}
                chain = {tuple(v) for v in path.chain.tolist()}
                rec["lowered_off_chain"] = len(changed - chain)
                m = tile.meta
                pr, pc = path.chain[path.placed]
                seed_xy = (m.x_min + pc * m.delta_x, m.y_max - pr * m.delta_y)
                cell_km2 = m.delta_x * m.delta_y / 1e6
                placed = round(float(r["a0"]) / cell_km2)
                rec |= {"placed_count": placed, "nodes": int(r["nodes"])}
                rec["oracle"] = placed == int(r["nodes"])
                rec |= {"monotone": r["monotone"], "drains": r["drains"]}
                up = substring(LineString(reach.line), 0.0, reach.at)  # the reach starts reach_up above P
                hits = tree.query(up, predicate="intersects")
                rec["lake_above_p"] = bool(len(hits))
            f = json.loads((batch / f"{st.station}.geojson").read_text())
            feat = f["features"][0] if "features" in f else f
            poly, props = shape(feat["geometry"]), feat["properties"]
            rec["simple"] = bool(poly.is_valid and poly.exterior.is_simple)
            rec["fine_area_m2"] = props["fine_area_m2"]
            rec["area_diff_m2"] = props["reduced_area_m2"] - props["fine_area_m2"]
            rec["seed_inside"] = bool(poly.contains(Point(seed_xy)))
            w.writerow(rec)
            fh.flush()
            print(st.station, rec, file=sys.stderr)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
