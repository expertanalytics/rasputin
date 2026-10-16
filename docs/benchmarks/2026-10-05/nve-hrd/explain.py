"""Step 4's measurements for each finding (@perf). No figure here is timed.

Usage: explain.py DATA DEM HERE HALF_M STATION [STATION ...]

For each station:
* ours against NVE's polygon (our reduced outline from HERE/catchments):
  the area only in ours and only in NVE's, km², and the largest connected
  piece of each with its distance from the station;
* where the DEM's water runs near the gauge: `_core.accumulate` on the RAW
  DEM (no burn) in a window of P +-HALF_M; the node of largest count
  within 60, 300 and 1000 m of the mapped position P, its count in km² (cells x 100 m2,
  a lower bound when the window cuts the catchment), and its distance from
  P and from the mapped line; each a lower bound when the window cuts its catchment;
* the lakes within 1 km of P (from NVE's lakes file).
Prints one block per station.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Point, shape

from tin_engine._core import accumulate
from tin_engine.catchment import _plan
from tin_engine.crs import reprojector
from tin_engine.dem_input import repository_for
from tin_engine.gauge import Gauge, place
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_lakes, read_references, read_stations
from tin_engine.mosaic import Bounds, assemble
from tin_engine.raster import to_core


def pieces(g):  # type: ignore[no-untyped-def]
    return sorted(getattr(g, "geoms", [g]), key=lambda p: -p.area) if not g.is_empty else []


def main(data: Path, dem: Path, here: Path, ids: list[str], half: float = 3000.0) -> None:
    stations, s_crs = read_stations(data / "stations.geojson")
    segments, crs, _ = read_segments(data / "rivers.geojson")
    refs, _ = read_references(data / "reference.geojson")
    lakes, _ = read_lakes(data / "lakes.geojson")
    repository, _ = repository_for((dem,))
    footprints = repository.footprints()
    move = reprojector(s_crs, crs)
    by = {s.station: s for s in stations}
    for sid in ids:
        st = by[sid]
        ((x, y),) = move([(st.x, st.y)])
        station = Point(x, y)
        p = place(Gauge(x=x, y=y, watercourse=st.watercourse, river=st.river), segments)
        print(f"== {sid} {st.name} (watercourse {st.watercourse}, river {st.river})")
        ref = refs[sid]
        f = here / "catchments" / f"{sid}.geojson"
        if f.exists():
            feat = json.loads(f.read_text())
            feat = feat["features"][0] if "features" in feat else feat
            ours = shape(feat["geometry"])
            for label, d in (("only ours", ours.difference(ref)), ("only NVE's", ref.difference(ours))):
                ps = pieces(d)
                big = ps[0] if ps else None
                line = f"  {label}: {d.area / 1e6:.2f} km2 in {len(ps)} pieces"
                if big is not None:
                    c = big.representative_point()
                    line += (f"; largest {big.area / 1e6:.2f} km2, {big.distance(station) / 1e3:.1f} km "
                             f"from the station, a point in it ({c.x:.0f}, {c.y:.0f})")
                print(line)
            print(f"  station inside NVE's polygon: {ref.contains(station)}, "
                  f"{ref.exterior.distance(station) if hasattr(ref, 'exterior') else ref.boundary.distance(station):.0f} m from its edge;"
                  f" inside ours: {ours.contains(station)}")
        if p is None:
            print("  no placement")
            continue
        px, py = p.position
        print(f"  placed on {p.placed_on}, line {p.objectid} ({p.elvid}), {p.distance_m:.0f} m; "
              f"P = ({px:.0f}, {py:.0f}); lake line {p.lake}; confluence near {p.confluence_near}")
        plan = _plan(footprints, Bounds(x_min=px - half, y_min=py - half, x_max=px + half, y_max=py + half))
        tile = assemble(plan, repository.load).tile
        m = tile.meta
        acc = accumulate(to_core(tile))
        count = np.asarray(acc.count)
        rows, cols = np.indices(count.shape)
        xs, ys = m.x_min + cols * m.delta_x, m.y_max - rows * m.delta_y
        dist = np.hypot(xs - px, ys - py)
        mapped = LineString(p.reach.line)
        for radius in (60.0, 300.0, 1000.0):
            near = dist <= radius
            k = np.argmax(np.where(near, count, 0))
            r, c = np.unravel_index(k, count.shape)
            q = Point(xs[r, c], ys[r, c])
            print(f"  raw DEM, largest count within {radius:.0f} m of P: {count[r, c] * m.delta_x * m.delta_y / 1e6:.2f} km2 "
                  f"at ({q.x:.0f}, {q.y:.0f}), {q.distance(Point(px, py)):.0f} m from P, "
                  f"{q.distance(mapped):.0f} m from the mapped line; inside NVE's polygon {ref.contains(q)}")
        near_lakes = [lk for lk in lakes if lk.polygon.distance(Point(px, py)) <= 1000.0]
        for lk in near_lakes:
            print(f"  lake {lk.number} {lk.name}: {lk.polygon.area / 1e6:.2f} km2, "
                  f"{lk.polygon.distance(Point(px, py)):.0f} m from P, {lk.polygon.distance(station):.0f} m from the station")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[5:], float(sys.argv[4]))
