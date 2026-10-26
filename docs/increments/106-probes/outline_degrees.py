"""Fix 5 probe: an outline in degrees (`catchment --dem glo30` writes EPSG:4326)
simplified in a local equal-area CRS, the Lambert azimuthal equal-area centred
on the outline's centroid, at band = h / 8 with A its area there, then moved
back vertex by vertex. Prints vertices, band, the geodesic area before and
after, validity and seconds; writes OUT/STEM-N.geojson (EPSG:4326).
Run:  .venv/bin/python -I docs/increments/106-probes/outline_degrees.py IN.geojson OUT N
"""

import json
import sys
import time
from pathlib import Path

import shapely
from pyproj import CRS, Transformer
from shapely.geometry import mapping, shape

from tin_engine import sub_catchments
from tin_engine.geodesy import area_m2

src, out, n = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
given = shape(json.loads(src.read_text())["features"][0]["geometry"])
c = given.centroid
laea = CRS.from_proj4(f"+proj=laea +lat_0={c.y} +lon_0={c.x} +ellps=GRS80 +units=m")
there = Transformer.from_crs("EPSG:4326", laea, always_xy=True)
back = Transformer.from_crs(laea, "EPSG:4326", always_xy=True)
t = time.perf_counter()
flat = shapely.transform(given, lambda xy: shapely.get_coordinates(
    shapely.points(*there.transform(xy[:, 0], xy[:, 1]))))
band = sub_catchments.band_m(flat.area, n)
(kept,) = sub_catchments.simplify([flat], band)
kept = shapely.transform(kept, lambda xy: shapely.get_coordinates(
    shapely.points(*back.transform(xy[:, 0], xy[:, 1]))))
s = time.perf_counter() - t
a0, a1 = area_m2(given, "EPSG:4326"), area_m2(kept, "EPSG:4326")
print(f"{src.stem}: N {n} band {band:.1f} m vertices {len(given.exterior.coords) - 1} -> "
      f"{len(kept.exterior.coords) - 1}, geodesic area {a0 / 1e6:.4f} -> {a1 / 1e6:.4f} km2 "
      f"(relative {(a1 - a0) / a0:.2g}), valid {kept.is_valid}, {s:.2f} s")
doc = {"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
       "features": [{"type": "Feature", "properties": {}, "geometry": mapping(kept)}]}
(out / f"{src.stem}-{n}.geojson").write_text(json.dumps(doc))
