"""Read a converted GeoPackage through tin_engine's existing vector readers."""

import sys
import time
from contextlib import closing
from pathlib import Path

import numpy as np
from pyproj import Transformer

from tin_engine.feature_input import read_source
from tin_engine.io.geopackage import distinct_values, layer_info, layer_summary
from tin_engine.io.repository import open_geopackage

path = Path(sys.argv[1])
with closing(open_geopackage(path)) as conn:
    info = layer_info(conn, None)
    s = layer_summary(conn, info)
    print(
        info.table,
        info.crs,
        info.pk,
        info.rtree is not None,
        s.features,
        s.geometry,
        s.extent,
        s.attributes,
    )
    for a in ("arealtype", "treslag", "skogbonitet"):
        print(a, distinct_values(conn, info, a, 100))
x0, y0 = 600000.0, 6800000.0
t = time.perf_counter()
rows = read_source(path, None, "arealtype", lambda own: (x0, y0, x0 + 10000, y0 + 10000))
print(
    "read_source 10 km box:",
    len(rows.rows),
    "rows",
    rows.crs,
    f"{time.perf_counter() - t:.2f}s",
    rows.rows[0][2],
)
# straight 4.4 km edge in UTM 32 near the zone's east edge, moved to UTM 33: chord deviation
tr = Transformer.from_crs(25832, 25833, always_xy=True)
for xa in (500000.0, 700000.0):
    a, b = np.array([xa, 6800000.0]), np.array([xa + 3111.0, 6803111.0])
    pts = np.array([a + (b - a) * k for k in np.linspace(0, 1, 101)])
    X, Y = tr.transform(pts[:, 0], pts[:, 1])
    q = np.c_[X, Y]
    d = q[-1] - q[0]
    n = np.array([-d[1], d[0]]) / np.hypot(*d)
    print(
        "UTM32 x",
        xa,
        "edge 4.4 km, max deviation from chord in UTM33 (mm):",
        1000 * np.abs((q - q[0]) @ n).max(),
    )
