# ruff: noqa: E402, E501
"""Increment 44 probe 3: GLO-30 tiles chosen, today's nominal spacing against the design's."""
import math
import sys

import tin_engine

assert tin_engine.__file__.startswith(sys.argv[1])
from tin_engine.fetch.plan import FetchRequest, source_box
from tin_engine.io.models import Bounds, RasterMeta


def chosen(box, dx, dy):
    meta = RasterMeta(x_min=0.0, y_max=0.0, delta_x=dx, delta_y=dy, cols=1, rows=1, epsg=4326,
                      nodata=None, nodata_source="absent", pixel_is_area=False, vertical_unit_assumed=True)
    b = source_box(FetchRequest(source="glo30", box=Bounds(**box)), meta)
    return b, [(la, lo) for la in range(math.floor(b.y_min), math.ceil(b.y_max))
               for lo in range(math.floor(b.x_min), math.ceil(b.x_max))]
cases = {
    "S87, east edge 0.0063 deg west of 166 E (10\" band)": (dict(x_min=165.2, y_min=-87.8, x_max=165.9937, y_max=-87.2), 10),
    "N82, east edge 0.0063 deg west of 21 E (5\" band)": (dict(x_min=20.2, y_min=82.2, x_max=20.9937, y_max=82.8), 5),
    "N60, every edge 0.003 deg inside 10..11 E, 60..61 N": (dict(x_min=10.003, y_min=60.003, x_max=10.997, y_max=60.997), 2),
}
for name, (box, s) in cases.items():
    out = []
    for label, dx, dy in [("today", 0.001, 0.001), ("own spacing", s / 3600, 1 / 3600), ("design", 0.0027778, 0.00027778)]:
        b, t = chosen(box, dx, dy)
        out.append(f"{label} {len(t)} tiles (x {b.x_max - box['x_max']:.6f}, y {b.y_max - box['y_max']:.6f})")
    print(name + ": " + "; ".join(out))
