"""Why `mesh --simplify-domain` keeps 198 of Ljungan's vertices where the
design's probe kept 195: the same ring simplified clockwise, as in the file,
and counter-clockwise, as `read_domain` orients it.

Run:  .venv/bin/python -I docs/benchmarks/2026-10-26/106e/orient_probe.py \
      ~/rasputin-output/opal-se2-2026-10-21/catchments-original/ljungan.geojson
"""

import json
import sys
from pathlib import Path

from shapely.geometry import Polygon, shape
from shapely.geometry.polygon import orient

from tin_engine import sub_catchments
from tin_engine.io.domain_file import read_domain


def kept(polygon: Polygon, band: float) -> int:
    """The outer ring's vertices after the simplification, closing repeat out."""
    return len(sub_catchments.simplify([polygon], band)[0].exterior.coords) - 1


path = Path(sys.argv[1])
raw = shape(json.loads(path.read_text())["features"][0]["geometry"])
band = sub_catchments.band_m(raw.area, 500)
print("as in the file, ccw", raw.exterior.is_ccw, "->", kept(raw, band))
for sign in (1.0, -1.0):
    turned = orient(raw, sign)
    print("orient", sign, "ccw", turned.exterior.is_ccw, "->", kept(turned, band))
read = read_domain(path).polygon
same_start = read.exterior.coords[0] == raw.exterior.coords[0]
print("read_domain, ccw", read.exterior.is_ccw, "same start", same_start, "->", kept(read, band))
