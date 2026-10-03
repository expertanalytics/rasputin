"""Throwaway probe for docs/increments/26-land-cover.md: how many MapBiomas
blocks and bytes a domain's box needs, from the COG header alone (one range
request), and which level-3 unit holds the big reservoirs.

    python docs/increments/26-probes/fetch_size.py $WORK OUTLINE.geojson [OUTLINE.geojson ...]
"""

from __future__ import annotations

import io
import json
import logging
import sys
import urllib.request
from pathlib import Path

import numpy as np
import shapely
import tifffile

logging.disable(logging.CRITICAL)
URL = (
    "https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection11/"
    "lulc/coverage/brazil_coverage/brazil_coverage-col11_2025.tif"
)
RESERVOIRS = {
    "Sobradinho": (-41.5, -9.6),
    "Tres Marias": (-45.25, -18.2),
    "Itaparica": (-38.5, -9.1),
}

work = Path(sys.argv[1])
head = work / "c11_2025_head.bin"
if not head.exists():
    req = urllib.request.Request(URL, headers={"Range": f"bytes=0-{(16 << 20) - 1}"})
    head.write_bytes(urllib.request.urlopen(req).read())
page = tifffile.TiffFile(io.BytesIO(head.read_bytes())).pages[0]
sx, sy = page.tags["ModelPixelScaleTag"].value[:2]
x0, y0 = page.tags["ModelTiepointTag"].value[3:5]
tw, th = page.tilewidth, page.tilelength
ncols = -(-page.shape[1] // tw)
counts = np.asarray(page.databytecounts)
print(
    json.dumps(
        {
            "tiles": int(counts.size),
            "sparse": int((counts == 0).sum()),
            "header_bytes_needed": int(min(o for o in page.dataoffsets if o)),
        }
    )
)
for path in sys.argv[2:]:
    gj = json.loads(Path(path).read_text())
    shape = shapely.union_all([shapely.geometry.shape(f["geometry"]) for f in gj["features"]])
    lon0, lat0, lon1, lat1 = shape.bounds
    c0, c1 = int((lon0 - x0) / sx) // tw, int((lon1 - x0) / sx) // tw
    r0, r1 = int((y0 - lat1) / sy) // th, int((y0 - lat0) / sy) // th
    idx = np.array([r * ncols + c for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)])
    b = counts[idx]
    # tiles meeting the outline itself (the outline-aware alternative 23a-2 did not take)
    boxes = shapely.box(
        x0 + (idx % ncols) * tw * sx,
        y0 - (idx // ncols + 1) * th * sy,
        x0 + (idx % ncols + 1) * tw * sx,
        y0 - (idx // ncols) * th * sy,
    )
    meet = shapely.intersects(shape, boxes)
    print(
        json.dumps(
            {
                "outline": Path(path).name,
                "box_tiles": int(idx.size),
                "box_bytes": int(b.sum()),
                "box_sparse": int((b == 0).sum()),
                "outline_tiles": int(meet.sum()),
                "outline_bytes": int(b[meet].sum()),
                "box_cells": int(idx.size * tw * th),
                "holds": [
                    k for k, (x, y) in RESERVOIRS.items() if shape.contains(shapely.Point(x, y))
                ],
            }
        )
    )
