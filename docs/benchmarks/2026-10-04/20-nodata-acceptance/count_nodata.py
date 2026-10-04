"""NoData nodes in the 1 m benchmark's DEM tile, in the whole tile and inside the quarter
domain (@perf). Reads the tile with tifffile; does not import tin_engine. Nodes at pixel
centres (PixelIsArea tiepoint + half a pixel). python count_nodata.py TILE QUARTER.geojson"""
import sys

import numpy as np
import shapely
import tifffile

tif, quarter = sys.argv[1], sys.argv[2]
with tifffile.TiffFile(tif) as t:
    p = t.pages[0]
    z = p.asarray()
    nd = float(str(p.tags["GDAL_NODATA"].value).strip("\x00"))
    tie = p.tags["ModelTiepointTag"].value
    sc = p.tags["ModelPixelScaleTag"].value
nod = z == nd
print(f"tile {z.shape[0]} x {z.shape[1]} nodes, NoData value {nd}: {int(nod.sum())} NoData nodes")
r, c = np.nonzero(nod)
x = tie[3] + (c + 0.5) * sc[0]
y = tie[4] - (r + 0.5) * sc[1]
dom = shapely.from_geojson(open(quarter).read())
if dom.geom_type == "GeometryCollection":
    dom = shapely.union_all(list(shapely.get_parts(dom)))
inside = shapely.contains_xy(dom, x, y) if len(x) else np.zeros(0, bool)
print(f"inside the quarter domain: {int(inside.sum())}")
if len(r):
    print(f"NoData rows {r.min()}-{r.max()}, cols {c.min()}-{c.max()}")
