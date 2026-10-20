"""Stream a GML 3.2 file with iterparse; build shapely polygons; report time, peak RSS."""

import resource
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import shapely

G = "{http://www.opengis.net/gml/3.2}"
path = sys.argv[1]
chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 0
t0 = time.perf_counter()
n = rings = holes = verts = 0
batch = []
biggest = 0
ctx = ET.iterparse(path, events=("start", "end"))
_, root = next(ctx)
for ev, el in ctx:
    if ev != "end" or el.tag != G + "featureMember":
        continue
    feat = el[0]
    poly = feat.find(f".//{G}Polygon")

    def ring(r):
        a = np.array(r.find(f"{G}LinearRing/{G}posList").text.split(), dtype=float).reshape(-1, 2)
        return a

    ext = ring(poly.find(f"{G}exterior"))
    ins = [ring(i) for i in poly.findall(f"{G}interior")]
    g = shapely.Polygon(ext, ins)
    k = len(ext) + sum(len(i) for i in ins)
    verts += k
    biggest = max(biggest, k)
    holes += len(ins)
    attrs = {c.tag.rpartition("}")[2]: c.text for c in feat if len(c) == 0}
    batch.append((n, g, attrs))
    n += 1
    if chunk and len(batch) >= chunk:
        batch.clear()
    root.clear()
dt = time.perf_counter() - t0
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
print(
    f"features {n} holes {holes} vertices {verts} biggest {biggest} "
    f"seconds {dt:.2f} peak_rss_MiB {rss:.0f} kept {len(batch)}"
)
