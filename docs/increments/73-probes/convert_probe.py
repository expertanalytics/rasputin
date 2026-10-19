"""Stream GML 3.2 into a GeoPackage through tin_engine's writer; time, peak RSS, longest edge."""

import resource
import sqlite3
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import shapely

from tin_engine.crs import srs_entry
from tin_engine.io.geopackage_write import create_features_table, create_geopackage, insert_features

G = "{http://www.opengis.net/gml/3.2}"
src, dst, chunk = sys.argv[1], sys.argv[2], int(sys.argv[3])
t0 = time.perf_counter()
srs = srs_entry("EPSG:25832")
conn = sqlite3.connect(dst)
create_geopackage(conn, srs)
cols = None
rows = []
n = 0
longest = 0.0
ctx = ET.iterparse(src, events=("start", "end"))
_, root = next(ctx)


def ring(r):
    return np.array(r.find(f"{G}LinearRing/{G}posList").text.split(), dtype=float).reshape(-1, 2)


for ev, el in ctx:
    if ev != "end" or el.tag != G + "featureMember":
        continue
    feat = el[0]
    poly = feat.find(f".//{G}Polygon")
    ext = ring(poly.find(f"{G}exterior"))
    longest = max(longest, float(np.hypot(*np.diff(ext, axis=0).T).max()))
    g = shapely.Polygon(ext, [ring(i) for i in poly.findall(f"{G}interior")])
    leaves = {c.tag.rpartition("}")[2]: c.text for c in feat if len(c) == 0}
    if cols is None:
        cols = [("gml_id", "esriFieldTypeString")] + [
            (k, "esriFieldTypeInteger" if (v or "").strip().isdigit() else "esriFieldTypeString")
            for k, v in leaves.items()
        ]
        create_features_table(
            conn, feat.tag.rpartition("}")[2], srs.srs_id, "MULTIPOLYGON", "fid", cols
        )
        name = feat.tag.rpartition("}")[2]
    n += 1
    rows.append((n, g, [feat.get(G + "id")] + [leaves.get(k) for k, _ in cols[1:]]))
    if len(rows) >= chunk:
        insert_features(conn, name, srs.srs_id, rows)
        rows = []
    root.clear()
insert_features(conn, name, srs.srs_id, rows)
conn.commit()
conn.close()
dt = time.perf_counter() - t0
print(
    f"features {n} seconds {dt:.2f} peak_rss_MiB "
    f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} "
    f"longest_edge_m {longest:.0f} cols {cols}"
)
