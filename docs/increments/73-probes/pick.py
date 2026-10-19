"""Pick fixture candidates from a GML 3.2 slice: figures per feature."""

import sys
import xml.etree.ElementTree as ET

import numpy as np
import shapely

G = "{http://www.opengis.net/gml/3.2}"
ctx = ET.iterparse(sys.argv[1], events=("start", "end"))
_, root = next(ctx)
k = 0
shown_hole = 0
small_hole = None


def ring(r):
    return np.array(r.find(f"{G}LinearRing/{G}posList").text.split(), float).reshape(-1, 2)


for ev, el in ctx:
    if ev != "end" or el.tag != G + "featureMember":
        continue
    f = el[0]
    p = f.find(f".//{G}Polygon")
    ins = p.findall(f"{G}interior")
    e = ring(p.find(f"{G}exterior"))
    g = shapely.Polygon(e, [ring(i) for i in ins])
    leaf = {c.tag.rpartition("}")[2]: c.text for c in f if len(c) == 0}
    nbytes = len(ET.tostring(el))
    row = (
        k,
        f.get(G + "id"),
        len(ins),
        len(e),
        sum(len(ring(i)) for i in ins),
        round(g.area, 1),
        g.is_valid,
        leaf["arealtype"],
        leaf["treslag"],
        leaf["skogbonitet"],
        nbytes,
    )
    if k < 5:
        print("first", row)
    if ins and (small_hole is None or nbytes < small_hole[-1]) and len(ins) >= 2:
        small_hole = row
    k += 1
    root.clear()
print("smallest with >=2 holes", small_hole)
