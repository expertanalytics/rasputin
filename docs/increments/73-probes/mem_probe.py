"""tracemalloc peak of a clearing iterparse over N synthetic GML 3.2 features."""

import io
import tracemalloc
import xml.etree.ElementTree as ET

import numpy as np
import shapely

G = "{http://www.opengis.net/gml/3.2}"
HEAD = (
    b'<?xml version="1.0"?><gml:FeatureCollection'
    b' xmlns:gml="http://www.opengis.net/gml/3.2" xmlns:app="urn:x">'
)
FEAT = (
    b'<gml:featureMember><app:F gml:id="f%d"><app:a>60</app:a><app:g>'
    b'<gml:Polygon srsName="EPSG:25832">'
    b"<gml:exterior><gml:LinearRing><gml:posList>0 0 10 0 10 10 0 10 0 0</gml:posList>"
    b"</gml:LinearRing></gml:exterior>"
    b"</gml:Polygon></app:g></app:F></gml:featureMember>"
)


def doc(n):
    return io.BytesIO(HEAD + b"".join(FEAT % k for k in range(n)) + b"</gml:FeatureCollection>")


for n in (2_000, 20_000):
    s = doc(n)
    size = len(s.getvalue())
    tracemalloc.start()
    ctx = ET.iterparse(s, events=("start", "end"))
    _, root = next(ctx)
    k = 0
    for ev, el in ctx:
        if ev == "end" and el.tag == G + "featureMember":
            r = el[0].find(f".//{G}posList")
            shapely.Polygon(np.array(r.text.split(), float).reshape(-1, 2))
            k += 1
            root.clear()
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    print(n, "features", size, "bytes", "peak KiB", peak // 1024)
