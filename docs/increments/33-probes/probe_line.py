"""Design probe for increment 33: the Bergen Line in Bane NOR's Banenettverk (NLOD).

Run in a directory holding the national GML in EPSG:25833; see README.md.
"""

import collections
import glob
import pickle
import time
import xml.etree.ElementTree as ET

import numpy as np
import shapely
from shapely.geometry import LineString, MultiLineString, box

GML = "Samferdsel_0000_Norge_25833_Banenettverk_GML.gml"
DTM10 = "/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/*.tfw"
NS = {
    "app": "http://skjema.geonorge.no/SOSI/produktspesifikasjon/Banenettverk/1.0",
    "gml": "http://www.opengis.net/gml/3.2",
}


Link = tuple[str | None, str | None, LineString]


def read_links() -> dict[str, list[Link]]:
    links: dict[str, list[Link]] = collections.defaultdict(list)
    for _event, el in ET.iterparse(GML):
        if el.tag.endswith("}Banelenke"):
            name = el.findtext(".//app:banenavn", namespaces=NS) or ""
            status = el.findtext(".//app:banestatus", namespaces=NS)
            medium = el.findtext(".//app:medium", namespaces=NS)
            pos = el.find(".//gml:posList", NS)
            assert pos is not None and pos.text is not None
            xy = np.array(pos.text.split(), float).reshape(-1, 2)
            links[name].append((status, medium, LineString(xy)))
            el.clear()
    return links


def tiles() -> list[shapely.Polygon]:
    out = []
    for tfw in glob.glob(DTM10):
        with open(tfw) as f:
            a = [float(x) for x in f]
        out.append(box(a[4] - 5, a[5] - 50505, a[4] + 50505, a[5] + 5))
    return out


def main() -> None:
    t0 = time.perf_counter()
    links = read_links()
    print("parse s", round(time.perf_counter() - t0, 1), flush=True)
    counts = collections.Counter((s, m) for s, m, _ in links["Bergensbanen"])
    print("Bergensbanen (status, medium) link counts", counts, flush=True)
    for names in (["Bergensbanen"], ["Bergensbanen", "Randsfjordbanen", "Drammenbanen"]):
        sel = [(m, line) for n in names for s, m, line in links[n] if s == "I"]
        km = sum(line.length for _, line in sel) / 1e3
        tunnel = sum(line.length for m, line in sel if m == "U") / 1e3
        v = sum(len(line.coords) for _, line in sel)
        mean = round(km * 1e3 / (v - len(sel)), 1)
        print(names, "in operation: links", len(sel), "km", round(km, 1), "tunnel km",
              round(tunnel, 1), "vertices", v, "mean seg m", mean, flush=True)  # fmt: skip
        ml = MultiLineString([line for _, line in sel])
        t = time.perf_counter()
        simp = ml.simplify(1.0, preserve_topology=False)
        print("  simplify 1.0 m: vertices", len(shapely.get_coordinates(simp)), "s",
              round(time.perf_counter() - t, 2), flush=True)  # fmt: skip
        t = time.perf_counter()
        corr = shapely.buffer(ml.simplify(5.0), 5000, quad_segs=4)
        print("  5 km corridor km2", round(corr.area / 1e6), "s",
              round(time.perf_counter() - t, 1), flush=True)  # fmt: skip
        band3 = shapely.buffer(ml.simplify(5.0), 3000, quad_segs=4)
        print("  3 km band km2", round(band3.area / 1e6), flush=True)
        for r in (100, 500, 1000, 2000):
            area = shapely.buffer(ml.simplify(5.0), r, quad_segs=4).area / 1e6
            print("    band", r, "m km2", round(area), flush=True)
        with open(f"line_{len(names)}.pkl", "wb") as f:
            pickle.dump((ml, corr), f)
        hit = [t for t in tiles() if t.intersects(corr)]
        missing = corr.difference(shapely.union_all(hit)).area / 1e6
        print("  DTM10 tiles touching corridor", len(hit), "uncovered corridor km2",
              round(missing, 3), flush=True)  # fmt: skip


if __name__ == "__main__":
    main()
