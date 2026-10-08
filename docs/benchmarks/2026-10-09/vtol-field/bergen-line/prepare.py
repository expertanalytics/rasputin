"""Prepare increment 33's Bergen Line inputs from Bane NOR's Banenettverk (NLOD 1.0).

Run from the directory holding `Samferdsel_0000_Norge_25833_Banenettverk_GML.gml`
(README.md beside this script says where it comes from). Writes, in that directory:

- `bergensbanen.geojson`: Bergensbanen's links in operation (banestatus I), tunnels kept;
- `bergen_line3.geojson`: the same for Bergensbanen, Randsfjordbanen and Drammenbanen;
- `section_domain.geojson`: Geilo to Ål, as `docs/increments/33-probes/section.py`;
- `corridor_domain.geojson`: Hokksund to Bergen, the largest part of the three lines'
  5 km buffer simplified by 20 m, as the same probe.

All in EPSG:25833, each with a `crs` member. The selection and the domains repeat
`33-probes/probe_line.py` and `section.py`; only the file writing is new.
"""

import collections
import json
import xml.etree.ElementTree as ET

import numpy as np
import shapely
from shapely.geometry import LineString, MultiLineString, box, mapping

GML = "Samferdsel_0000_Norge_25833_Banenettverk_GML.gml"
NS = {
    "app": "http://skjema.geonorge.no/SOSI/produktspesifikasjon/Banenettverk/1.0",
    "gml": "http://www.opengis.net/gml/3.2",
}
CRS = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::25833"}}


def read_links() -> dict[str, list[tuple[str | None, LineString]]]:
    links: dict[str, list[tuple[str | None, LineString]]] = collections.defaultdict(list)
    for _event, el in ET.iterparse(GML):
        if el.tag.endswith("}Banelenke"):
            name = el.findtext(".//app:banenavn", namespaces=NS) or ""
            status = el.findtext(".//app:banestatus", namespaces=NS)
            pos = el.find(".//gml:posList", NS)
            assert pos is not None and pos.text is not None
            links[name].append((status, LineString(np.array(pos.text.split(), float).reshape(-1, 2))))
            el.clear()
    return links


def write(path: str, geoms: list[shapely.Geometry]) -> None:
    feats = [{"type": "Feature", "properties": {}, "geometry": mapping(g)} for g in geoms]
    with open(path, "w") as f:
        json.dump({"type": "FeatureCollection", "crs": CRS, "features": feats}, f)


def main() -> None:
    links = read_links()
    one = [line for s, line in links["Bergensbanen"] if s == "I"]
    three = [line for n in ("Bergensbanen", "Randsfjordbanen", "Drammenbanen") for s, line in links[n] if s == "I"]
    write("bergensbanen.geojson", one)
    write("bergen_line3.geojson", three)
    print("links", len(one), len(three))
    ml = MultiLineString(one)
    sec = shapely.line_merge(ml.intersection(box(127000, 6720000, 150000, 6745000)))
    grown = shapely.buffer(sec.simplify(5.0), 5000, quad_segs=4)
    dom = grown.intersection(box(127000, 6700000, 150000, 6770000)).simplify(20.0)
    print("section km2", round(dom.area / 1e6, 1), "vertices", len(shapely.get_coordinates(dom)))
    write("section_domain.geojson", [dom])
    corr = shapely.buffer(MultiLineString(three).simplify(5.0), 5000, quad_segs=4)
    big = max(corr.simplify(20.0).geoms, key=lambda g: g.area)
    print("corridor km2", round(big.area / 1e6), "vertices", len(shapely.get_coordinates(big)),
          "holes", len(big.interiors))  # fmt: skip
    write("corridor_domain.geojson", [big])


if __name__ == "__main__":
    main()
