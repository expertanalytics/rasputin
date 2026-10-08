"""Design probe for increment 33: the Geilo-Ål section, its 5 km domain, and the corridor."""

import json
import pickle

import shapely
from shapely.geometry import box, mapping


def main() -> None:
    with open("line_1.pkl", "rb") as f:
        ml, _corr = pickle.load(f)
    sec = shapely.line_merge(ml.intersection(box(127000, 6720000, 150000, 6745000)))
    print("section km (all links)", round(sec.length / 1e3, 1))
    grown = shapely.buffer(sec.simplify(5.0), 5000, quad_segs=4)
    dom = grown.intersection(box(127000, 6700000, 150000, 6770000)).simplify(20.0)
    print("domain km2", round(dom.area / 1e6, 1), "vertices",
          len(shapely.get_coordinates(dom)), dom.geom_type)  # fmt: skip
    write_domain("section_domain.geojson", dom)
    with open("section.pkl", "wb") as f:
        pickle.dump((sec, dom), f)
    # The corridor: the largest part of the three lines' 5 km buffer, simplified by 20 m.
    with open("line_3.pkl", "rb") as f:
        _ml3, corr = pickle.load(f)
    big = max(corr.simplify(20.0).geoms, key=lambda g: g.area)
    print("corridor km2", round(big.area / 1e6), "vertices", len(shapely.get_coordinates(big)),
          "holes", len(big.interiors))  # fmt: skip
    write_domain("corridor_domain.geojson", big)


def write_domain(path: str, polygon: shapely.Polygon) -> None:
    doc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::25833"}},
        "features": [{"type": "Feature", "properties": {}, "geometry": mapping(polygon)}],
    }
    with open(path, "w") as f:
        json.dump(doc, f)


if __name__ == "__main__":
    main()
