"""Lake and glacier shares of one NVE gauge catchment from NVE's lake database,
NVE's glacier units (Bre2/MapServer/9) and the cached CORINE 2018 extract.

Run: .venv/bin/python -I docs/increments/107-probes/prefer_probe.py 73.27.0 \
    ../rasputin_data/cache/clc2018/features.gpkg   (and 76.5.0). Reads NVE's map services.
"""

import json
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import shape

from tin_engine.feature_input import read_class_polygons
from tin_engine.legends import load

NVE = "https://kart.nve.no/enterprise/rest/services"
WATER = {512, 511, 521, 522}


def query(url: str, **params: str) -> dict:
    params |= {"f": "geojson", "outSR": "25833"}
    with urllib.request.urlopen(f"{url}/query?{urllib.parse.urlencode(params)}", timeout=120) as r:
        return json.load(r)


def union(features: list) -> shapely.Geometry:
    return shapely.make_valid(shapely.union_all([shape(f["geometry"]) for f in features]))


def main(station: str, clc: Path) -> None:
    fc = query(f"{NVE}/HydrologiskeData3/MapServer/38", where=f"stasjonnr='{station}'")
    c = union(fc["features"])
    a = c.area
    print(f"catchment {a / 1e6:.2f} km2")
    box = {
        "geometry": ",".join(map(str, c.bounds)),
        "geometryType": "esriGeometryEnvelope",
        "inSR": "25833",
        "spatialRel": "esriSpatialRelIntersects",
        "where": "1=1",
    }
    lakes = union(query(f"{NVE}/Innsjodatabase2/MapServer/5", **box)["features"]).intersection(c)
    glaciers = union(query(f"{NVE}/Bre2/MapServer/9", **box)["features"]).intersection(c)
    rings = [np.asarray(p.exterior.coords) for p in shapely.get_parts(c)]
    found, _ = read_class_polygons(clc, None, load("corine"), "EPSG:25833", np.vstack(rings))
    water = shapely.union_all([g for k, g in found if k in WATER]).intersection(c)
    ice = shapely.union_all([g for k, g in found if k == 335]).intersection(c)
    print(f"lake: nve {lakes.area / a:.2%}, clc {water.area / a:.2%}")
    print(f"(a) nve lake not clc water {lakes.difference(water).area / a:.2%}")
    print(f"(b) clc water not nve lake {water.difference(lakes).area / a:.2%}")
    print(f"glacier: nve {glaciers.area / a:.2%}, clc {ice.area / a:.2%}")
    print(f"(c) nve lake inside nve glacier {glaciers.intersection(lakes).area / a:.4%}")
    under: Counter[int] = Counter()
    for k, g in found:
        under[k] += g.intersection(c).intersection(lakes).area
    shares = {k: round(v / a * 100, 2) for k, v in under.most_common(6)}
    print("clc classes under nve lakes, %:", shares)


if __name__ == "__main__":
    main(sys.argv[1], Path(sys.argv[2]))
