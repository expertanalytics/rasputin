"""BHO 2017 level-3 outlines of the São Francisco (76k, k = 1..9), and their
target-grid sizes. A measurement script, not production code; nothing imports it.

    python fetch_level3.py DATA_DIR OUT_CRS_WKT_FILE OUT_JSON

ANA's ``Divisao_de_bacias`` service stops at level 2, so a level-3 outline is
the union of the BHO 50k drainage areas whose ``COBACIA`` starts with ``76k``
(``SPR/BHO2017_50K_AREADRENAGEM``), fetched in EPSG:4674 as ``fetch_bho.py``
(2026-10-01/basin-piece) does for 5k. Writes ``DATA_DIR/76k_epsg4674.geojson``
and, per unit, its area and the target grid ``rasputin mesh`` would plan in
the given ``--out-crs`` (``target_grid_for`` at ANADEM's 30 m).
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import shapely
from shapely.geometry import mapping, shape

from tin_engine.crs import crs_label
from tin_engine.domain import DomainPolygon
from tin_engine.target_grid import target_grid_for

URL = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_50K_AREADRENAGEM/FeatureServer/0/query"


def _get(params: dict[str, str], tries: int = 5) -> dict:
    url = URL + "?" + urllib.parse.urlencode(params)
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                d = json.load(r)
            if "error" in d:
                raise RuntimeError(d["error"])
            return d
        except Exception as e:  # noqa: BLE001 - a flaky public server; retry, then raise
            if k == tries - 1:
                raise
            print(f"retry {k + 1}: {e}", file=sys.stderr)
            time.sleep(5 * (k + 1))
    raise AssertionError


def outline(prefix: str) -> tuple[shapely.Geometry, int]:
    ids = sorted(_get({"where": f"COBACIA LIKE '{prefix}%'", "returnIdsOnly": "true", "f": "json"})["objectIds"])
    parts = []
    for i in range(0, len(ids), 200):
        d = _get({"objectIds": ",".join(map(str, ids[i : i + 200])), "outFields": "OBJECTID",
                  "returnGeometry": "true", "outSR": "4674", "f": "geojson"})  # fmt: skip
        parts += [shape(f["geometry"]) for f in d["features"]]
    u = shapely.union_all(shapely.make_valid(shapely.GeometryCollection(parts)).geoms)
    if u.geom_type == "MultiPolygon":  # keep the largest part; slivers are join artefacts
        u = max(u.geoms, key=lambda g: g.area)
    return shapely.Polygon(u.exterior), len(ids)  # holes dropped: a sub-basin has none


def main() -> None:
    data, out = Path(sys.argv[1]), Path(sys.argv[3])
    wkt = crs_label(Path(sys.argv[2]).read_text().strip())
    rows = []
    for k in range(1, 10):
        p = f"76{k}"
        poly, n = outline(p)
        path = data / f"{p}_epsg4674.geojson"
        fc = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {"COBACIA": p},
              "geometry": mapping(poly)}]}  # fmt: skip
        path.write_text(json.dumps(fc))
        dom = DomainPolygon(polygon=poly, crs="EPSG:4674").to_crs(wkt)
        g = target_grid_for(dom, wkt, 30)
        rows.append({"prefix": p, "features": n, "vertices": len(poly.exterior.coords) - 1,
                     "area_km2": dom.polygon.area / 1e6, "grid_rows": g.rows, "grid_cols": g.cols,
                     "grid_nodes": g.rows * g.cols, "file": path.name})  # fmt: skip
        print(rows[-1], file=sys.stderr)
    out.write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
