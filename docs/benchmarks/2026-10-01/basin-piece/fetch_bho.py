"""Fetch BHO 2017 (5k) drainage areas from ANA's ArcGIS REST service.

A stopgap, not production code: rasputin should fetch its own data one day
(ROADMAP, the basin's own inputs); until then this script is how the basin
piece's outline was obtained. Nothing imports it.

    python fetch_bho.py table  DATA_DIR            # attributes of every São Francisco
                                                   # elementary catchment (prefix 76)
    python fetch_bho.py outline DATA_DIR PREFIX    # dissolve one ottobasin into a polygon

The service is ``BHO2017_5K_AREADRENAGEM`` (ANA, "Base Hidrográfica Ottocodificada
2017 5K - área de drenagem"), served in EPSG:3857; geometry is requested in
EPSG:4674 (SIRGAS 2000, BHO's own CRS) so no web-Mercator round trip is made
on our side. The server refuses count and statistics queries; ``COBACIA LIKE``
with ``returnIdsOnly`` works, then features are fetched by object id.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import shapely
from shapely.geometry import shape

URL = "https://www.snirh.gov.br/arcgis/rest/services/SPR/BHO2017_5K_AREADRENAGEM/FeatureServer/0/query"
FIELDS = "OBJECTID,COBACIA,NUAREACONT,NUNIVOTTO4,NUNIVOTTO5,NUNIVOTTO6,COCURSODAG,DSVERSAO"


def _post(params: dict[str, str], tries: int = 5) -> dict:
    # GET, not POST: the same query as a POST times out at the gateway (504).
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


def ids(prefix: str) -> list[int]:
    d = _post({"where": f"COBACIA LIKE '{prefix}%'", "returnIdsOnly": "true", "f": "json"})
    return sorted(d.get("objectIds") or [])  # a prefix with no catchments gives null


def table(data: Path) -> None:
    # One query for all of "76" times out at the gateway; ten level-3 queries do not.
    all_ids = sorted(i for k in range(10) for i in ids(f"76{k}"))
    out = data / "bho2017_5k_sao_francisco_attributes.csv"
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(FIELDS.split(","))
        for i in range(0, len(all_ids), 500):
            chunk = ",".join(map(str, all_ids[i : i + 500]))
            d = _post({"objectIds": chunk, "outFields": FIELDS, "returnGeometry": "false", "f": "json"})
            for feat in d["features"]:
                a = feat["attributes"]
                w.writerow([a[k] for k in FIELDS.split(",")])
            print(f"{min(i + 500, len(all_ids))}/{len(all_ids)}", file=sys.stderr)
    print(out)


def outline(data: Path, prefix: str) -> None:
    got = []
    pids = ids(prefix)
    for i in range(0, len(pids), 200):
        chunk = ",".join(map(str, pids[i : i + 200]))
        d = _post({"objectIds": chunk, "outFields": FIELDS, "returnGeometry": "true",
                   "outSR": "4674", "f": "geojson"})  # fmt: skip
        got += d["features"]
        print(f"{len(got)}/{len(pids)}", file=sys.stderr)
    raw = data / f"bho2017_5k_{prefix}_elementary.geojson"
    raw.write_text(json.dumps({"type": "FeatureCollection", "features": got}))
    geoms = [shape(f["geometry"]) for f in got]
    merged = shapely.union_all(geoms)
    parts = list(getattr(merged, "geoms", [merged]))
    parts.sort(key=lambda g: g.area, reverse=True)
    main = parts[0]
    print(f"{len(got)} catchments; union has {len(parts)} part(s); areas (deg²) "
          f"{[round(p.area, 8) for p in parts[:5]]}; main interiors {len(main.interiors)}",
          file=sys.stderr)  # fmt: skip
    # Slivers between neighbouring elementary polygons leave tiny holes; drop
    # holes under 1e-6 deg² (about 0.012 km²) and report how many were dropped.
    small = [r for r in main.interiors if shapely.Polygon(r).area < 1e-6]
    main = shapely.Polygon(main.exterior, [r for r in main.interiors if shapely.Polygon(r).area >= 1e-6])
    print(f"dropped {len(small)} sliver holes; kept {len(main.interiors)}", file=sys.stderr)
    area = sum(float(f["properties"]["NUAREACONT"]) for f in got)
    out = data / f"bho2017_5k_{prefix}_outline_epsg4674.geojson"
    out.write_text(json.dumps({
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::4674"}},
        "features": [{"type": "Feature", "properties": {
            "ottobasin": prefix, "elementary_catchments": len(got),
            "sum_NUAREACONT_km2": round(area, 3), "source": URL,
        }, "geometry": shapely.geometry.mapping(main)}],
    }))  # fmt: skip
    print(out)


if __name__ == "__main__":
    cmd, data = sys.argv[1], Path(sys.argv[2])
    data.mkdir(parents=True, exist_ok=True)
    {"table": lambda: table(data), "outline": lambda: outline(data, sys.argv[3])}[cmd]()
