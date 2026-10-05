"""`rasputin fetch-stations nve-hrd`: NVE's HRD stations, catchments and rivers (29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "The station set" and
"Data use". The packaged list says which stations; NVE's map services give
their points (`HydrologiskeData3` layer 0), their catchment polygons (layer
38), the river lines round each (`Elvenett1` layer 2) and, for PR 4's lake
gauges, the lakes near each (`Innsjodatabase2` layer 5). Every request names
an explicit field allow-list, never `*`, and asks for GeoJSON in EPSG:25833.

The network call is injected (`get_text`), and the requests go out one at a
time. Nothing here writes a file: :func:`fetch_station_set` returns each
file's bytes, the manifest last, and `cli.py` writes them.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.resources
import io
import json
from collections.abc import Callable, Mapping, Sequence
from datetime import UTC, datetime
from typing import Any

from tin_engine.sources import StationSource, notice

from .http import FetchError, query_url

#: "Data use": the fields asked of each layer, and the only ones copied.
FIELDS: Mapping[int, tuple[str, ...]] = {
    0: (
        "stasjonnr",
        "stasjonnavn",
        "totalt_feltareal_km2",
        "stasjonstatus",
        "vassdragsnr",
        "elvenavnhierarki",
    ),
    38: ("stasjonnr", "nedborfeltaareal_km2", "oppdateringsdato", "objectid"),
    2: ("objectid", "objekttype", "strekninglnr", "elvid", "vassdragsnr", "elvenavn", "vatnlnr"),
    5: ("objectid", "vatnlnr", "navn", "areal_km2"),
}
LAYERS = {
    0: "HydrologiskeData3/MapServer/0",
    38: "HydrologiskeData3/MapServer/38",
    2: "Elvenett1/MapServer/2",
    5: "Innsjodatabase2/MapServer/5",
}
#: Station numbers per `where stasjonnr in (...)` query.
CHUNK = 40
#: `map_radius + reach_up + 500 m`, each way round the station point.
ENVELOPE_HALF = 2000.0
#: Each way round the station point: any lake within `gauge.LAKE_GAP_M` meets it.
LAKE_ENVELOPE_HALF = 100.0
FILES = ("stations.geojson", "reference.geojson", "rivers.geojson", "lakes.geojson", "NOTICE.txt")


def read_list(source: StationSource) -> list[dict[str, str]]:
    """The packaged list's rows, in file order; `#` lines are comments."""
    text = (importlib.resources.files("tin_engine") / "data" / source.list_file).read_text(
        encoding="utf-8"
    )
    body = "\n".join(line for line in text.splitlines() if not line.startswith("#"))
    return list(csv.DictReader(io.StringIO(body)))


def _url(source: StationSource, layer: int, params: Mapping[str, str]) -> str:
    epsg = source.crs.split(":")[-1]
    fields = {"outFields": ",".join(FIELDS[layer]), "outSR": epsg, "f": "geojson"}
    return query_url(f"{source.service_url}/{LAYERS[layer]}/query", {**params, **fields})


def station_urls(source: StationSource, layer: int, numbers: Sequence[str]) -> list[str]:
    """Layer 0's or 38's queries for `numbers`, `CHUNK` stations each."""
    chunks = [numbers[k : k + CHUNK] for k in range(0, len(numbers), CHUNK)]
    where = ["stasjonnr in (" + ",".join(f"'{n}'" for n in c) + ")" for c in chunks]
    return [_url(source, layer, {"where": w}) for w in where]


def envelope_url(source: StationSource, layer: int, x: float, y: float, half: float) -> str:
    """Layer 2's or 5's query of the envelope `half` metres round `(x, y)`."""
    box = (x - half, y - half, x + half, y + half)
    epsg = source.crs.split(":")[-1]
    return _url(
        source,
        layer,
        {
            "where": "1=1",
            "geometry": ",".join(str(v) for v in box),
            "geometryType": "esriGeometryEnvelope",
            "inSR": epsg,
            "spatialRel": "esriSpatialRelIntersects",
        },
    )


def newest(versions: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    """The newest `oppdateringsdato` wins; a tie goes to the larger `objectid`."""
    return max(
        versions,
        key=lambda f: (f["properties"]["oppdateringsdato"] or 0, f["properties"]["objectid"]),
    )


def _features(text: str, where: str) -> list[dict[str, Any]]:
    doc = json.loads(text)
    if doc.get("exceededTransferLimit") or (doc.get("properties") or {}).get(
        "exceededTransferLimit"
    ):
        raise FetchError(f"{where}: the answer was truncated (exceededTransferLimit)")
    return list(doc["features"])


def _collection(crs: str, features: list[dict[str, Any]]) -> bytes:
    doc = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": crs}},
        "features": features,
    }
    return (json.dumps(doc, ensure_ascii=False, indent=1) + "\n").encode("utf-8")


def _feature(geometry: Any, properties: Mapping[str, Any]) -> dict[str, Any]:
    return {"type": "Feature", "geometry": geometry, "properties": dict(properties)}


def fetch_station_set(source: StationSource, get_text: Callable[[str], str]) -> dict[str, bytes]:
    """Each output file's name and bytes, `manifest.json` last. Refuses, with a
    `FetchError` naming the station, a listed station with no point or no
    polygon, and a river or lake answer flagged as truncated."""
    rows = read_list(source)
    numbers = [r["station"] for r in rows]
    asked: list[str] = []

    def get(url: str, where: str) -> list[dict[str, Any]]:
        asked.append(url)
        return _features(get_text(url), where)

    points: dict[str, dict[str, Any]] = {}
    for url in station_urls(source, 0, numbers):
        points |= {str(f["properties"]["stasjonnr"]): f for f in get(url, "station points")}
    polygons: dict[str, list[dict[str, Any]]] = {}
    for url in station_urls(source, 38, numbers):
        for f in get(url, "catchment polygons"):
            polygons.setdefault(str(f["properties"]["stasjonnr"]), []).append(f)
    for n in numbers:
        for found, what in ((points, "point"), (polygons, "polygon")):
            if n not in found:
                raise FetchError(f"station {n}: NVE serves no {what} for it")
    segments: dict[int, dict[str, Any]] = {}
    lakes: dict[int, dict[str, Any]] = {}
    for n in numbers:
        x, y = points[n]["geometry"]["coordinates"][:2]
        for f in get(envelope_url(source, 2, x, y, ENVELOPE_HALF), f"station {n}'s rivers"):
            segments.setdefault(int(f["properties"]["objectid"]), f)
        for f in get(envelope_url(source, 5, x, y, LAKE_ENVELOPE_HALF), f"station {n}'s lakes"):
            lakes.setdefault(int(f["properties"]["objectid"]), f)

    stations, references = [], []
    for row in rows:
        n, p = row["station"], points[row["station"]]["properties"]
        hierarchy = (p.get("elvenavnhierarki") or "").split("/")[0].strip()
        stations.append(
            _feature(
                points[n]["geometry"],
                {
                    "station": n,
                    "name": p.get("stasjonnavn"),
                    "series": [f"1001.{row['series_version']}"],
                    "nve_area_km2": p.get("totalt_feltareal_km2"),
                    "hrd_start_daily": int(row["hrd_start_daily"]),
                    "watercourse": p.get("vassdragsnr"),
                    "river": hierarchy or None,
                },
            )
        )
        chosen = newest(polygons[n])
        ms = chosen["properties"]["oppdateringsdato"]
        day = None if ms is None else datetime.fromtimestamp(ms / 1000, UTC).date().isoformat()
        references.append(
            _feature(
                chosen["geometry"],
                {
                    "station": n,
                    "reference_area_km2": chosen["properties"]["nedborfeltaareal_km2"],
                    "reference_updated": day,
                    "versions": len(polygons[n]),
                },
            )
        )
    rivers, lake_features = (
        [
            _feature(f["geometry"], {k: f["properties"].get(k) for k in FIELDS[layer]})
            for _, f in sorted(found.items())
        ]
        for layer, found in ((2, segments), (5, lakes))
    )
    files = dict(
        zip(
            FILES,
            (
                _collection(source.crs, stations),
                _collection(source.crs, references),
                _collection(source.crs, rivers),
                _collection(source.crs, lake_features),
                notice(source).encode("utf-8"),
            ),
            strict=True,
        )
    )
    manifest = {
        "source": source.id,
        "fetched": datetime.now(UTC).isoformat(timespec="seconds"),
        "requests": asked,
        "files": {name: hashlib.sha256(data).hexdigest() for name, data in files.items()},
    }
    files["manifest.json"] = (json.dumps(manifest, indent=1) + "\n").encode("utf-8")
    return files


__all__ = ["CHUNK", "ENVELOPE_HALF", "FIELDS", "FILES", "fetch_station_set", "read_list"]
