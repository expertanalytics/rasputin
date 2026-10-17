"""`rasputin fetch` of an `arcgis-features` source: a map-service layer's polygons meeting a box.

The layer is described first (`<url>?f=json`) and checked against the
catalogue entry; then the ids of the features meeting the box are asked
(`returnIdsOnly`, which the ArcGIS REST API answers with no record limit),
and the ids not cached are asked by `objectIds` in chunks, one request at a
time. A chunk the service refuses is halved until one id is refused alone:
a single answer may be too large for the server (two large sea polygons of
CORINE in one answer are an HTTP 500 every time; each alone is answered).

Every answer is checked: the ids returned are the ids asked, its CRS is the
one asked, and each geometry is a polygon or multipolygon with finite
coordinates. Polygons are stored as served: no repair, no clipping. Each
chunk is one SQLite transaction of `features.gpkg`, so a rerun resumes from
the ids present. The network call is injected (`get_text`), as in
`fetch/nve.py`; nothing here opens a file.
"""

from __future__ import annotations

import datetime
import json
import time
from collections.abc import Callable, Mapping, Sequence
from contextlib import closing
from typing import Any, Literal

import numpy as np
import shapely
from pydantic import BaseModel, ConfigDict
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from tin_engine import installed_version
from tin_engine.crs import same_crs, srs_entry
from tin_engine.fetch.http import FetchError, query_url
from tin_engine.fetch.plan import FetchRequest, box_in_crs, region_sha256
from tin_engine.io.geopackage_write import insert_features
from tin_engine.io.repository import FEATURES, CachedRequest, CacheWriter, FeatureManifest
from tin_engine.sources import RemoteSource, notice

#: The most ids asked in one request, a chosen value: NVE dropped a GET of
#: 2000 ids (10,906 characters) and answered 1000; 500 keeps NVE's and
#: CORINE's ids under the 8,000-character URL CORINE answered.
MAX_IDS = 500

GetText = Callable[[str], str]


class LayerFacts(BaseModel):
    """What a layer states of itself (`<url>?f=json`)."""

    model_config = ConfigDict(frozen=True)

    geometry: Literal["esriGeometryPolygon"]
    id_field: str
    #: Every field, name to esri field type.
    fields: Mapping[str, str]
    max_records: int
    epsg: int


class FeatureReport(BaseModel):
    """What a fetch did, or, for a dry run, what the real run will do."""

    model_config = ConfigDict(frozen=True)

    in_box: int
    present: int
    fetched: int
    requests: int
    bytes: int
    seconds: float


def _answer(get_text: GetText, url: str, where: str) -> Any:
    """The JSON of `url`. ArcGIS answers some refusals with HTTP 200 and an
    `error` member; that is a `FetchError` with the service's own message."""
    doc = json.loads(get_text(url))
    if isinstance(doc, dict) and "error" in doc:
        error = doc["error"] if isinstance(doc["error"], dict) else {"message": doc["error"]}
        said = "; ".join(str(m) for m in (error.get("message"), *(error.get("details") or ())) if m)
        raise FetchError(f"{where}: the service refused the request: {said or 'no reason given'}")
    return doc


def layer_facts(source: RemoteSource, get_text: GetText) -> LayerFacts:
    """The layer's facts, refused (a `FetchError` naming the source and URL)
    unless it is a polygon layer in the entry's CRS, holding every field the
    entry asks for, whose ids are not hash values (which a query by
    `objectIds` must avoid)."""
    url, crs = source.url_for(), source.crs
    assert crs is not None
    where = f"{source.id}: {url}"
    doc = _answer(get_text, query_url(url, {"f": "json"}), where)
    if doc.get("geometryType") != "esriGeometryPolygon":
        raise FetchError(
            f"{where}: a {doc.get('geometryType')} layer; only polygon layers are fetched"
        )
    unique = doc.get("uniqueIdInfo") or {}
    if unique.get("OIDFieldContainsHashValue") or unique.get("oidFieldContainsHashValue"):
        raise FetchError(f"{where}: its ids are hash values (OIDFieldContainsHashValue)")
    fields = {str(f["name"]): str(f["type"]) for f in doc.get("fields") or ()}
    ids = [n for n, t in fields.items() if t == "esriFieldTypeOID"]
    reference = (doc.get("extent") or {}).get("spatialReference") or {}
    epsg = reference.get("latestWkid") or reference.get("wkid")
    if len(ids) != 1 or epsg is None:
        raise FetchError(f"{where}: the layer states no single id field or no CRS")
    if not same_crs(f"EPSG:{epsg}", crs):
        raise FetchError(f"{where}: the layer's CRS is EPSG:{epsg}; the catalogue's is {crs}")
    if missing := [f for f in source.fields if f not in fields]:
        raise FetchError(f"{where}: the layer has no field {', '.join(missing)}")
    return LayerFacts(
        geometry="esriGeometryPolygon",
        id_field=ids[0],
        fields=fields,
        max_records=int(doc.get("maxRecordCount") or MAX_IDS),
        epsg=int(epsg),
    )


def chunks(ids: Sequence[int], size: int) -> list[list[int]]:
    """`ids` in order, in runs of at most `size`."""
    return [list(ids[k : k + size]) for k in range(0, len(ids), size)]


def _rows(
    doc: Any, asked: Sequence[int], facts: LayerFacts, fields: Sequence[str], where: str
) -> list[tuple[int, BaseGeometry, list[Any]]]:
    """A features answer's rows `(id, polygon, values of fields)`, checked."""
    stated = ((doc.get("crs") or {}).get("properties") or {}).get("name")
    if stated != f"EPSG:{facts.epsg}":
        raise FetchError(f"{where}: the answer is in {stated}, not EPSG:{facts.epsg} as asked")
    rows = []
    for feature in doc.get("features") or ():
        properties = feature.get("properties") or {}
        key = properties.get(facts.id_field, feature.get("id"))
        if not isinstance(key, int):
            raise FetchError(f"{where}: a feature without its id ({facts.id_field})")
        given = feature.get("geometry") or {}
        try:
            if given.get("type") not in ("Polygon", "MultiPolygon"):
                raise ValueError(f"a {given.get('type')} geometry")
            polygon = shape(given)
            if not np.isfinite(shapely.get_coordinates(polygon)).all():
                raise ValueError("a coordinate that is not finite")
        except (ValueError, TypeError, AttributeError) as exc:
            raise FetchError(f"{where}: feature {key} has no polygon: {exc}") from exc
        rows.append((key, polygon, [properties.get(f) for f in fields]))
    served = [k for k, _, _ in rows]
    if sorted(served) != sorted(asked):
        raise FetchError(f"{where}: the layer changed during the fetch; run again")
    return rows


def fetch_features(
    request: FetchRequest, source: RemoteSource, get_text: GetText, cache: CacheWriter
) -> FeatureReport:
    """Copy into `cache` the polygons of `source`'s layer meeting the
    request's box (or its domain's bounds), moved into the layer's CRS with
    no growth: a feature meeting the box comes whole."""
    t0, sent, got = time.perf_counter(), [0], [0]
    try:
        url = source.url_for(request.year)
    except ValueError as exc:
        raise FetchError(str(exc)) from exc
    where = f"{source.id}: {url}"

    def get(link: str) -> str:
        sent[0] += 1
        text = get_text(link)
        got[0] += len(text.encode())
        return text

    facts = layer_facts(source, get)
    epsg = str(facts.epsg)
    b = box_in_crs(request, f"EPSG:{facts.epsg}")
    query = f"{url}/query"
    found = _answer(
        get,
        query_url(
            query,
            {
                "where": "1=1",
                "geometry": ",".join(str(v) for v in (b.x_min, b.y_min, b.x_max, b.y_max)),
                "geometryType": "esriGeometryEnvelope",
                "inSR": epsg,
                "spatialRel": "esriSpatialRelIntersects",
                "returnIdsOnly": "true",
                "f": "json",
            },
        ),
        where,
    )
    if found.get("exceededTransferLimit"):
        raise FetchError(f"{where}: the ids answer was truncated (exceededTransferLimit)")
    in_box = {int(i) for i in found.get("objectIds") or ()}
    size = min(facts.max_records, MAX_IDS)

    def report(present: int, fetched: int, planned: int = 0) -> FeatureReport:
        return FeatureReport(
            in_box=len(in_box),
            present=present,
            fetched=fetched,
            requests=sent[0] + planned,
            bytes=got[0],
            seconds=time.perf_counter() - t0,
        )

    if request.dry_run:
        present = set() if request.refresh else in_box & cache.feature_ids(facts.id_field)
        return report(len(present), 0, len(chunks(sorted(in_box - present), size)))
    asked = CachedRequest(
        region_sha256=region_sha256(request, source.crs), date=datetime.date.today()
    )
    columns = [(f, facts.fields[f]) for f in source.fields]
    srs = srs_entry(f"EPSG:{facts.epsg}")
    with cache:
        if request.refresh:
            cache.discard()
        manifest = cache.feature_manifest() or FeatureManifest(
            source=source.id,
            url=url,
            crs=f"EPSG:{facts.epsg}",
            id_field=facts.id_field,
            fields=source.fields,
            rasputin_version=installed_version(),
            first_fetched=datetime.date.today(),
        )
        if (manifest.id_field, manifest.fields) != (facts.id_field, source.fields):
            raise FetchError(f"{where}: the cache holds other fields; --refresh fetches anew")
        with closing(cache.features(srs, facts.id_field, columns)) as conn:
            cache.put_manifest(manifest)
            cache.put_notice(
                notice(source, (source.credit_for(url=url, date=manifest.first_fetched),))
            )
            present = in_box & cache.feature_ids(facts.id_field)

            def take(ids: list[int]) -> None:
                link = query_url(
                    query,
                    {
                        "objectIds": ",".join(map(str, ids)),
                        "outFields": ",".join((facts.id_field, *source.fields)),
                        "outSR": epsg,
                        "returnGeometry": "true",
                        "f": "geojson",
                    },
                )
                try:
                    doc = _answer(get, link, where)
                except FetchError as exc:
                    if len(ids) == 1:
                        raise FetchError(f"{where}: feature {ids[0]} is refused: {exc}") from exc
                    for half in (ids[: len(ids) // 2], ids[len(ids) // 2 :]):
                        take(half)
                    return
                insert_features(
                    conn, FEATURES, srs.srs_id, _rows(doc, ids, facts, source.fields, where)
                )

            for chunk in chunks(sorted(in_box - present), size):
                take(chunk)
        if asked not in manifest.requests:
            cache.put_manifest(
                manifest.model_copy(update={"requests": (*manifest.requests, asked)})
            )
    return report(len(present), len(in_box - present))


__all__ = ["MAX_IDS", "FeatureReport", "LayerFacts", "chunks", "fetch_features", "layer_facts"]
