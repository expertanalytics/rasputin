"""A fake ArcGIS map-service layer and two hand-made layers of four polygons. Test support only.

`docs/increments/63-fetch-vectors.md`, *6. Tests*. No test reaches the network:
`FakeLayer.get_text` stands in for `RangeClient.get_text`, answers from memory
and logs every URL it is asked, so a test reads the requests back.

The fake answers in the shapes measured from `kart.nve.no` and
`image.discomap.eea.europa.eu` (2026-10-17, `curl`):

- `<url>?f=json`: the layer, with `geometryType`, `fields` (name and esri
  type), `maxRecordCount` and `extent.spatialReference` (`wkid`, `latestWkid`).
  Neither real layer states `uniqueIdInfo`; the fake states it only when a
  test asks for a hashed id field, and then under both spellings
  (`OIDFieldContainsHashValue`, as the API's caution names it, and
  `oidFieldContainsHashValue`), since no answer carrying it was measured.
- `query?...&returnIdsOnly=true&f=json`: `{"objectIdFieldName", "objectIds"}`,
  the ids in no particular order (the fake gives them descending), and
  `"objectIds": null` for a box meeting nothing (both services do).
- `query?objectIds=...&outFields=...&outSR=...&f=geojson`: a FeatureCollection
  with a `crs` member `{"type": "name", "properties": {"name": "EPSG:<outSR>"}}`.
  A feature carries `id` and its id field among `properties` **only when the id
  field is in `outFields`** (both services leave both out otherwise). The fake
  serves features in descending id order, as NVE served two ids out of the
  order asked. It honours `outFields` exactly and never reprojects: it serves
  only `outSR` equal to the layer's own CRS.

A refusal is what `RangeClient.get_text` raises after its retries: a
`FetchError`. `refuse(ids)` decides which chunk requests are refused.

The two layers (ids 3, 7, 11, 12 meet `BOX_ALL`; id 20 lies far outside every
box). Every coordinate is an integer, so the GeoPackage R-tree's 32-bit floats
hold it exactly (all are below 2**24):

- id 3: a 100 m square; id 7: a 200 m square with a 100 m square hole;
- id 11: a bow-tie (self-intersecting ring), stored as served, unrepaired;
- id 12: a MultiPolygon of two 100 m squares 300 m apart;
- id 20: a 100 m square, 100 km away.

`NVE` sits in EPSG:25833 near (200 000, 6 600 000), as NVE's lake layer;
`CORINE` the same shapes in EPSG:3035 near (4 300 000, 4 000 000), with CORINE
codes. Each serves fields the catalogue entry does not ask for (NVE's
`objekttype`, `hoyde`, `magasinnr`; CORINE's `Remark`, `Area_Ha`), so a cache
that keeps whatever it is given is caught.
"""

from __future__ import annotations

import json
import urllib.parse
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import shapely
from shapely.geometry import MultiPolygon, Polygon, box, mapping

from tin_engine.fetch.http import FetchError

NVE_URL = "https://kart.nve.no/enterprise/rest/services/Innsjodatabase2/MapServer/5"
CORINE_URL = (
    "https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2018_LAEA/MapServer/0"
)
IDS = (3, 7, 11, 12)
FAR = 20


def _shapes(x0: int, y0: int) -> dict[int, Polygon | MultiPolygon]:
    """The five shapes, laid out from the corner (`x0`, `y0`)."""
    square = box(x0, y0, x0 + 100, y0 + 100)
    holed = Polygon(
        [(x0 + 1000, y0), (x0 + 1200, y0), (x0 + 1200, y0 + 200), (x0 + 1000, y0 + 200)],
        [
            [
                (x0 + 1050, y0 + 50),
                (x0 + 1150, y0 + 50),
                (x0 + 1150, y0 + 150),
                (x0 + 1050, y0 + 150),
            ]
        ],
    )
    bowtie = Polygon(
        [(x0 + 2000, y0), (x0 + 2100, y0 + 100), (x0 + 2100, y0), (x0 + 2000, y0 + 100)]
    )
    pair = MultiPolygon(
        [box(x0 + 3000, y0, x0 + 3100, y0 + 100), box(x0 + 3300, y0, x0 + 3400, y0 + 100)]
    )
    far = box(x0 + 100_000, y0 + 100_000, x0 + 100_100, y0 + 100_100)
    return {3: square, 7: holed, 11: bowtie, 12: pair, FAR: far}


@dataclass(frozen=True)
class LayerData:
    """One layer as the service holds it."""

    url: str
    epsg: int
    id_field: str
    #: Every field of the layer, name to esri type, in the layer's order.
    fields: Mapping[str, str]
    shapes: Mapping[int, Polygon | MultiPolygon]
    properties: Mapping[int, Mapping[str, Any]]
    #: Boxes in the layer's CRS, x0 y0 x1 y1: all four ids, the first two,
    #: the last three (id 7 in both halves).
    box_all: tuple[int, int, int, int]
    box_left: tuple[int, int, int, int]
    box_right: tuple[int, int, int, int]
    #: The fields the catalogue entry asks for.
    asked: tuple[str, ...]

    def extent(self, ids: Iterable[int]) -> tuple[float, float, float, float]:
        """The bounds of the union of the shapes `ids`, min_x min_y max_x max_y."""
        x0, y0, x1, y1 = shapely.total_bounds([self.shapes[i] for i in ids])
        return float(x0), float(y0), float(x1), float(y1)


def _boxes(x0: int, y0: int) -> dict[str, tuple[int, int, int, int]]:
    return {
        "box_all": (x0 - 1000, y0 - 1000, x0 + 4000, y0 + 1000),
        "box_left": (x0 - 1000, y0 - 1000, x0 + 1500, y0 + 1000),
        "box_right": (x0 + 1100, y0 - 1000, x0 + 4000, y0 + 1000),
    }


_NX, _NY = 200_000, 6_600_000
NVE = LayerData(
    NVE_URL,
    25833,
    "objectid",
    {
        "objectid": "esriFieldTypeOID",
        "shape": "esriFieldTypeGeometry",
        "objekttype": "esriFieldTypeString",
        "vatnlnr": "esriFieldTypeInteger",
        "navn": "esriFieldTypeString",
        "hoyde": "esriFieldTypeInteger",
        "areal_km2": "esriFieldTypeDouble",
        "magasinnr": "esriFieldTypeInteger",
    },
    _shapes(_NX, _NY),
    {
        i: {
            "objekttype": "Innsjo",
            "vatnlnr": 1000 + i,
            "navn": None if i == 7 else f"Vatn {i}",
            "hoyde": 400 + i,
            "areal_km2": i / 100,
            "magasinnr": None,
        }
        for i in (*IDS, FAR)
    },
    **_boxes(_NX, _NY),
    asked=("vatnlnr", "navn", "areal_km2"),
)

_CX, _CY = 4_300_000, 4_000_000
CODES = {3: "311", 7: "512", 11: "211", 12: "112", FAR: "231"}
CORINE = LayerData(
    CORINE_URL,
    3035,
    "OBJECTID",
    {
        "OBJECTID": "esriFieldTypeOID",
        "Shape": "esriFieldTypeGeometry",
        "Code_18": "esriFieldTypeString",
        "ID": "esriFieldTypeString",
        "Remark": "esriFieldTypeString",
        "Area_Ha": "esriFieldTypeDouble",
        "Shape_Length": "esriFieldTypeDouble",
        "Shape_Area": "esriFieldTypeDouble",
    },
    _shapes(_CX, _CY),
    {
        i: {
            "Code_18": CODES[i],
            "ID": f"EU-{i}",
            "Remark": None,
            "Area_Ha": 1.0 + i,
            "Shape_Length": 400.0,
            "Shape_Area": 10_000.0,
        }
        for i in (*IDS, FAR)
    },
    **_boxes(_CX, _CY),
    asked=("ID", "Code_18"),
)


#: An ArcGIS refusal as it comes, in an HTTP 200 body; its words are the service's.
REFUSAL: Mapping[str, Any] = {
    "error": {
        "code": 400,
        "message": "Unable to complete operation.",
        "details": ["Unable to perform query operation."],
    }
}


#: A maintenance page as a proxy serves it, in an HTTP 200 body: not JSON at all.
MAINTENANCE = "<!DOCTYPE html><html><body><h1>Down for maintenance</h1></body></html>"


@dataclass
class FakeLayer:
    """The service of one layer; every knob below is off by default."""

    data: LayerData
    #: The layer's `maxRecordCount`: the chunk size when below 500.
    max_records: int = 2000
    geometry_type: str = "esriGeometryPolygon"
    #: The CRS the layer states; None is the layer's own.
    stated_epsg: int | None = None
    #: State only `wkid`, no `latestWkid`.
    wkid_only: bool = False
    hashed_ids: bool = False
    #: Fields left out of the layer description.
    missing_fields: tuple[str, ...] = ()
    #: A chunk request whose ids make this true is refused with a `FetchError`.
    refuse: Callable[[list[int]], bool] = field(default=lambda ids: False)
    #: Ids left out of, or added to, every features answer.
    drop: frozenset[int] = frozenset()
    extra: frozenset[int] = frozenset()
    #: The `crs` member's name; None is `EPSG:<outSR>`.
    crs_name: str | None = None
    #: Geometry members served instead of the shape's, by id (None: JSON null).
    raw_geometry: Mapping[int, Any] = field(default_factory=dict)
    #: Flag the ids answer `exceededTransferLimit`.
    exceeded: bool = False
    #: Answer the ids query with `REFUSAL`, HTTP 200, as ArcGIS refuses some queries.
    ids_refused: bool = False
    #: Answer the chunks `refuse` picks with `REFUSAL`, HTTP 200, instead of an HTTP error.
    refuse_in_body: bool = False
    #: Answer every chunk request with `MAINTENANCE`, HTTP 200, as a proxy in front may.
    html_chunks: bool = False
    urls: list[str] = field(default_factory=list)
    served: int = 0

    def get_text(self, url: str) -> str:
        self.urls.append(url)
        parts = urllib.parse.urlsplit(url)
        base = urllib.parse.urlunsplit((parts.scheme, parts.netloc, parts.path, "", ""))
        query = {k: v[-1] for k, v in urllib.parse.parse_qs(parts.query).items()}
        if base == self.data.url:
            assert query.get("f") == "json", url
            text = json.dumps(self._layer())
        elif base == self.data.url + "/query" and query.get("returnIdsOnly") == "true":
            text = json.dumps(REFUSAL if self.ids_refused else self._ids(query, url))
        elif base == self.data.url + "/query" and "objectIds" in query:
            ids = [int(v) for v in query["objectIds"].split(",")]
            if self.html_chunks:
                self.served += len(MAINTENANCE.encode())
                return MAINTENANCE
            if self.refuse(ids) and not self.refuse_in_body:
                raise FetchError(f"{url}: HTTP 500, after 4 tries")
            refused = self.refuse(ids) and self.refuse_in_body
            text = json.dumps(REFUSAL if refused else self._features(query, ids, url))
        else:
            raise AssertionError(f"the fake serves no such request: {url}")
        self.served += len(text.encode())
        return text

    def _layer(self) -> dict[str, Any]:
        epsg = self.stated_epsg or self.data.epsg
        reference = {"wkid": epsg} if self.wkid_only else {"wkid": epsg, "latestWkid": epsg}
        x0, y0, x1, y1 = self.data.extent((*IDS, FAR))
        doc: dict[str, Any] = {
            "id": 0,
            "name": "layer",
            "type": "Feature Layer",
            "geometryType": self.geometry_type,
            "fields": [
                {"name": n, "type": t, "alias": n}
                for n, t in self.data.fields.items()
                if n not in self.missing_fields
            ],
            "maxRecordCount": self.max_records,
            "supportedQueryFormats": "JSON, geoJSON, PBF",
            "extent": {
                "xmin": x0, "ymin": y0, "xmax": x1, "ymax": y1, "spatialReference": reference
            },
        }  # fmt: skip
        if self.hashed_ids:
            doc["uniqueIdInfo"] = {
                "type": "simple",
                "fields": [self.data.id_field],
                "OIDFieldContainsHashValue": True,
                "oidFieldContainsHashValue": True,
            }
        return doc

    def _ids(self, query: Mapping[str, str], url: str) -> dict[str, Any]:
        assert query.get("f") == "json", url
        assert query.get("geometryType") == "esriGeometryEnvelope", url
        assert query.get("spatialRel") == "esriSpatialRelIntersects", url
        assert query.get("inSR") == str(self.data.epsg), url
        envelope = box(*(float(v) for v in query["geometry"].split(",")))
        found = sorted(
            (i for i, g in self.data.shapes.items() if g.intersects(envelope)), reverse=True
        )
        doc: dict[str, Any] = {"objectIdFieldName": self.data.id_field, "objectIds": found or None}
        if self.exceeded:
            doc["exceededTransferLimit"] = True
        return doc

    def _features(self, query: Mapping[str, str], ids: list[int], url: str) -> dict[str, Any]:
        assert query.get("f") == "geojson", url
        assert query.get("outSR") == str(self.data.epsg), url
        fields = query.get("outFields", "").split(",")
        keyed = self.data.id_field in fields
        served = sorted(({*ids} - self.drop) | self.extra, reverse=True)
        features = []
        for i in served:
            given = {self.data.id_field: i, **self.data.properties[i]}
            feature: dict[str, Any] = {"type": "Feature"}
            if keyed:
                feature["id"] = i
            if query.get("returnGeometry") != "true":
                feature["geometry"] = None
            elif i in self.raw_geometry:
                feature["geometry"] = self.raw_geometry[i]
            else:
                feature["geometry"] = mapping(self.data.shapes[i])
            feature["properties"] = {k: given[k] for k in fields if k in given}
            features.append(feature)
        name = self.crs_name or f"EPSG:{query['outSR']}"
        return {
            "type": "FeatureCollection",
            "crs": {"type": "name", "properties": {"name": name}},
            "features": features,
        }

    # ------------------------------------------------------------ reading the log

    def queries(self) -> list[dict[str, str]]:
        """Every request's query parameters, in the order asked."""
        return [
            {k: v[-1] for k, v in urllib.parse.parse_qs(urllib.parse.urlsplit(u).query).items()}
            for u in self.urls
        ]

    def chunk_queries(self) -> list[dict[str, str]]:
        """The feature requests' parameters (those naming `objectIds`)."""
        return [q for q in self.queries() if "objectIds" in q]

    def chunks_asked(self) -> list[list[int]]:
        """The ids of each feature request, as asked."""
        return [[int(v) for v in q["objectIds"].split(",")] for q in self.chunk_queries()]

    def ids_queries(self) -> list[dict[str, str]]:
        return [q for q in self.queries() if q.get("returnIdsOnly") == "true"]


def holds(*wanted: int) -> Callable[[list[int]], bool]:
    """A refusal of every chunk holding all of `wanted`."""
    return lambda ids: set(wanted) <= set(ids)
