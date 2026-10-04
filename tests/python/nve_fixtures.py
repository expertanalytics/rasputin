"""A fake NVE map service and hand-written station, polygon and river files (29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "The station set" and "Data
use". No test touches the network: `FakeNve.get_text` stands in for
`RangeClient.get_text` (the one network call `fetch/nve.py` makes, injected
by `monkeypatch` on the class), answers from tables held in memory, and logs
every URL it is asked, so a test reads the requests back.

The fake answers as the real service does (checked against `kart.nve.no`,
2026-10-04): `f=geojson` gives a FeatureCollection with no `crs` member, each
feature with an `id` and a `properties` object; a truncated answer carries
`"exceededTransferLimit": true` both at the top level and in a top-level
`properties` object. It honours `outFields`, except that it **leaks** one
forbidden field per layer (`stasjoneier`, `qnormal6190_m3s`, `oppdatertav`)
in every answer, so a writer that copies whatever it is given is caught.

The tables are built from the packaged list (`tin_engine/data/nve_hrd_2025.csv`),
so the fetch sees the real 140 stations. Station `i` (list order) sits on a
20 km grid; the special cases sit at named stations:

- `NEWEST` (2.11.0): three polygon versions, the newest by date served
  second (areas 119.74, 119.43, 119.74; Narsjø's real ones);
- `TIE` (2.32.0): two versions with one date, the larger `objectid` served
  first; it sits 1500 m east of `NEWEST`, and segment `SHARED` lies between
  them, so both envelopes return it;
- `COPIES` (2.142.0): segments `COPY_LOW` and `COPY_HIGH` are one geometry
  under two `objectid`s and one `elvid` (served high first), and `SAME_NUMBER`
  shares `COPY_LOW`'s `strekninglnr` with another geometry (79.3.0's case);
- `MULTI` (19.79.0): a two-part polygon (Gravå), and the `objekttype` cases
  `NULL_TYPE` (null, `vatnlnr` set), `BLANK_TYPE` (" ", `vatnlnr` 0, which
  is how NVE says "no lake": both blank-type features of the real layer have it),
  `ODD_LAKE` ("InnsjoMidtlin") and `STRAY` ("SK");
- `NO_HIERARCHY` (2.265.0): `elvenavnhierarki` is null;
- `RENAMED` (311.4.0): layer 0 serves `stasjonnavn` as `LAYER_0_NAME`, a
  synthetic name no list row can hold, so a fetch that copies the list's
  name fails whatever the packaged list writes for this station (the real
  "Femundsenden (Femunden)" would pass if the list held it in full).
  Every other station's layer 0 name is its list row's.
"""

from __future__ import annotations

import csv
import importlib.resources
import io
import json
import re
import threading
import time
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlsplit

from shapely.geometry import box, mapping, shape

CRS = "EPSG:25833"
SERVICES = "https://kart.nve.no/enterprise/rest/services"
LAYER_0 = f"{SERVICES}/HydrologiskeData3/MapServer/0/query"
LAYER_38 = f"{SERVICES}/HydrologiskeData3/MapServer/38/query"
LAYER_2 = f"{SERVICES}/Elvenett1/MapServer/2/query"
ENDPOINTS = {LAYER_0: 0, LAYER_38: 38, LAYER_2: 2}

#: "Data use": the explicit allow-list per layer, never `outFields=*`.
ALLOW: Mapping[int, frozenset[str]] = {
    0: frozenset(
        {
            "stasjonnr",
            "stasjonnavn",
            "totalt_feltareal_km2",
            "stasjonstatus",
            "vassdragsnr",
            "elvenavnhierarki",
        }
    ),
    38: frozenset({"stasjonnr", "nedborfeltaareal_km2", "oppdateringsdato", "objectid"}),
    2: frozenset(
        {"objectid", "objekttype", "strekninglnr", "elvid", "vassdragsnr", "elvenavn", "vatnlnr"}
    ),
}
#: Never collected ("Data use"); `qnormal` prefixes layer 38's discharge normals.
FORBIDDEN = ("stasjoneier", "oppdatertav", "globalid", "qnormal")
OWNER = "Eier Kraftverk AS"
LEAKED: Mapping[int, Mapping[str, Any]] = {
    0: {"stasjoneier": OWNER, "globalid": "{0F0F}"},
    38: {"qnormal6190_m3s": 6.12, "globalid": "{1E1E}"},
    2: {"oppdatertav": "ABC", "globalid": "{2D2D}"},
}

#: `map_radius + reach_up + 500 m` each way ("The station set").
ENVELOPE_HALF = 2000.0
CHUNK = 40

NEWEST, TIE, COPIES, MULTI, NO_HIERARCHY = "2.11.0", "2.32.0", "2.142.0", "19.79.0", "2.265.0"
#: Ola's ruling of 2026-10-04: the station's name is layer 0's `stasjonnavn`.
#: Synthetic, so the test does not depend on how the packaged list spells it.
RENAMED, LAYER_0_NAME = "311.4.0", "Layer-0 name of 311.4.0"
SHARED = 900001
COPY_LOW, COPY_HIGH, SAME_NUMBER = 900010, 900011, 900012
NULL_TYPE, BLANK_TYPE, ODD_LAKE, STRAY = 900020, 900021, 900022, 900023
#: Epoch milliseconds, as layer 38 serves `oppdateringsdato`.
MS_2023_04_11, MS_2025_03_19, MS_2026_05_28 = 1681171200000, 1742342400000, 1779926400000

STATION_NUMBER = re.compile(r"^\d+\.\d+\.\d+$")


# --------------------------------------------------------------------------
# The packaged list, read as data (the test's own reader, not the package's)
# --------------------------------------------------------------------------

LIST_COLUMNS = ["station", "series_version", "name", "hrd_start_daily"]


def list_text() -> str:
    """The packaged HRD list as text, through `importlib.resources`."""
    resource = importlib.resources.files("tin_engine") / "data" / "nve_hrd_2025.csv"
    return resource.read_text(encoding="utf-8")


def list_comments() -> list[str]:
    return [line for line in list_text().splitlines() if line.startswith("#")]


def list_rows() -> list[dict[str, str]]:
    """The data rows, in file order: comment lines (`#`) skipped, then a
    header naming the four columns."""
    body = "\n".join(line for line in list_text().splitlines() if not line.startswith("#"))
    reader = csv.DictReader(io.StringIO(body))
    assert reader.fieldnames == LIST_COLUMNS, reader.fieldnames
    return list(reader)


# --------------------------------------------------------------------------
# The fake service
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Call:
    """One request: its endpoint's layer and its query parameters (last value each)."""

    url: str
    layer: int
    params: Mapping[str, str]

    @property
    def out_fields(self) -> list[str]:
        return [f.strip() for f in self.params.get("outFields", "").split(",") if f.strip()]

    @property
    def stations(self) -> list[str]:
        """The station numbers of a `where ... in (...)` clause, in order."""
        return re.findall(r"'(\d+\.\d+\.\d+)'", self.params.get("where", ""))

    @property
    def envelope(self) -> tuple[float, float, float, float]:
        """`geometry` as `(xmin, ymin, xmax, ymax)`, given either as four
        comma-separated numbers or as an ArcGIS envelope object."""
        text = self.params["geometry"]
        if text.lstrip().startswith("{"):
            env = json.loads(text)
            return (env["xmin"], env["ymin"], env["xmax"], env["ymax"])
        x0, y0, x1, y1 = (float(v) for v in text.split(","))
        return (x0, y0, x1, y1)


@dataclass
class Station:
    """One station as layer 0 serves it (`name` is `stasjonnavn`)."""

    number: str
    name: str
    x: float
    y: float
    area_km2: float
    watercourse: str
    hierarchy: str | None


@dataclass
class FakeNve:
    """The three layers in memory; `get_text(url)` answers one query."""

    stations: list[Station]
    polygons: list[dict[str, Any]]
    segments: list[dict[str, Any]]
    calls: list[Call] = field(default_factory=list)
    unexpected: list[str] = field(default_factory=list)
    drop_point: set[str] = field(default_factory=set)
    drop_polygon: set[str] = field(default_factory=set)
    truncate: set[str] = field(default_factory=set)
    #: `(layer, field)` pairs the fake leaves out of every answer of that layer.
    omit: set[tuple[int, str]] = field(default_factory=set)
    max_active: int = 0
    _active: int = 0
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def station(self, number: str) -> Station:
        return next(s for s in self.stations if s.number == number)

    def get_text(self, url: str) -> str:
        with self._lock:
            self._active += 1
            self.max_active = max(self.max_active, self._active)
        try:
            time.sleep(0.0005)  # widens any overlap of two calls; adds nothing when sequential
            return self._answer(url)
        finally:
            with self._lock:
                self._active -= 1

    def calls_to(self, layer: int) -> list[Call]:
        return [c for c in self.calls if c.layer == layer]

    def _answer(self, url: str) -> str:
        parts = urlsplit(url)
        base = f"{parts.scheme}://{parts.netloc}{parts.path}"
        if base not in ENDPOINTS:
            self.unexpected.append(url)
            from tin_engine.fetch.http import FetchError

            raise FetchError(f"{url}: HTTP 404 Not Found (not a layer the design names)")
        params = {k: v[-1] for k, v in parse_qs(parts.query, keep_blank_values=True).items()}
        call = Call(url, ENDPOINTS[base], params)
        self.calls.append(call)
        wanted = set(call.out_fields)
        if call.layer == 0:
            features = [self._point(s, wanted) for s in self._listed(call, self.drop_point)]
        elif call.layer == 38:
            numbers = {s.number for s in self._listed(call, self.drop_polygon)}
            features = [_project(p, wanted, 38) for p in self.polygons if _number(p) in numbers]
        else:
            area = box(*call.envelope)
            features = [
                _project(s, wanted, 2)
                for s in self.segments
                if shape(s["geometry"]).intersects(area)
            ]
        omitted = {name for layer, name in self.omit if layer == call.layer}
        for f in features:
            f["properties"] = {k: v for k, v in f["properties"].items() if k not in omitted}
        doc: dict[str, Any] = {"type": "FeatureCollection", "features": features}
        if call.layer == 2 and any(area.contains(_xy(self.station(n))) for n in self.truncate):
            doc["exceededTransferLimit"] = True
            doc["properties"] = {"exceededTransferLimit": True}
        return json.dumps(doc)

    def _listed(self, call: Call, dropped: set[str]) -> list[Station]:
        asked = set(call.stations)
        return [s for s in self.stations if s.number in asked and s.number not in dropped]

    def _point(self, s: Station, wanted: set[str]) -> dict[str, Any]:
        full = {
            "stasjonnr": s.number,
            "stasjonnavn": s.name,
            "totalt_feltareal_km2": s.area_km2,
            "stasjonstatus": 1,
            "vassdragsnr": s.watercourse,
            "elvenavnhierarki": s.hierarchy,
            "objectid": 10_000 + self.stations.index(s),
        }
        feature = {
            "type": "Feature",
            "id": full["objectid"],
            "geometry": {"type": "Point", "coordinates": [s.x, s.y]},
            "properties": full,
        }
        return _project(feature, wanted, 0)


def _number(feature: Mapping[str, Any]) -> str:
    return str(feature["properties"]["stasjonnr"])


def _xy(s: Station) -> Any:
    from shapely.geometry import Point

    return Point(s.x, s.y)


def _project(feature: Mapping[str, Any], wanted: set[str], layer: int) -> dict[str, Any]:
    """The feature with the asked fields only, plus the layer's leaked ones."""
    props = {k: v for k, v in feature["properties"].items() if k in wanted}
    return {**feature, "properties": {**props, **LEAKED[layer]}}


def build_fake() -> FakeNve:
    """The fake service for the packaged list (see the module docstring)."""
    rows = list_rows()
    stations: list[Station] = []
    for i, row in enumerate(rows):
        regine = int(row["station"].split(".")[0])
        x, y = 300_000.0 + 20_000.0 * (i % 10), 6_600_000.0 + 20_000.0 * (i // 10)
        hierarchy = None if row["station"] == NO_HIERARCHY else f"Elv{i}/Hovedvassdrag{i}"
        stations.append(
            Station(row["station"], row["name"], x, y, 10.0 + i, f"{regine:03d}.A{i}", hierarchy)
        )
    by = {s.number: s for s in stations}
    by[RENAMED].name = LAYER_0_NAME
    newest, tie = by[NEWEST], by[TIE]
    tie.x, tie.y = newest.x + 1500.0, newest.y

    polygons = [
        _polygon(s, 3000 + i, MS_2025_03_19, 4.0 + i / 100)
        for i, s in enumerate(stations)
        if s.number not in (NEWEST, TIE, MULTI)
    ]
    polygons += [
        _polygon(newest, 7001, MS_2023_04_11, 119.74),
        _polygon(newest, 7002, MS_2026_05_28, 119.43),
        _polygon(newest, 7003, MS_2025_03_19, 119.74),
        _polygon(tie, 7102, MS_2026_05_28, 51.0),
        _polygon(tie, 7101, MS_2026_05_28, 50.0),
        _polygon(by[MULTI], 7201, MS_2025_03_19, 8.0, parts=2),
    ]

    segments = [_own_segment(s, 100_000 + i) for i, s in enumerate(stations)]
    mx, my = newest.x + 750.0, newest.y
    segments.append(_segment(SHARED, [(mx, my - 300), (mx, my + 300)], elvid="2-11-77"))
    c = by[COPIES]
    copy_line = [(c.x - 200, c.y + 60), (c.x, c.y + 70), (c.x + 200, c.y + 80)]
    segments.append(_segment(COPY_HIGH, copy_line, elvid="2-142-9", strekning="2142001"))
    segments.append(_segment(COPY_LOW, copy_line, elvid="2-142-9", strekning="2142001"))
    segments.append(
        _segment(
            SAME_NUMBER,
            [(c.x - 200, c.y - 400), (c.x + 200, c.y - 380)],
            elvid="2-142-9",
            strekning="2142001",
            objekttype="InnsjøMidtlinje",
            vatnlnr=77,
        )
    )
    m = by[MULTI]
    for k, (oid, kind, lake) in enumerate(
        [
            (NULL_TYPE, None, 495),
            (BLANK_TYPE, " ", 0),
            (ODD_LAKE, "InnsjoMidtlin", 12),
            (STRAY, "SK", None),
        ]
    ):
        line = [(m.x - 300, m.y - 100 * (k + 1)), (m.x + 300, m.y - 100 * (k + 1))]
        segments.append(_segment(oid, line, elvid=f"19-79-{k}", objekttype=kind, vatnlnr=lake))
    return FakeNve(stations, polygons, segments)


def _square(x: float, y: float, half: float) -> list[list[float]]:
    return [
        [x - half, y - half],
        [x + half, y - half],
        [x + half, y + half],
        [x - half, y + half],
        [x - half, y - half],
    ]


def _polygon(s: Station, oid: int, ms: int, area: float, parts: int = 1) -> dict[str, Any]:
    if parts == 1:
        geometry: dict[str, Any] = {"type": "Polygon", "coordinates": [_square(s.x, s.y, 1000.0)]}
    else:
        geometry = {
            "type": "MultiPolygon",
            "coordinates": [[_square(s.x, s.y, 800.0)], [_square(s.x + 3000.0, s.y, 300.0)]],
        }
    props = {
        "stasjonnr": s.number,
        "nedborfeltaareal_km2": area,
        "oppdateringsdato": ms,
        "objectid": oid,
    }
    return {"type": "Feature", "id": oid, "geometry": geometry, "properties": props}


def _own_segment(s: Station, oid: int) -> dict[str, Any]:
    line = [(s.x - 300.0, s.y + 20.0), (s.x + 300.0, s.y + 20.0)]
    return _segment(
        oid, line, elvid=f"own-{oid}", vassdragsnr=s.watercourse, name=f"Elv {s.number}"
    )


def _segment(
    oid: int,
    line: Sequence[tuple[float, float]],
    *,
    elvid: str,
    objekttype: str | None = "ElvBekk",
    vatnlnr: int | None = None,
    strekning: str | None = None,
    vassdragsnr: str | None = "000.X",
    name: str | None = None,
) -> dict[str, Any]:
    props = {
        "objectid": oid,
        "objekttype": objekttype,
        "strekninglnr": strekning or str(oid),
        "elvid": elvid,
        "vassdragsnr": vassdragsnr,
        "elvenavn": name,
        "vatnlnr": vatnlnr,
    }
    geometry = {"type": "LineString", "coordinates": [list(p) for p in line]}
    return {"type": "Feature", "id": oid, "geometry": geometry, "properties": props}


# --------------------------------------------------------------------------
# Hand-written files for the readers
# --------------------------------------------------------------------------


def crs_member(name: str = CRS) -> dict[str, Any]:
    return {"type": "name", "properties": {"name": name}}


def collection(features: Iterable[Mapping[str, Any]], crs: str | None = CRS) -> dict[str, Any]:
    doc: dict[str, Any] = {"type": "FeatureCollection", "features": list(features)}
    if crs is not None:
        doc["crs"] = crs_member(crs)
    return doc


def feature(geometry: Any, **properties: Any) -> dict[str, Any]:
    """A GeoJSON feature from a shapely geometry or a GeoJSON geometry dict."""
    geo = geometry if isinstance(geometry, dict) else mapping(geometry)
    return {"type": "Feature", "geometry": json.loads(json.dumps(geo)), "properties": properties}


def point(x: float, y: float, **properties: Any) -> dict[str, Any]:
    return feature({"type": "Point", "coordinates": [x, y]}, **properties)


def line(coords: Sequence[tuple[float, float]], **properties: Any) -> dict[str, Any]:
    return feature({"type": "LineString", "coordinates": [list(c) for c in coords]}, **properties)


def river(
    objectid: int,
    coords: Sequence[tuple[float, float]],
    *,
    elvid: str = "1-1-1",
    objekttype: str | None = "ElvBekk",
    vatnlnr: int | None = None,
    strekninglnr: str | None = None,
    vassdragsnr: str | None = "002.A",
    elvenavn: str | None = "Elva",
) -> dict[str, Any]:
    """One `rivers.geojson` feature with the seven fields `fetch-stations` writes."""
    return line(
        coords,
        objectid=objectid,
        objekttype=objekttype,
        strekninglnr=strekninglnr or str(objectid),
        elvid=elvid,
        vassdragsnr=vassdragsnr,
        elvenavn=elvenavn,
        vatnlnr=vatnlnr,
    )


def write(path: Path, doc: Mapping[str, Any]) -> Path:
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path
