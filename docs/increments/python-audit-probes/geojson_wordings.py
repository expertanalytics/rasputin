"""What every GeoJSON reader says for every malformed geometry shape
(python-audit.md, section 10, "Refusal wordings that change").

Run against one revision's `tin_engine`, imported from a scratch copy:

    python3 tools/scratch_copy.py <rev> <dir>        # prints PYTHONPATH=...
    PYTHONPATH=<that value> .venv/bin/python \\
        docs/increments/python-audit-probes/geojson_wordings.py <dir>

The scratch copy's sitecustomize drops the editable finder, and the probe
asserts that `tin_engine` was imported from `<dir>/src_python`, so a run
cannot measure the work tree by mistake. The probe writes its inputs to a
temporary directory and removes it.

One line per (reader, shape): `reader | shape | outcome`. The outcome is
`read: <summary>`, `refused: <message>` (an exception the command line
catches for that reader, so the user sees the message), or
`UNCAUGHT <class>: <message>` (a traceback from the command line). What the
command line catches: `--domain` a DomainError, `--features` and
`catchment --lakes` a FeatureError, the station, reference, NVE lake and
river readers an OSError or a ValueError. The temporary directory is
printed as `<dir>`. The full list of changed wordings between two revisions
is `diff` of their outputs.

Shapes: each bad geometry value (missing, null, the empty values `""`, 0,
false, [], {}, the non-empty non-objects 7, "Point", [1, 2], an object with
no `type`) and each empty `type` value (null, "", 0, false, [], {}) in a
geometry otherwise right for the reader, as a top-level Feature, as the
first feature of a collection and as the second (after a good one); a bare
geometry of the right type and of a wrong one; an object with only a `crs`
member; and a good Feature and a good collection, as controls. All of these
carry a `crs` member naming EPSG:25833. Then the file-level shapes, around
one good feature: each `crs` member rule (null, `{}`, a name that is null or
missing, an unreadable name, no member at all), a collection with no `type`,
`"features": {}`, `"features": [7]`, a JSON list, a file that is not JSON,
and a good collection behind a UTF-8 byte order mark.
"""

from __future__ import annotations

import json
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

CRS_MEMBER = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::25833"}}
X, Y = 500_000.0, 6_600_000.0
SQUARE = [[[X - 100, Y - 100], [X + 100, Y - 100], [X + 100, Y + 100], [X - 100, Y + 100],
           [X - 100, Y - 100]]]  # fmt: skip
GOOD = {"Point": [X, Y], "LineString": [[X - 50, Y], [X + 50, Y]], "Polygon": SQUARE}
MISSING = object()
BAD_GEOMETRY: list[tuple[str, Any]] = [
    ("geometry missing", MISSING), ("geometry null", None), ('geometry ""', ""),
    ("geometry 0", 0), ("geometry false", False), ("geometry []", []), ("geometry {}", {}),
    ("geometry 7", 7), ('geometry "Point"', "Point"), ("geometry [1, 2]", [1, 2]),
]  # fmt: skip
EMPTY_TYPES: list[tuple[str, Any]] = [
    ("null", None), ('""', ""), ("0", 0), ("false", False), ("[]", []), ("{}", {}),
]  # fmt: skip


def geometries(right: str) -> list[tuple[str, Any]]:
    """Every bad geometry value, for a reader whose right type is `right`."""
    coords = GOOD[right]
    typed = [(f"type {n}", {"type": t, "coordinates": coords}) for n, t in EMPTY_TYPES]
    return [*BAD_GEOMETRY, ("geometry without type", {"coordinates": coords}), *typed]


def feature(geometry: Any, k: int, prop: str) -> dict[str, Any]:
    props = {"station": f"1.2.{k}", "objectid": k + 1, "property": prop}
    f: dict[str, Any] = {"type": "Feature", "properties": props}
    if geometry is not MISSING:
        f["geometry"] = geometry
    return f


def shapes(right: str, wrong: str, prop: str) -> list[tuple[str, Any]]:
    """(shape name, the JSON document, or the file's bytes) for a reader whose
    right type is `right`."""
    good = feature({"type": right, "coordinates": GOOD[right]}, 0, prop)

    def collection(*features: dict[str, Any]) -> dict[str, Any]:
        return {"type": "FeatureCollection", "crs": CRS_MEMBER, "features": list(features)}

    out: list[tuple[str, Any]] = []
    for name, g in geometries(right):
        out.append((f"Feature, {name}", {**feature(g, 0, prop), "crs": CRS_MEMBER}))
        out.append((f"first feature, {name}", collection(feature(g, 0, prop))))
        out.append((f"second feature, {name}", collection(good, feature(g, 1, prop))))
    out.append(("bare geometry, right type", {"type": right, "coordinates": GOOD[right],
                                              "crs": CRS_MEMBER}))  # fmt: skip
    out.append(("bare geometry, wrong type", {"type": wrong, "coordinates": GOOD[wrong],
                                              "crs": CRS_MEMBER}))  # fmt: skip
    out.append(("crs member only", {"crs": CRS_MEMBER}))
    out.append(("control: good Feature", {**good, "crs": CRS_MEMBER}))
    out.append(("control: good collection", collection(good)))
    named: list[tuple[str, Any]] = [
        ("crs null", None), ("crs {}", {}), ("crs name null", {"name": None}),
        ("crs without name", {}), ("crs unreadable name", {"name": "EPSG:1"}),
    ]  # fmt: skip
    for name, props in named:
        crs = props if name in ("crs null", "crs {}") else {"type": "name", "properties": props}
        out.append((name, {**collection(good), "crs": crs}))
    out.append(("no crs member", {"type": "FeatureCollection", "features": [good]}))
    out.append(("collection without type", {"crs": CRS_MEMBER, "features": [good]}))
    out.append(('"features": {}', {**collection(), "features": {}}))
    out.append(('"features": [7]', {**collection(), "features": [7]}))
    out.append(("a JSON list", [good]))
    out.append(("not JSON", b"{not json"))
    out.append(("byte order mark", b"\xef\xbb\xbf" + json.dumps(collection(good)).encode()))
    return out


def readers() -> list[tuple[str, str, str, Callable[[Path], str], tuple[type, ...]]]:
    """(name, right type, wrong type, read and summarise, what the CLI catches)."""
    from shapely.geometry import box

    from tin_engine import feature_input as fi
    from tin_engine.domain import DomainError, DomainPolygon
    from tin_engine.io import rivers, station_set

    try:
        from tin_engine.io.domain_file import read_domain
    except ImportError:  # before PR C
        from tin_engine.domain import read_domain  # type: ignore[attr-defined,no-redef]
    nve = getattr(station_set, "read_nve_lakes", None) or station_set.read_lakes
    lakes = getattr(fi, "read_lake_polygons", None) or fi.read_lakes
    cmap = fi.CLASS_MAPS["property"]
    region = DomainPolygon(polygon=box(X - 1000, Y - 1000, X + 1000, Y + 1000), crs="EPSG:25833")
    caught = (OSError, ValueError)

    def domain(p: Path) -> str:
        return f"a {read_domain(p).polygon.geom_type}"

    def source(p: Path) -> str:
        rows = fi.read_source(p, None, cmap.attribute, lambda own: (0.0, 0.0, 0.0, 0.0)).rows
        # read_source's rows as it returns them: a geometry by its type, anything else as is
        return str([(fid, getattr(g, "geom_type", g)) for fid, g, _ in rows])

    def features(p: Path) -> str:
        request = fi.FeatureRequest(sources=(fi.FeatureSource(path=p, class_map=cmap),))
        got = fi.open_features(request, region, "EPSG:25833")
        return f"{len(got.features)} features, {got.empty} empty, {got.outside} outside"

    def catchment(p: Path) -> str:
        return f"{len(lakes(p, None, (X, Y), 'EPSG:25833')[0])} lakes"

    def stations(p: Path) -> str:
        return str([s.station for s in station_set.read_stations(p)[0]])

    def references(p: Path) -> str:
        return str(sorted(station_set.read_references(p)[0]))

    def nve_lakes(p: Path) -> str:
        return f"{len(nve(p)[0])} lakes"

    def segments(p: Path) -> str:
        return str([s.objectid for s in rivers.read_segments(p)[0]])

    return [
        ("--domain", "Polygon", "Point", domain, (DomainError,)),
        ("read_source", "LineString", "Point", source, (fi.FeatureError,)),
        ("--features", "LineString", "Point", features, (fi.FeatureError,)),
        ("catchment --lakes", "Polygon", "Point", catchment, (fi.FeatureError,)),
        ("stations", "Point", "LineString", stations, caught),
        ("references", "Polygon", "Point", references, caught),
        ("NVE lakes", "Polygon", "Point", nve_lakes, caught),
        ("rivers", "LineString", "Point", segments, caught),
    ]


def main(root: Path) -> int:
    import tin_engine

    where = Path(tin_engine.__file__).resolve()
    assert where.is_relative_to((root / "src_python").resolve()), f"tin_engine from {where}"
    from tin_engine.features import DEFAULT_VOCABULARY

    prop = DEFAULT_VOCABULARY.properties[0].name
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "probe.geojson"
        for name, right, wrong, read, caught in readers():
            for shape, doc in shapes(right, wrong, prop):
                path.write_bytes(doc if isinstance(doc, bytes) else json.dumps(doc).encode())
                try:
                    outcome = f"read: {read(path)}"
                except Exception as exc:  # the probe's subject: every exception
                    text = " ".join(str(exc).replace(tmp, "<dir>").split())
                    kind = f"UNCAUGHT {type(exc).__name__}"
                    outcome = f"{'refused' if isinstance(exc, caught) else kind}: {text}"
                print(f"{name} | {shape} | {outcome}")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
