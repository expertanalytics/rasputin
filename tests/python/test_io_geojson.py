"""`tin_engine.io.geojson.catchment_geojson`: the catchment file's bytes (increment 29, PR 4).

`docs/increments/29-nve-reference-catchments.md`, "The batch" (the last
paragraph: "The GeoJSON writer moves"), and `project_structure.md` ("The
catchment GeoJSON writer is in `cli.py`"). In PR 4 `station-catchments`
becomes the second command to write a catchment file, so the writer moves
from `cli.py` to `io/geojson.py` as `catchment_geojson(polygon, crs,
properties) -> bytes`: no path, nothing opened, the shape `io/ply.py` and
`io/vtk_legacy.py` have. Both commands call it; `cli.py` keeps only the
write.

Pinned here beyond the design: the bytes are UTF-8 JSON of a
`FeatureCollection` with a `crs` member `{"type": "name", "properties":
{"name": crs}}` and one `Feature` whose geometry is the polygon's exterior
ring (as `rasputin catchment` wrote it in increment 22); the moved code
leaves no `FeatureCollection` in `cli.py`, and `io/geojson.py` opens no file.
"""

from __future__ import annotations

import importlib
import inspect
import json
from pathlib import Path
from types import ModuleType

import pytest
from shapely.geometry import Polygon

from tin_engine.crs import parse_crs
from tin_engine.domain import read_domain

SRC = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine"
EPSG = "EPSG:25833"
#: Coordinates that need every digit of their repr to survive.
RING = [
    (500000.1, 6600000.3),
    (500123.45678901234, 6600000.3),
    (500123.45678901234, 6600210.000000001),
    (500000.1, 6600210.000000001),
]


@pytest.fixture(scope="module")
def gj() -> ModuleType:
    return importlib.import_module("tin_engine.io.geojson")


def test_it_takes_no_path(gj: ModuleType) -> None:
    params = list(inspect.signature(gj.catchment_geojson).parameters)
    assert params == ["polygon", "crs", "properties"]


def test_the_bytes_are_one_feature_with_a_crs_member(gj: ModuleType) -> None:
    props = {"station": "2.32.0", "nodes": 12, "causes": ["swing"], "swing": 0.25}
    out = gj.catchment_geojson(Polygon(RING), EPSG, props)
    assert isinstance(out, bytes)
    doc = json.loads(out.decode("utf-8"))
    assert doc["type"] == "FeatureCollection"
    assert doc["crs"] == {"type": "name", "properties": {"name": EPSG}}
    (feature,) = doc["features"]
    assert feature["type"] == "Feature"
    assert feature["properties"] == props
    assert feature["geometry"]["type"] == "Polygon"
    (ring,) = feature["geometry"]["coordinates"]
    assert [tuple(p) for p in ring] == list(Polygon(RING).exterior.coords)


def test_non_ascii_names_survive(gj: ModuleType) -> None:
    out = gj.catchment_geojson(Polygon(RING), EPSG, {"name": "Atnasjø"})
    assert json.loads(out.decode("utf-8"))["features"][0]["properties"]["name"] == "Atnasjø"


def test_read_domain_reads_the_bytes_back_exactly(gj: ModuleType, tmp_path: Path) -> None:
    path = tmp_path / "c.geojson"
    path.write_bytes(gj.catchment_geojson(Polygon(RING), EPSG, {}))
    domain = read_domain(path)
    assert parse_crs(domain.crs) == parse_crs(EPSG)
    # `read_domain` may turn the ring round; every vertex survives bit for bit.
    assert set(domain.polygon.exterior.coords) == set(RING)


def test_the_writer_has_moved_out_of_cli(gj: ModuleType) -> None:
    """`cli.py` no longer builds the document; `io/geojson.py` builds it and
    opens nothing."""
    cli = (SRC / "cli.py").read_text(encoding="utf-8")
    assert "FeatureCollection" not in cli
    assert "catchment_geojson" in cli
    own = (SRC / "io" / "geojson.py").read_text(encoding="utf-8")
    for opener in ("open(", "write_text", "write_bytes", "read_text", "Path("):
        assert opener not in own, opener
