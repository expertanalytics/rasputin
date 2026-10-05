"""`tin_engine.io.geojson.catchment_geojson`: the catchment file's bytes (increment 29, PR 4).

`docs/increments/29-nve-reference-catchments.md`, "The batch" (the last
paragraph: "The GeoJSON writer moves"), and `project_structure.md` ("The
catchment file's bytes come from `io/geojson.py`'s `catchment_geojson`, which
opens nothing"). In PR 4 `station-catchments`
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

Audit PR C (`docs/increments/python-audit.md`, section 10, red test 1): the
module also holds the one GeoJSON reading path, `read_collection(doc, *,
default_crs) -> (features, crs_text)`, which takes the parsed document (it
opens nothing; the caller adds the file's name to a refusal), and
`feature_collection(crs, features)`, the document both writers build. The
`crs` rule and the shape rule are the design's "The one reading path"; each
refusal below is its wording, whole.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any

import pytest
from shapely.geometry import Polygon

import tin_engine.io.geojson as gj
from tin_engine.crs import parse_crs
from tin_engine.io.domain_file import read_domain

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src_python" / "tin_engine"
#: The bench's domain: a bare Polygon whose `crs` member names this URN.
QUARTER = ROOT / "docs" / "benchmarks" / "2026-09-26" / "quarter.geojson"
EPSG = "EPSG:25833"
#: Coordinates that need every digit of their repr to survive.
RING = [
    (500000.1, 6600000.3),
    (500123.45678901234, 6600000.3),
    (500123.45678901234, 6600210.000000001),
    (500000.1, 6600210.000000001),
]


def test_it_takes_no_path() -> None:
    params = list(inspect.signature(gj.catchment_geojson).parameters)
    assert params == ["polygon", "crs", "properties"]


def test_the_bytes_are_one_feature_with_a_crs_member() -> None:
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


def test_non_ascii_names_survive() -> None:
    out = gj.catchment_geojson(Polygon(RING), EPSG, {"name": "Atnasjø"})
    assert json.loads(out.decode("utf-8"))["features"][0]["properties"]["name"] == "Atnasjø"


def test_read_domain_reads_the_bytes_back_exactly(tmp_path: Path) -> None:
    path = tmp_path / "c.geojson"
    path.write_bytes(gj.catchment_geojson(Polygon(RING), EPSG, {}))
    domain = read_domain(path)
    assert parse_crs(domain.crs) == parse_crs(EPSG)
    # `read_domain` may turn the ring round; every vertex survives bit for bit.
    assert set(domain.polygon.exterior.coords) == set(RING)


def test_the_writer_has_moved_out_of_cli() -> None:
    """`cli.py` no longer builds the document; `io/geojson.py` builds it and
    opens nothing."""
    cli = (SRC / "cli.py").read_text(encoding="utf-8")
    assert "FeatureCollection" not in cli
    assert "catchment_geojson" in cli
    own = (SRC / "io" / "geojson.py").read_text(encoding="utf-8")
    for opener in ("open(", "write_text", "write_bytes", "read_text", "Path("):
        assert opener not in own, opener


# ------------------------------------------------- audit PR C: read_collection

NOT_AN_OBJECT = (
    "not a GeoJSON object; the file must hold a FeatureCollection, a Feature or a geometry"
)
NO_MEMBER = "no crs member; the file must name its CRS"
NULL_MEMBER = "the crs member is null; the file must name its CRS"
NO_NAME = "the crs member has no name; it must name the CRS"
NO_FEATURES = "no features list; the file is not a FeatureCollection"
WGS84 = "EPSG:4326"
URN = "urn:ogc:def:crs:EPSG::25833"
POINT = {"type": "Point", "coordinates": [500000.5, 6600000.5]}
FEATURE = {"type": "Feature", "properties": {"station": "2.11.0"}, "geometry": POINT}
NAMED = {"type": "name", "properties": {"name": URN}}


def collection(crs: Any = NAMED, **members: Any) -> dict[str, Any]:
    """A FeatureCollection holding `FEATURE`; `crs=...` (Ellipsis) leaves the
    member out, anything else is written as the member."""
    doc: dict[str, Any] = {"type": "FeatureCollection", "features": [FEATURE]} | members
    if crs is not ...:
        doc["crs"] = crs
    return doc


def refusal(doc: object, default_crs: str | None) -> str:
    with pytest.raises(ValueError) as info:
        gj.read_collection(doc, default_crs=default_crs)
    return str(info.value)


class TestReadCollectionTakesTheDocument:
    def test_the_signature(self) -> None:
        """The parsed document, and the default by keyword with no default of
        its own: every caller says whether it has one."""
        params = inspect.signature(gj.read_collection).parameters
        assert list(params) == ["doc", "default_crs"]
        assert params["default_crs"].kind is inspect.Parameter.KEYWORD_ONLY
        assert params["default_crs"].default is inspect.Parameter.empty

    def test_a_path_is_not_a_document(self, tmp_path: Path) -> None:
        """It opens nothing: a path to a good file is refused as no object."""
        path = tmp_path / "c.geojson"
        path.write_text(json.dumps(collection()), encoding="utf-8")
        assert refusal(path, WGS84) == NOT_AN_OBJECT

    def test_one_copy_of_the_suffixes_and_the_default(self) -> None:
        """Moved here from `domain.py`, which keeps none of the three constants."""
        import tin_engine.domain as domain

        assert gj.GEOJSON_SUFFIXES == (".geojson", ".json")
        assert gj.RFC7946_CRS == WGS84
        for gone in ("GEOJSON_SUFFIXES", "GEOJSON_DEFAULT_CRS", "WKT_SUFFIXES"):
            assert not hasattr(domain, gone), gone


DEFAULTS = pytest.mark.parametrize("default_crs", [WGS84, None], ids=["default", "no_default"])


class TestTheCrsRule:
    """Rules 1-5 of "The one reading path", each with and without a default."""

    @DEFAULTS
    @pytest.mark.parametrize(
        "doc", [[collection()], "FeatureCollection", 7, None], ids=["list", "str", "int", "null"]
    )
    def test_rule_1_not_an_object(self, doc: object, default_crs: str | None) -> None:
        assert refusal(doc, default_crs) == NOT_AN_OBJECT

    def test_rule_2_no_member_is_the_default(self) -> None:
        features, crs = gj.read_collection(collection(crs=...), default_crs=WGS84)
        assert (features, crs) == ([FEATURE], WGS84)

    def test_rule_2_no_member_and_no_default_is_refused(self) -> None:
        assert refusal(collection(crs=...), None) == NO_MEMBER

    @DEFAULTS
    def test_rule_3_a_null_member_is_refused_whatever_the_default(
        self, default_crs: str | None
    ) -> None:
        """GeoJSON 2008: "If the value of CRS is null, no CRS can be assumed"."""
        assert refusal(collection(crs=None), default_crs) == NULL_MEMBER

    @DEFAULTS
    @pytest.mark.parametrize(
        "member",
        [
            {},
            "EPSG:25833",
            {"type": "name"},
            {"type": "name", "properties": None},
            {"type": "name", "properties": "EPSG:25833"},
            {"type": "name", "properties": {}},
            {"type": "name", "properties": {"name": None}},
            {"type": "link", "properties": {"href": "http://example.com/crs", "type": "proj4"}},
            [NAMED],
        ],
        ids=[
            "empty_object",
            "a_string",
            "no_properties",
            "null_properties",
            "string_properties",
            "empty_properties",
            "null_name",
            "link",
            "a_list",
        ],
    )
    def test_rule_4_a_member_with_no_name_is_refused(
        self, member: object, default_crs: str | None
    ) -> None:
        """A linked CRS among them: nothing here dereferences a URL."""
        assert refusal(collection(crs=member), default_crs) == NO_NAME

    @DEFAULTS
    def test_rule_4_the_name_is_the_text(self, default_crs: str | None) -> None:
        assert gj.read_collection(collection(), default_crs=default_crs)[1] == URN

    @DEFAULTS
    def test_rule_4_a_number_is_read_as_its_digits(self, default_crs: str | None) -> None:
        member = {"type": "name", "properties": {"name": 4326}}
        crs = gj.read_collection(collection(crs=member), default_crs=default_crs)[1]
        assert crs == "4326"
        assert isinstance(crs, str)

    @DEFAULTS
    @pytest.mark.parametrize("name", ["EPSG:999999", "not a crs"])
    def test_rule_5_an_unreadable_name_is_parse_crs_s_refusal(
        self, name: str, default_crs: str | None
    ) -> None:
        member = {"type": "name", "properties": {"name": name}}
        assert refusal(collection(crs=member), default_crs).startswith(
            f"cannot read the CRS {name!r}"
        )


class TestTheShapeRule:
    """A Feature is one feature, a collection its features, anything else a
    geometry wrapped as one feature with properties `{}`."""

    @DEFAULTS
    def test_a_collection_gives_its_features(self, default_crs: str | None) -> None:
        other = {"type": "Feature", "properties": {}, "geometry": None}
        doc = collection(features=[FEATURE, other])
        features, crs = gj.read_collection(doc, default_crs=default_crs)
        assert (features, crs) == ([FEATURE, other], URN)

    @DEFAULTS
    def test_an_empty_collection_gives_no_features(self, default_crs: str | None) -> None:
        assert gj.read_collection(collection(features=[]), default_crs=default_crs)[0] == []

    @DEFAULTS
    def test_a_collection_without_type_is_read(self, default_crs: str | None) -> None:
        doc = {"crs": NAMED, "features": [FEATURE]}
        assert gj.read_collection(doc, default_crs=default_crs) == ([FEATURE], URN)

    @DEFAULTS
    def test_a_feature_is_one_feature(self, default_crs: str | None) -> None:
        doc = FEATURE | {"crs": NAMED}
        features, crs = gj.read_collection(doc, default_crs=default_crs)
        assert crs == URN
        (feature,) = features
        assert (feature["type"], feature["properties"], feature["geometry"]) == (
            "Feature",
            FEATURE["properties"],
            POINT,
        )

    @DEFAULTS
    def test_a_bare_geometry_is_wrapped(self, default_crs: str | None) -> None:
        doc = POINT | {"crs": NAMED}
        features, crs = gj.read_collection(doc, default_crs=default_crs)
        assert crs == URN
        (feature,) = features
        assert feature["type"] == "Feature"
        assert feature["properties"] == {}
        assert feature["geometry"] == doc

    @DEFAULTS
    @pytest.mark.parametrize(
        "members",
        [{"features": {}}, {"features": [7]}, {"features": [FEATURE, None]}, {"features": None}],
        ids=["an_object", "a_number", "a_null_feature", "null"],
    )
    def test_features_must_be_a_list_of_objects(
        self, members: dict[str, Any], default_crs: str | None
    ) -> None:
        assert refusal(collection(**members), default_crs) == NO_FEATURES

    @DEFAULTS
    def test_a_feature_collection_without_features_is_refused(
        self, default_crs: str | None
    ) -> None:
        doc = {"type": "FeatureCollection", "crs": NAMED}
        assert refusal(doc, default_crs) == NO_FEATURES

    def test_the_benchs_domain(self) -> None:
        """`quarter.geojson`: a bare Polygon naming UTM 33 by URN, read as its
        member's text and one feature whose geometry is the file's own object."""
        doc = json.loads(QUARTER.read_text(encoding="utf-8"))
        assert doc["type"] == "Polygon"
        features, crs = gj.read_collection(doc, default_crs=gj.RFC7946_CRS)
        assert crs == doc["crs"]["properties"]["name"] == URN
        (feature,) = features
        assert feature["geometry"] == doc
        assert feature["properties"] == {}


class TestFeatureCollection:
    def test_the_document(self) -> None:
        """Key order is part of it: both writers' bytes stay as they were."""
        doc = gj.feature_collection("EPSG:25833", [])
        assert doc == {
            "type": "FeatureCollection",
            "crs": {"type": "name", "properties": {"name": "EPSG:25833"}},
            "features": [],
        }
        assert list(doc) == ["type", "crs", "features"]

    def test_the_features_are_the_given_list(self) -> None:
        doc = gj.feature_collection(URN, [FEATURE])
        assert doc["features"] == [FEATURE]

    def test_it_reads_back(self) -> None:
        doc = gj.feature_collection(URN, [FEATURE])
        assert gj.read_collection(doc, default_crs=None) == ([FEATURE], URN)
