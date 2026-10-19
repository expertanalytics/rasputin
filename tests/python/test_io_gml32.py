"""`tin_engine.io.gml.iter_gml32`: GML 3.2 surfaces, streamed one feature at a time.

`docs/increments/73-gml-vector-source.md`, sections 4.1 and 9, tests 1 to 4:
fixture A (six features of NIBIO's AR50 map) decodes to the measured figures;
``gml:Surface`` and ``gml:MultiSurface`` map to shapely; every refusal of
section 4.1 is a `GmlError` naming the feature's ``gml:id``; and memory holds
one feature at a time.

The decoder is imported inside each test, so each fails on its own while
`iter_gml32` does not exist.
"""

from __future__ import annotations

import importlib
import io
import tracemalloc
from collections.abc import Iterator
from types import ModuleType
from typing import Any

import pytest
from shapely.geometry import MultiPolygon, Polygon

from gml32_fixtures import (
    AR50_FIGURES,
    AR50_SAMPLE,
    AR50_VALUES,
    GML32_NS,
    HOLE,
    SQUARE,
    SQUARE_2,
    URI_25832,
    document,
    feature,
    polygon,
    ring,
    rings,
    surface,
    values,
)


@pytest.fixture(scope="module")
def gml() -> ModuleType:
    return importlib.import_module("tin_engine.io.gml")


def decoded(gml: ModuleType, stream: io.BytesIO) -> list[Any]:
    return list(gml.iter_gml32(stream))


def one(gml: ModuleType, *geometries: str, props: str | None = None) -> Any:
    (found,) = decoded(gml, document(feature("f1", *geometries, props=props)))
    return found


def refused(gml: ModuleType, stream: io.BytesIO, *words: str) -> str:
    """The `GmlError` text of decoding ``stream``, which must contain each of ``words``."""
    with pytest.raises(gml.GmlError) as caught:
        decoded(gml, stream)
    message = str(caught.value)
    for word in words:
        assert word in message, f"{word!r} not in {message!r}"
    return message


# ---------------------------------------------------------------- test 1


@pytest.fixture(scope="module")
def features(gml: ModuleType) -> list[Any]:
    """Fixture A, decoded once."""
    with AR50_SAMPLE.open("rb") as stream:
        return list(gml.iter_gml32(stream))


class TestFixtureA:
    def test_six_features_in_document_order(self, features: list[Any]) -> None:
        assert [f.gml_id for f in features] == [f.gml_id for f in AR50_FIGURES]

    @pytest.mark.parametrize("index", range(6), ids=[f.gml_id[:10] for f in AR50_FIGURES])
    def test_each_feature_has_the_measured_rings_and_area(
        self, features: list[Any], index: int
    ) -> None:
        want, got = AR50_FIGURES[index], features[index]
        shape = got.geometry
        assert isinstance(shape, Polygon) and shape.is_valid
        assert len(shape.interiors) == want.holes
        assert len(shape.exterior.coords) == want.exterior_positions
        assert sum(len(r.coords) for r in shape.interiors) == want.interior_positions
        # Areas are pinned to 0.1 m2 as measured; the largest is 353,418.5 m2.
        assert round(shape.area, 1) == want.area_m2
        assert got.values["arealtype"] == want.arealtype
        assert got.values["treslag"] == want.treslag
        assert got.values["skogbonitet"] == want.skogbonitet

    def test_type_and_srs_name(self, features: list[Any]) -> None:
        assert {f.type_name for f in features} == {"ArealressursFlate"}
        assert {f.srs_name for f in features} == {URI_25832}

    def test_values_are_the_leaf_properties_in_order(self, features: list[Any]) -> None:
        for f in features:
            assert tuple(f.values) == AR50_VALUES, f.gml_id

    def test_values_are_stripped_text(self, features: list[Any]) -> None:
        first = features[0].values
        assert (
            first["informasjon"]
            == "AR5 Årsversjon 2024. ARFJELL2. N50 per 15.01.2025. Produsert 13.3.2025"
        )
        assert first["oppdateringsdato"] == "2025-04-16T11:01:33.554Z"
        assert first["kartstandard"] == "AR50"

    def test_the_feature_is_frozen(self, features: list[Any]) -> None:
        with pytest.raises(AttributeError):
            features[0].gml_id = "other"


# ---------------------------------------------------------------- test 2


class TestSurfacesMap:
    def test_a_surface_of_two_patches_is_a_multipolygon(self, gml: ModuleType) -> None:
        found = one(gml, surface(SQUARE, SQUARE_2))
        assert isinstance(found.geometry, MultiPolygon)
        assert len(found.geometry.geoms) == 2
        assert found.srs_name == URI_25832

    def test_a_surface_of_one_patch_is_a_polygon(self, gml: ModuleType) -> None:
        found = one(gml, surface(SQUARE))
        assert isinstance(found.geometry, Polygon)
        assert found.geometry.area == 100.0

    def test_a_multisurface_keeps_its_members_and_their_holes(self, gml: ModuleType) -> None:
        members = "".join(
            f"<gml:surfaceMember>{p}</gml:surfaceMember>"
            for p in (polygon(SQUARE, HOLE, srs_name=None), polygon(SQUARE_2, srs_name=None))
        )
        found = one(gml, f'<gml:MultiSurface srsName="{URI_25832}">{members}</gml:MultiSurface>')
        assert isinstance(found.geometry, MultiPolygon)
        assert len(found.geometry.geoms) == 2
        assert sum(len(p.interiors) for p in found.geometry.geoms) == 1
        assert found.geometry.area == 196.0
        assert found.srs_name == URI_25832  # the members inherit it

    def test_surface_members_and_a_members_patches_are_flattened(self, gml: ModuleType) -> None:
        many = polygon(SQUARE, srs_name=None) + surface(
            "40 0 50 0 50 10 40 10 40 0", SQUARE_2, srs_name=None
        )
        found = one(
            gml,
            f'<gml:MultiSurface srsName="{URI_25832}">'
            f"<gml:surfaceMembers>{many}</gml:surfaceMembers></gml:MultiSurface>",
        )
        assert isinstance(found.geometry, MultiPolygon)
        assert len(found.geometry.geoms) == 3

    def test_holes_are_kept_in_order(self, gml: ModuleType) -> None:
        second = "6 6 6 8 8 8 8 6 6 6"
        found = one(gml, polygon(SQUARE, HOLE, second))
        assert [r.coords[0] for r in found.geometry.interiors] == [(2.0, 2.0), (6.0, 6.0)]

    def test_a_three_dimensional_pos_list_drops_the_height(self, gml: ModuleType) -> None:
        three = "0 0 5 10 0 6 10 10 7 0 10 8 0 0 5"
        found = one(gml, f"<gml:Polygon>{rings(ring(three, dimension=3))}</gml:Polygon>")
        assert list(found.geometry.exterior.coords) == [
            (0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)
        ]  # fmt: skip
        assert not found.geometry.has_z

    def test_the_geometrys_srs_dimension_applies_when_the_pos_list_has_none(
        self, gml: ModuleType
    ) -> None:
        three = "0 0 5 10 0 6 10 10 7 0 10 8 0 0 5"
        found = one(gml, f'<gml:Polygon srsDimension="3">{rings(ring(three))}</gml:Polygon>')
        assert found.geometry.area == 100.0
        assert not found.geometry.has_z

    def test_coordinates_are_read_bit_for_bit(self, gml: ModuleType) -> None:
        exact = "620989.187637389 6769027.11937103 620997.111133458 6769022.8280574 621012.958184295 6769014.24548286 620989.187637389 6769027.11937103"  # noqa: E501
        found = one(gml, polygon(exact))
        assert found.geometry.exterior.coords[1] == (620997.111133458, 6769022.8280574)

    def test_rings_are_stored_as_served(self, gml: ModuleType) -> None:
        bowtie = "0 0 10 10 10 0 0 10 0 0"  # self-crossing: not repaired or refused
        found = one(gml, polygon(bowtie))
        assert not found.geometry.is_valid
        assert len(found.geometry.exterior.coords) == 5

    def test_a_feature_with_no_geometry_is_yielded_with_none(self, gml: ModuleType) -> None:
        found = one(gml)
        assert found.geometry is None
        assert found.srs_name is None
        assert dict(found.values) == {"a": "60"}

    def test_an_empty_or_nil_value_is_none(self, gml: ModuleType) -> None:
        found = one(gml, polygon(), props=values(a="  60 ", b="", c=None))
        assert dict(found.values) == {"a": "60", "b": None, "c": None}

    def test_no_srs_name_is_none(self, gml: ModuleType) -> None:
        assert one(gml, polygon(srs_name=None)).srs_name is None

    def test_several_features_keep_their_order_and_types(self, gml: ModuleType) -> None:
        stream = document(
            feature("f1", polygon()), feature("f2", type_name="G"), feature("f3", polygon())
        )
        found = decoded(gml, stream)
        assert [(f.gml_id, f.type_name) for f in found] == [("f1", "F"), ("f2", "G"), ("f3", "F")]


# ---------------------------------------------------------------- test 3


def ring_with(inner: str) -> str:
    return f'<gml:Polygon srsName="{URI_25832}"><gml:exterior>{inner}</gml:exterior></gml:Polygon>'


POS_RING = "<gml:LinearRing>" + "".join(
    f"<gml:pos>{p}</gml:pos>" for p in ("0 0", "10 0", "10 10", "0 0")
) + "</gml:LinearRing>"  # fmt: skip
COORDINATES_RING = (
    "<gml:LinearRing><gml:coordinates>0,0 10,0 10,10 0,0</gml:coordinates></gml:LinearRing>"
)
CURVE_RING = (
    "<gml:Ring><gml:curveMember><gml:LineString><gml:posList>0 0 10 0 10 10 0 0"
    "</gml:posList></gml:LineString></gml:curveMember></gml:Ring>"
)


TRIANGLE = "0 0 10 0 10 10 0 0"


def dimensioned(on: str, dimension: str) -> str:
    """The square as a ``gml:Polygon`` whose ``srsDimension`` is ``dimension``,
    given ``on`` the geometry element or its ``gml:posList``; the positions
    have that many numbers each when it is 3, else 2."""
    positions = "0 0 5 10 0 6 10 10 7 0 10 8 0 0 5" if dimension == "3" else SQUARE
    attribute = f' srsDimension="{dimension}"'
    if on == "geometry":
        return f"<gml:Polygon{attribute}>{rings(ring(positions))}</gml:Polygon>"
    pos_list = f"<gml:posList{attribute}>{positions}</gml:posList>"
    return f"<gml:Polygon>{rings(f'<gml:LinearRing>{pos_list}</gml:LinearRing>')}</gml:Polygon>"


#: A 10 by 5 rectangle with a vertex mid-edge, in three dimensions: 18
#: numbers, which read as two dimensions are nine other positions, so a
#: reader that misses the dimension decodes a wrong ring instead of failing.
RECTANGLE_3D = "0 0 1 5 0 2 10 0 3 10 5 4 0 5 5 0 0 1"
RECTANGLE_XY = [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0), (0.0, 0.0)]
#: The same rectangle moved by (2, 2), as a hole of the 20 by 20 square.
HOLE_3D = "2 2 1 7 2 2 12 2 3 12 7 4 2 7 5 2 2 1"
HOLE_XY = [(x + 2.0, y + 2.0) for x, y in RECTANGLE_XY]
OUTER_20 = "0 0 20 0 20 20 0 20 0 0"


def dimension_inside(on: str, dimension: str) -> str:
    """A surface whose ``srsDimension`` is ``dimension`` on the element ``on``
    names, every ``gml:posList`` without one. Its ring there holds
    `RECTANGLE_3D` (or `HOLE_3D`) when ``dimension`` is 3, else the 2-D
    rectangle (or hole); the other rings are two-dimensional."""
    three = dimension == "3"
    rectangle = RECTANGLE_3D if three else " ".join(f"{x:g} {y:g}" for x, y in RECTANGLE_XY)
    hole = HOLE_3D if three else " ".join(f"{x:g} {y:g}" for x, y in HOLE_XY)
    attribute = f' srsDimension="{dimension}"'
    bare = f"<gml:LinearRing><gml:posList>{rectangle}</gml:posList></gml:LinearRing>"
    marked = f"<gml:LinearRing{attribute}><gml:posList>{rectangle}</gml:posList></gml:LinearRing>"
    if on == "exterior LinearRing":
        return f"<gml:Polygon>{rings(marked)}</gml:Polygon>"
    if on == "interior LinearRing":
        inner = f"<gml:LinearRing{attribute}><gml:posList>{hole}</gml:posList></gml:LinearRing>"
        return f"<gml:Polygon>{rings(ring(OUTER_20), inner)}</gml:Polygon>"
    if on == "MultiSurface":
        member = f"<gml:Polygon>{rings(bare)}</gml:Polygon>"
        return f"<gml:MultiSurface{attribute}><gml:surfaceMember>{member}</gml:surfaceMember></gml:MultiSurface>"  # noqa: E501
    if on == "Polygon in a MultiSurface":
        member = f"<gml:Polygon{attribute}>{rings(bare)}</gml:Polygon>"
    else:  # "Surface in a MultiSurface"
        patch = f"<gml:PolygonPatch>{rings(bare)}</gml:PolygonPatch>"
        member = f"<gml:Surface{attribute}><gml:patches>{patch}</gml:patches></gml:Surface>"
    return f"<gml:MultiSurface><gml:surfaceMember>{member}</gml:surfaceMember></gml:MultiSurface>"


#: Every element of what the reader accepts that may carry ``srsDimension``
#: (GML 3.2's SRSReferenceGroup: a geometry, or a ``gml:posList``), other than
#: the outermost surface and the posList, which `dimensioned` covers.
INNER_DIMENSION_ON = [
    "MultiSurface",
    "Polygon in a MultiSurface",
    "Surface in a MultiSurface",
    "exterior LinearRing",
    "interior LinearRing",
]


class TestRefusals:
    @pytest.mark.parametrize(
        ("geometry", "name"),
        [
            pytest.param(
                f'<gml:LineString srsName="{URI_25832}">'
                "<gml:posList>0 0 1 1</gml:posList></gml:LineString>",
                "gml:LineString", id="line",
            ),
            pytest.param(
                f'<gml:Point srsName="{URI_25832}"><gml:pos>0 0</gml:pos></gml:Point>',
                "gml:Point", id="point",
            ),
            pytest.param(
                f'<gml:MultiCurve srsName="{URI_25832}"><gml:curveMember><gml:LineString>'
                "<gml:posList>0 0 1 1</gml:posList></gml:LineString></gml:curveMember>"
                "</gml:MultiCurve>",
                "gml:MultiCurve", id="multicurve",
            ),
        ],
    )  # fmt: skip
    def test_a_curve_or_point_geometry(self, gml: ModuleType, geometry: str, name: str) -> None:
        refused(gml, document(feature("f.curve", geometry)), "f.curve", "surfaces only", name)

    @pytest.mark.parametrize(
        ("geometry", "name"),
        [
            pytest.param(
                f'<gml:CompositeSurface srsName="{URI_25832}"><gml:surfaceMember>'
                f"{polygon(srs_name=None)}</gml:surfaceMember></gml:CompositeSurface>",
                "gml:CompositeSurface", id="composite",
            ),
            pytest.param(
                f'<gml:OrientableSurface srsName="{URI_25832}" orientation="-">'
                f"<gml:baseSurface>{polygon(srs_name=None)}</gml:baseSurface>"
                "</gml:OrientableSurface>",
                "gml:OrientableSurface", id="orientable",
            ),
            pytest.param(
                f'<gml:PolyhedralSurface srsName="{URI_25832}"><gml:patches>'
                f"<gml:PolygonPatch>{rings(ring(SQUARE))}</gml:PolygonPatch>"
                "</gml:patches></gml:PolyhedralSurface>",
                "gml:PolyhedralSurface", id="polyhedral",
            ),
            pytest.param(
                f'<gml:TriangulatedSurface srsName="{URI_25832}"><gml:patches>'
                f"<gml:Triangle>{rings(ring(TRIANGLE))}</gml:Triangle>"
                "</gml:patches></gml:TriangulatedSurface>",
                "gml:TriangulatedSurface", id="triangulated",
            ),
            pytest.param(
                f'<gml:Tin srsName="{URI_25832}"><gml:patches>'
                f"<gml:Triangle>{rings(ring(TRIANGLE))}</gml:Triangle></gml:patches>"
                "<gml:controlPoint><gml:posList>0 0 10 0 10 10</gml:posList>"
                "</gml:controlPoint></gml:Tin>",
                "gml:Tin", id="tin",
            ),
            pytest.param(
                f'<gml:MultiSurface srsName="{URI_25832}"><gml:surfaceMember>'
                f"<gml:CompositeSurface><gml:surfaceMember>{polygon(srs_name=None)}"
                "</gml:surfaceMember></gml:CompositeSurface></gml:surfaceMember>"
                "</gml:MultiSurface>",
                "gml:CompositeSurface", id="composite-in-multisurface",
            ),
        ],
    )  # fmt: skip
    def test_a_surface_kind_not_read_is_refused_not_dropped(
        self, gml: ModuleType, geometry: str, name: str
    ) -> None:
        """A surface kind the reader does not map is refused in words, at the
        top or inside a MultiSurface, never yielded with no geometry."""
        refused(gml, document(feature("f.kind", geometry)), "f.kind", name)

    @pytest.mark.parametrize("bad", ["0", "1", "x", "-2", "2.5"])
    @pytest.mark.parametrize("on", ["geometry", "posList"])
    def test_an_srs_dimension_not_a_whole_number_of_two_or_more(
        self, gml: ModuleType, bad: str, on: str
    ) -> None:
        """A bad ``srsDimension`` on the geometry or its ``gml:posList`` is
        refused naming the feature; ``0`` would otherwise divide by zero."""
        refused(gml, document(feature("f.dim", dimensioned(on, bad))), "f.dim", "srsDimension")

    @pytest.mark.parametrize("good", ["2", "3"])
    @pytest.mark.parametrize("on", ["geometry", "posList"])
    def test_an_srs_dimension_of_two_or_three_is_read(
        self, gml: ModuleType, good: str, on: str
    ) -> None:
        found = one(gml, dimensioned(on, good))
        assert found.geometry.area == 100.0
        assert not found.geometry.has_z

    @pytest.mark.parametrize("on", INNER_DIMENSION_ON)
    def test_an_srs_dimension_of_three_inside_the_surface_is_read(
        self, gml: ModuleType, on: str
    ) -> None:
        """``srsDimension="3"`` on any element between the outermost surface
        and a bare ``gml:posList`` makes its positions three numbers each; the
        height is dropped and the (x, y) points are the rectangle's."""
        found = one(gml, dimension_inside(on, "3"))
        (shape,) = getattr(found.geometry, "geoms", [found.geometry])
        if on == "interior LinearRing":
            assert shape.area == 350.0
            assert [list(r.coords) for r in shape.interiors] == [HOLE_XY]
        else:
            assert shape.area == 50.0
            assert list(shape.exterior.coords) == RECTANGLE_XY
        assert not shape.has_z

    @pytest.mark.parametrize("bad", ["0", "1", "x"])
    @pytest.mark.parametrize("on", INNER_DIMENSION_ON)
    def test_an_srs_dimension_inside_the_surface_not_a_whole_number_of_two_or_more(
        self, gml: ModuleType, bad: str, on: str
    ) -> None:
        """A bad ``srsDimension`` is refused wherever in the surface it
        stands, naming the feature; ``0`` would otherwise divide by zero."""
        geometry = dimension_inside(on, bad)
        refused(gml, document(feature("f.dim", geometry)), "f.dim", "srsDimension")

    def test_two_surface_properties(self, gml: ModuleType) -> None:
        refused(gml, document(feature("f.two", polygon(), polygon(SQUARE_2))), "f.two")

    @pytest.mark.parametrize(
        ("inner", "name"),
        [
            pytest.param(POS_RING, "gml:pos", id="pos"),
            pytest.param(COORDINATES_RING, "gml:coordinates", id="coordinates"),
            pytest.param(CURVE_RING, "gml:Ring", id="ring"),
        ],
    )
    def test_a_ring_not_of_one_pos_list(self, gml: ModuleType, inner: str, name: str) -> None:
        refused(gml, document(feature("f.ring", ring_with(inner))), "f.ring", name)

    def test_a_count_not_a_multiple_of_the_dimension(self, gml: ModuleType) -> None:
        refused(gml, document(feature("f.odd", polygon("0 0 10 0 10 10 0 10 0"))), "f.odd")

    def test_a_count_not_a_multiple_of_three_dimensions(self, gml: ModuleType) -> None:
        eight = "0 0 5 10 0 6 10 10 7 0 10 8 0 0 5 1"
        geometry = f"<gml:Polygon>{rings(ring(eight, dimension=3))}</gml:Polygon>"
        refused(gml, document(feature("f.odd3", geometry)), "f.odd3")

    def test_a_ring_of_three_positions(self, gml: ModuleType) -> None:
        refused(gml, document(feature("f.short", polygon("0 0 10 0 0 0"))), "f.short")

    def test_a_short_hole(self, gml: ModuleType) -> None:
        refused(gml, document(feature("f.hole", polygon(SQUARE, "2 2 4 4 2 2"))), "f.hole")

    @pytest.mark.parametrize("bad", ["NaN", "nan", "inf", "-INF", "Infinity"])
    def test_a_coordinate_that_is_not_finite(self, gml: ModuleType, bad: str) -> None:
        positions = f"0 0 10 0 10 {bad} 0 10 0 0"
        refused(gml, document(feature("f.nan", polygon(positions))), "f.nan")

    def test_a_coordinate_that_is_not_a_number(self, gml: ModuleType) -> None:
        refused(gml, document(feature("f.text", polygon("0 0 10 0 10 x 0 10 0 0"))), "f.text")

    def test_malformed_xml_names_the_parsers_line(self, gml: ModuleType) -> None:
        text = document(feature("f1", polygon())).getvalue().replace(b"</app:F>", b"</app:G>")
        line = text[: text.index(b"</app:G>")].count(b"\n") + 1
        refused(gml, io.BytesIO(text), f"line {line}")

    def test_a_document_cut_short(self, gml: ModuleType) -> None:
        text = document(feature("f1", polygon()), feature("f2", polygon())).getvalue()
        with pytest.raises(gml.GmlError):
            decoded(gml, io.BytesIO(text[: len(text) // 2]))

    def test_a_gml2_document(self, gml: ModuleType) -> None:
        gml2 = (
            b'<?xml version="1.0"?><ogr:FeatureCollection xmlns:ogr="http://ogr.maptools.org/"'
            b' xmlns:gml="http://www.opengis.net/gml"><gml:featureMember><ogr:f fid="1">'
            b"<ogr:a>1</ogr:a></ogr:f></gml:featureMember></ogr:FeatureCollection>"
        )
        refused(gml, io.BytesIO(gml2), "no GML 3.2 gml:featureMember")

    def test_a_wfs_member_document(self, gml: ModuleType) -> None:
        wfs = (
            '<?xml version="1.0"?><wfs:FeatureCollection'
            ' xmlns:wfs="http://www.opengis.net/wfs/2.0"'
            f' xmlns:gml="{GML32_NS}" xmlns:app="urn:test">'
            f"<wfs:member>{feature('f1', polygon())}</wfs:member></wfs:FeatureCollection>"
        )
        refused(gml, io.BytesIO(wfs.encode()), "no GML 3.2 gml:featureMember")

    def test_an_empty_collection(self, gml: ModuleType) -> None:
        refused(gml, document(), "no GML 3.2 gml:featureMember")

    def test_the_error_is_the_readers_value_error(self, gml: ModuleType) -> None:
        assert issubclass(gml.GmlError, ValueError)


# ---------------------------------------------------------------- test 4

#: The generator of `docs/increments/73-probes/mem_probe.py`.
MEM_HEAD = (
    b'<?xml version="1.0"?><gml:FeatureCollection'
    b' xmlns:gml="http://www.opengis.net/gml/3.2" xmlns:app="urn:x">'
)
MEM_FEATURE = (
    b'<gml:featureMember><app:F gml:id="f%d"><app:a>60</app:a><app:g>'
    b'<gml:Polygon srsName="EPSG:25832">'
    b"<gml:exterior><gml:LinearRing><gml:posList>0 0 10 0 10 10 0 10 0 0</gml:posList>"
    b"</gml:LinearRing></gml:exterior>"
    b"</gml:Polygon></app:g></app:F></gml:featureMember>"
)
#: A chosen bound, about 6 times the measured peak with the element cleared
#: (167 KiB at 20,000 features); a reader that keeps the tree peaks near
#: 30,000 KiB there.
PEAK_BOUND = 1 << 20


def synthetic(n: int) -> io.BytesIO:
    body = b"".join(MEM_FEATURE % k for k in range(n))
    return io.BytesIO(MEM_HEAD + body + b"</gml:FeatureCollection>")


def peak_while_consuming(features: Iterator[Any]) -> tuple[int, int]:
    """Features consumed one at a time, none kept, and the traced peak in bytes."""
    tracemalloc.start()
    try:
        count = sum(1 for _ in features)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    return count, peak


class TestMemory:
    def test_the_document_is_the_measured_one(self) -> None:
        assert len(synthetic(20_000).getvalue()) == 5_229_019

    def test_memory_holds_one_feature_at_a_time(self, gml: ModuleType) -> None:
        stream = synthetic(20_000)
        count, peak = peak_while_consuming(gml.iter_gml32(stream))
        assert count == 20_000
        assert peak < PEAK_BOUND, f"peak {peak // 1024} KiB"
