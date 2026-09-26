"""`tin_engine.domain`: reading the domain polygon. Increment 16, R1.

`docs/increments/16-domain-polygon.md`, R1, U1 (a) and U4 (a). The design
names `DomainPolygon` and `read_domain(path, crs)`; the extent check needs the
DEM's node rectangle, so the signature chosen here is::

    read_domain(path: Path, meta: RasterMeta, crs: str | None = None) -> DomainPolygon

`DomainPolygon` is frozen with `.polygon` (a shapely `Polygon`, oriented outer
counter-clockwise and holes clockwise) and `.epsg` (int). Every refusal raises
`DomainError`, a `ValueError` subclass (name chosen here), whose message the
CLI shows. `crs` is the `--domain-crs` text, `EPSG:n`; it is what a `.wkt`
file needs.

The module is fetched inside a fixture, so its absence fails these tests and
leaves the rest of the session collecting.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import shapely
from shapely.geometry import Polygon, mapping

from tin_engine.io.models import RasterMeta

UTM33 = "urn:ogc:def:crs:EPSG::25833"

# Nodes at x 500 000 .. 500 080 and y 6 599 970 .. 6 600 000.
META = RasterMeta(
    x_min=500_000.0,
    y_max=6_600_000.0,
    delta_x=10.0,
    delta_y=5.0,
    cols=9,
    rows=7,
    epsg=25833,
    nodata=None,
    nodata_source="absent",
    pixel_is_area=False,
    vertical_unit_assumed=True,
)

# Off-node, counter-clockwise, and a clockwise hole.
OUTER = [
    (500_012.3, 6_599_972.1),
    (500_071.7, 6_599_973.4),
    (500_066.2, 6_599_996.9),
    (500_004.1, 6_599_991.3),
]
HOLE = [
    (500_030.1, 6_599_980.2),
    (500_030.9, 6_599_988.8),
    (500_050.5, 6_599_987.7),
    (500_049.3, 6_599_981.1),
]


@pytest.fixture
def domain() -> ModuleType:
    return importlib.import_module("tin_engine.domain")


def geojson(
    geometry: dict[str, Any], crs: str | None = UTM33, wrap: str = "geometry"
) -> dict[str, Any]:
    feature = {"type": "Feature", "properties": {}, "geometry": geometry}
    doc: dict[str, Any]
    if wrap == "geometry":
        doc = dict(geometry)
    elif wrap == "feature":
        doc = feature
    else:
        doc = {"type": "FeatureCollection", "features": [feature] * (2 if wrap == "two" else 1)}
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    return doc


def write(tmp_path: Path, doc: dict[str, Any] | str, name: str = "d.geojson") -> Path:
    path = tmp_path / name
    path.write_text(doc if isinstance(doc, str) else json.dumps(doc))
    return path


def square(holes: bool = True) -> dict[str, Any]:
    return dict(mapping(Polygon(OUTER, [HOLE] if holes else [])))


def refused(domain: ModuleType, path: Path, *says: str, crs: str | None = None) -> None:
    with pytest.raises(domain.DomainError) as info:
        domain.read_domain(path, META, crs)
    assert isinstance(info.value, ValueError)
    for word in says:
        assert word in str(info.value), str(info.value)


class TestReading:
    @pytest.mark.parametrize("wrap", ["geometry", "feature", "collection"])
    def test_geojson_as_geometry_feature_or_one_feature_collection(
        self, domain: ModuleType, tmp_path: Path, wrap: str
    ) -> None:
        out = domain.read_domain(write(tmp_path, geojson(square(), wrap=wrap)), META)
        assert out.epsg == 25833
        assert isinstance(out.polygon, Polygon)
        assert len(out.polygon.interiors) == 1

    @pytest.mark.parametrize("crs", [UTM33, "EPSG:25833"])
    def test_both_spellings_of_the_crs_member(
        self, domain: ModuleType, tmp_path: Path, crs: str
    ) -> None:
        assert domain.read_domain(write(tmp_path, geojson(square(), crs=crs)), META).epsg == 25833

    @pytest.mark.parametrize("flag", ["EPSG:25833", "epsg:25833"])
    @pytest.mark.parametrize("member", [UTM33, "EPSG:25833"])
    def test_geojson_with_an_agreeing_domain_crs(
        self, domain: ModuleType, tmp_path: Path, member: str, flag: str
    ) -> None:
        path = write(tmp_path, geojson(square(), crs=member))
        assert domain.read_domain(path, META, flag).epsg == 25833

    def test_the_json_suffix_is_geojson(self, domain: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path, geojson(square()), name="d.json")
        assert domain.read_domain(path, META).epsg == 25833

    def test_wkt_with_a_crs(self, domain: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path, Polygon(OUTER, [HOLE]).wkt, name="d.wkt")
        out = domain.read_domain(path, META, "EPSG:25833")
        assert out.epsg == 25833
        assert len(out.polygon.interiors) == 1

    def test_vertices_are_used_as_given(self, domain: ModuleType, tmp_path: Path) -> None:
        """Directions 3 and 6: no densifying, simplifying or snapping."""
        out = domain.read_domain(write(tmp_path, geojson(square())), META)
        assert set(out.polygon.exterior.coords) == set(OUTER)
        assert set(out.polygon.interiors[0].coords) == set(HOLE)
        assert len(out.polygon.exterior.coords) == len(OUTER) + 1

    @pytest.mark.parametrize("reverse", [False, True])
    def test_orientation_is_outer_ccw_and_holes_cw(
        self, domain: ModuleType, tmp_path: Path, reverse: bool
    ) -> None:
        outer, hole = (OUTER[::-1], HOLE[::-1]) if reverse else (OUTER, HOLE)
        doc = geojson(dict(mapping(Polygon(outer, [hole]))))
        out = domain.read_domain(write(tmp_path, doc), META)
        assert out.polygon.exterior.is_ccw
        assert not out.polygon.interiors[0].is_ccw

    def test_a_vertex_on_the_border_of_the_node_rectangle_is_inside(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        corners = [
            (500_000.0, 6_599_970.0),
            (500_080.0, 6_599_970.0),
            (500_080.0, 6_600_000.0),
            (500_000.0, 6_600_000.0),
        ]
        out = domain.read_domain(write(tmp_path, geojson(dict(mapping(Polygon(corners))))), META)
        assert out.polygon.bounds == (500_000.0, 6_599_970.0, 500_080.0, 6_600_000.0)

    def test_the_model_is_frozen(self, domain: ModuleType, tmp_path: Path) -> None:
        out = domain.read_domain(write(tmp_path, geojson(square())), META)
        # Pydantic's ValidationError is a ValueError; a frozen dataclass's
        # FrozenInstanceError is an AttributeError.
        with pytest.raises((ValueError, AttributeError, TypeError)):
            out.epsg = 4326


class TestRefusals:
    def test_a_multipolygon(self, domain: ModuleType, tmp_path: Path) -> None:
        left = Polygon(
            [(500_005.0, 6_599_975.0), (500_020.0, 6_599_975.0), (500_020.0, 6_599_990.0)]
        )
        right = Polygon(
            [(500_050.0, 6_599_975.0), (500_070.0, 6_599_975.0), (500_070.0, 6_599_990.0)]
        )
        multi = dict(mapping(shapely.MultiPolygon([left, right])))
        refused(domain, write(tmp_path, geojson(multi)), "MultiPolygon")

    def test_a_collection_of_two_features(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, geojson(square(), wrap="two")))

    @pytest.mark.parametrize(
        "wkt", ["POINT (500010 6599980)", "LINESTRING (500010 6599980, 500020 6599990)"]
    )
    def test_not_a_polygon(self, domain: ModuleType, tmp_path: Path, wkt: str) -> None:
        refused(domain, write(tmp_path, wkt, name="d.wkt"), crs="EPSG:25833")

    def test_empty(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, "POLYGON EMPTY", name="d.wkt"), "empty", crs="EPSG:25833")

    def test_invalid_names_the_reason(self, domain: ModuleType, tmp_path: Path) -> None:
        bowtie = [
            (500_010.0, 6_599_975.0),
            (500_060.0, 6_599_995.0),
            (500_060.0, 6_599_975.0),
            (500_010.0, 6_599_995.0),
        ]
        doc = geojson({"type": "Polygon", "coordinates": [[*bowtie, bowtie[0]]]})
        refused(domain, write(tmp_path, doc), "Self-intersection")

    def test_geojson_without_crs_is_wgs84_by_the_standard(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        refused(domain, write(tmp_path, geojson(square(), crs=None)), "4326", "25833")

    def test_a_mismatched_epsg(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, geojson(square(), crs="EPSG:25832")), "25832", "25833")

    @pytest.mark.parametrize(
        ("member", "flag", "says"),
        [
            (UTM33, "EPSG:25832", ("25833", "25832", "--domain-crs")),
            ("EPSG:25832", "EPSG:25833", ("25832", "25833", "--domain-crs")),
            (None, "EPSG:25833", ("4326", "25833", "--domain-crs")),
        ],
    )
    def test_geojson_with_a_disagreeing_domain_crs(
        self,
        domain: ModuleType,
        tmp_path: Path,
        member: str | None,
        flag: str,
        says: tuple[str, ...],
    ) -> None:
        """The flag never overrides the file's own CRS: a disagreement is refused.

        The third row is a file with no ``crs`` member, which RFC 7946 makes
        WGS 84; ``--domain-crs`` does not supply one."""
        refused(domain, write(tmp_path, geojson(square(), crs=member)), *says, crs=flag)

    def test_wkt_without_a_crs(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, Polygon(OUTER).wkt, name="d.wkt"))

    def test_wkt_with_a_mismatched_crs(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, Polygon(OUTER).wkt, name="d.wkt"), "4326", crs="EPSG:4326")

    @pytest.mark.parametrize("dx", [1.0, 1e-6])
    def test_a_vertex_outside_the_node_rectangle(
        self, domain: ModuleType, tmp_path: Path, dx: float
    ) -> None:
        """U4 (a): refused, naming the vertex and saying it is outside."""
        outer = [*OUTER[:1], (500_080.0 + dx, 6_599_975.0), *OUTER[2:]]
        refused(domain, write(tmp_path, geojson(dict(mapping(Polygon(outer))))), "outside")

    def test_an_unknown_suffix(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, geojson(square()), name="d.shp"), ".shp")

    def test_malformed_json_is_a_domain_error(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, '{"type": "Polygon", '))

    def test_malformed_wkt_is_a_domain_error(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, "POLYGON ((1 2, 3", name="d.wkt"), crs="EPSG:25833")
