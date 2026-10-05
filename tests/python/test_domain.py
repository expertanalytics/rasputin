"""`tin_engine.domain`: reading the domain polygon, and moving it into the DEM's CRS.

Increment 16, R1, U1 (a) and U4 (a), as changed by increment 15b
(`docs/increments/15-dem-mosaic.md` R9): the domain keeps its own CRS, and
`DomainPolygon.to_crs` replaces 16's must-match rule. Signatures, as pinned by
the 15b red suite ("Pinned by the red suite (15b)")::

    read_domain(path: Path, crs: str | None = None) -> DomainPolygon
    DomainPolygon.to_crs(dst: str | pyproj.CRS) -> DomainPolygon
    check_extent(domain: DomainPolygon, meta: RasterMeta) -> None

`read_domain` no longer takes the DEM: the domain is read before the tiles are
planned, because its bounds in the DEM's CRS choose them (R6). The extent
check runs after the transform, in `check_extent`, against the mosaic's node
rectangle (and the coverage check, in `open_dem`).

`DomainPolygon` is frozen with `.polygon` (a shapely `Polygon`, oriented outer
counter-clockwise and holes clockwise) and `.crs` (`str`, text pyproj parses to
the domain's CRS; it was `.epsg: int` in 16). Every refusal raises
`DomainError`, a `ValueError` subclass, whose message the CLI shows. `crs` is
the `--domain-crs` text, anything pyproj parses; it is what a `.wkt` file
needs, and for GeoJSON it must agree with the file's own (16, unchanged).

The module is fetched inside a fixture, so its absence fails these tests and
leaves the rest of the session collecting.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely
from pyproj import CRS, Transformer
from shapely.geometry import Polygon, mapping

from crs_fixtures import proj4_of, refuse_point_moves
from tin_engine.io.models import RasterMeta

UTM33 = "urn:ogc:def:crs:EPSG::25833"
# UTM 33 with the easting axis pointing west: `to_epsg()` still answers 25833
# (at pyproj's default confidence), but it is not EPSG:25833, and the numbers
# differ in sign. A same-CRS test by EPSG code would skip its transform.
UTM33_WEST = "+proj=utm +zone=33 +ellps=GRS80 +units=m +axis=wnu +no_defs"

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
        domain.read_domain(path, crs)
    assert isinstance(info.value, ValueError)
    for word in says:
        assert word in str(info.value), str(info.value)


def transformed(src: str, dst: str, ring: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """pyproj's own `always_xy` transform of `ring`, vertex by vertex: the oracle."""
    t = Transformer.from_crs(src, dst, always_xy=True)
    xs, ys = t.transform(np.array([p[0] for p in ring]), np.array([p[1] for p in ring]))
    return [(float(x), float(y)) for x, y in zip(xs, ys, strict=True)]


def ring_set(ring: Any) -> set[tuple[float, float]]:
    return {(float(x), float(y)) for x, y in ring.coords}


def the_crs(out: Any) -> CRS:
    assert isinstance(out.crs, str)
    return CRS.from_user_input(out.crs)


@pytest.fixture
def no_transformer(monkeypatch: pytest.MonkeyPatch) -> None:
    """No point is moved: any point-moving `Transformer` method from here on
    fails the test (`crs_fixtures.refuse_point_moves`). Building one is
    allowed, since `crs.same_crs` builds one to compare (audit PR B)."""
    refuse_point_moves(monkeypatch)


class TestReading:
    @pytest.mark.parametrize("wrap", ["geometry", "feature", "collection"])
    def test_geojson_as_geometry_feature_or_one_feature_collection(
        self, domain: ModuleType, tmp_path: Path, wrap: str
    ) -> None:
        out = domain.read_domain(write(tmp_path, geojson(square(), wrap=wrap)))
        assert the_crs(out) == CRS.from_epsg(25833)
        assert isinstance(out.polygon, Polygon)
        assert len(out.polygon.interiors) == 1

    @pytest.mark.parametrize("crs", [UTM33, "EPSG:25833"])
    def test_both_spellings_of_the_crs_member(
        self, domain: ModuleType, tmp_path: Path, crs: str
    ) -> None:
        out = domain.read_domain(write(tmp_path, geojson(square(), crs=crs)))
        assert the_crs(out) == CRS.from_epsg(25833)

    @pytest.mark.parametrize("flag", ["EPSG:25833", "epsg:25833"])
    @pytest.mark.parametrize("member", [UTM33, "EPSG:25833"])
    def test_geojson_with_an_agreeing_domain_crs(
        self, domain: ModuleType, tmp_path: Path, member: str, flag: str
    ) -> None:
        path = write(tmp_path, geojson(square(), crs=member))
        assert the_crs(domain.read_domain(path, flag)) == CRS.from_epsg(25833)

    def test_the_json_suffix_is_geojson(self, domain: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path, geojson(square()), name="d.json")
        assert the_crs(domain.read_domain(path)) == CRS.from_epsg(25833)

    def test_wkt_with_a_crs(self, domain: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path, Polygon(OUTER, [HOLE]).wkt, name="d.wkt")
        out = domain.read_domain(path, "EPSG:25833")
        assert the_crs(out) == CRS.from_epsg(25833)
        assert len(out.polygon.interiors) == 1

    def test_geojson_without_crs_is_wgs84_by_the_standard(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        """RFC 7946: no `crs` member is WGS 84. 15b reads it as such, where 16
        refused it as a mismatch."""
        lon_lat = transformed("EPSG:25833", "EPSG:4326", OUTER)
        out = domain.read_domain(write(tmp_path, geojson(dict(mapping(Polygon(lon_lat))), None)))
        assert the_crs(out) == CRS.from_epsg(4326)
        assert ring_set(out.polygon.exterior) == set(lon_lat)

    def test_another_utm_zone_is_read_in_its_own_crs(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        out = domain.read_domain(write(tmp_path, geojson(square(), crs="EPSG:25832")))
        assert the_crs(out) == CRS.from_epsg(25832)

    def test_wkt_in_a_geographic_crs(self, domain: ModuleType, tmp_path: Path) -> None:
        lon_lat = transformed("EPSG:25833", "EPSG:4326", OUTER)
        out = domain.read_domain(write(tmp_path, Polygon(lon_lat).wkt, name="d.wkt"), "EPSG:4326")
        assert the_crs(out) == CRS.from_epsg(4326)

    def test_a_crs_without_an_epsg_code(self, domain: ModuleType, tmp_path: Path) -> None:
        """R9: `crs: str`, "since a domain may have no EPSG code". 16 refused it."""
        lcc = "+proj=lcc +lat_1=60 +lat_2=65 +lat_0=62 +lon_0=15 +ellps=GRS80 +units=m +no_defs"
        out = domain.read_domain(write(tmp_path, Polygon(OUTER).wkt, name="d.wkt"), lcc)
        assert the_crs(out) == CRS.from_user_input(lcc)

    def test_ogc_crs84_member(self, domain: ModuleType, tmp_path: Path) -> None:
        """Older GeoJSON writers name longitude-latitude WGS 84 this way."""
        lon_lat = transformed("EPSG:25833", "EPSG:4326", OUTER)
        doc = geojson(dict(mapping(Polygon(lon_lat))), "urn:ogc:def:crs:OGC:1.3:CRS84")
        out = domain.read_domain(write(tmp_path, doc))
        assert the_crs(out) == CRS.from_user_input("OGC:CRS84")

    def test_vertices_are_used_as_given(self, domain: ModuleType, tmp_path: Path) -> None:
        """Directions 3 and 6: no densifying, simplifying or snapping."""
        out = domain.read_domain(write(tmp_path, geojson(square())))
        assert set(out.polygon.exterior.coords) == set(OUTER)
        assert set(out.polygon.interiors[0].coords) == set(HOLE)
        assert len(out.polygon.exterior.coords) == len(OUTER) + 1

    @pytest.mark.parametrize("reverse", [False, True])
    def test_orientation_is_outer_ccw_and_holes_cw(
        self, domain: ModuleType, tmp_path: Path, reverse: bool
    ) -> None:
        outer, hole = (OUTER[::-1], HOLE[::-1]) if reverse else (OUTER, HOLE)
        doc = geojson(dict(mapping(Polygon(outer, [hole]))))
        out = domain.read_domain(write(tmp_path, doc))
        assert out.polygon.exterior.is_ccw
        assert not out.polygon.interiors[0].is_ccw

    def test_a_vertex_outside_any_dem_is_not_reading_s_business(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        """The extent is checked after the transform, against the tiles (R9)."""
        far = [(x + 1e6, y) for x, y in OUTER]
        out = domain.read_domain(write(tmp_path, geojson(dict(mapping(Polygon(far))))))
        assert ring_set(out.polygon.exterior) == set(far)

    def test_the_model_is_frozen(self, domain: ModuleType, tmp_path: Path) -> None:
        out = domain.read_domain(write(tmp_path, geojson(square())))
        # Pydantic's ValidationError is a ValueError; a frozen dataclass's
        # FrozenInstanceError is an AttributeError.
        with pytest.raises((ValueError, AttributeError, TypeError)):
            out.crs = "EPSG:4326"

    def test_the_must_match_rule_is_gone(self, domain: ModuleType) -> None:
        """R9: `check_crs` is replaced by the transform, and `epsg: int` by `crs: str`."""
        assert not hasattr(domain, "check_crs")
        assert "epsg" not in domain.DomainPolygon.model_fields
        assert "crs" in domain.DomainPolygon.model_fields


class TestToCrs:
    """R9: vertices only, each through pyproj's `always_xy` transform, exactly once."""

    def lon_lat_domain(self, domain: ModuleType, tmp_path: Path) -> Any:
        outer = transformed("EPSG:25833", "EPSG:4326", OUTER)
        hole = transformed("EPSG:25833", "EPSG:4326", HOLE)
        doc = geojson(dict(mapping(Polygon(outer, [hole]))), crs=None)
        return domain.read_domain(write(tmp_path, doc))

    def test_wgs84_to_utm33_is_pyprojs_transform_of_every_vertex(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        read = self.lon_lat_domain(domain, tmp_path)
        out = read.to_crs("EPSG:25833")
        assert the_crs(out) == CRS.from_epsg(25833)
        for before, after in (
            (read.polygon.exterior, out.polygon.exterior),
            (read.polygon.interiors[0], out.polygon.interiors[0]),
        ):
            expected = transformed("EPSG:4326", "EPSG:25833", list(before.coords))
            assert ring_set(after) == set(expected)
            # Vertices only: not densified, nothing dropped.
            assert len(after.coords) == len(before.coords)

    def test_the_result_lands_where_the_utm_polygon_was(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        """A round trip through longitude and latitude is not exact, but it is
        within a micrometre; a swapped axis would be thousands of km off."""
        out = self.lon_lat_domain(domain, tmp_path).to_crs("EPSG:25833")
        got = sorted(ring_set(out.polygon.exterior))
        for (x, y), (ex, ey) in zip(got, sorted(set(OUTER)), strict=True):
            assert x == pytest.approx(ex, abs=1e-6)
            assert y == pytest.approx(ey, abs=1e-6)

    def test_the_source_is_unchanged(self, domain: ModuleType, tmp_path: Path) -> None:
        read = self.lon_lat_domain(domain, tmp_path)
        before = list(read.polygon.exterior.coords)
        read.to_crs("EPSG:25833")
        assert list(read.polygon.exterior.coords) == before
        assert the_crs(read) == CRS.from_epsg(4326)

    @pytest.mark.parametrize("src", ["EPSG:25832", "EPSG:3035"])
    def test_projected_to_projected(self, domain: ModuleType, tmp_path: Path, src: str) -> None:
        ring = transformed("EPSG:25833", src, OUTER)
        read = domain.read_domain(write(tmp_path, Polygon(ring).wkt, name="d.wkt"), src)
        out = read.to_crs("EPSG:25833")
        assert ring_set(out.polygon.exterior) == set(transformed(src, "EPSG:25833", ring))

    def test_a_crs_object_as_the_destination(self, domain: ModuleType, tmp_path: Path) -> None:
        read = self.lon_lat_domain(domain, tmp_path)
        by_text = read.to_crs("EPSG:25833")
        by_object = read.to_crs(CRS.from_epsg(25833))
        assert list(by_object.polygon.exterior.coords) == list(by_text.polygon.exterior.coords)

    def test_the_winding_contract_survives_a_mirroring_transform(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        """From a west-pointing easting, the transform mirrors the ring; the
        result is re-oriented, outer counter-clockwise and holes clockwise."""
        outer = [(-x, y) for x, y in OUTER]
        hole = [(-x, y) for x, y in HOLE]
        path = write(tmp_path, Polygon(outer, [hole]).wkt, name="d.wkt")
        out = domain.read_domain(path, UTM33_WEST).to_crs("EPSG:25833")
        assert out.polygon.exterior.is_ccw
        assert not out.polygon.interiors[0].is_ccw
        assert ring_set(out.polygon.exterior) == set(transformed(UTM33_WEST, "EPSG:25833", outer))

    def test_the_same_epsg_code_is_not_the_same_crs(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        """Same-CRS is decided by CRS equality, not by `to_epsg()`: this CRS's
        `to_epsg()` is 25833, and skipping its transform would put the domain
        a thousand kilometres west of where it is."""
        assert CRS.from_user_input(UTM33_WEST).to_epsg() == 25833  # the trap is real
        outer = [(-x, y) for x, y in OUTER]
        out = domain.read_domain(write(tmp_path, Polygon(outer).wkt, name="d.wkt"), UTM33_WEST)
        moved = out.to_crs("EPSG:25833")
        assert ring_set(moved.polygon.exterior) == set(transformed(UTM33_WEST, "EPSG:25833", outer))
        assert all(x > 0 for x, _ in moved.polygon.exterior.coords)

    @pytest.mark.parametrize("dst", ["EPSG:25833", "urn:ogc:def:crs:EPSG::25833"])
    def test_the_same_crs_is_not_transformed(
        self, domain: ModuleType, tmp_path: Path, no_transformer: None, dst: str
    ) -> None:
        """A domain already in the DEM's CRS keeps its coordinates bit for bit,
        and no point is moved (the mesh must be 16's, bit for bit)."""
        read = domain.read_domain(write(tmp_path, geojson(square(), crs=UTM33)))
        out = read.to_crs(dst)
        assert list(out.polygon.exterior.coords) == list(read.polygon.exterior.coords)
        assert list(out.polygon.interiors[0].coords) == list(read.polygon.interiors[0].coords)
        assert the_crs(out) == CRS.from_epsg(25833)

    def test_a_proj_string_of_the_targets_epsg_code_is_not_transformed(
        self, domain: ModuleType, no_transformer: None
    ) -> None:
        """Audit PR B, red test 8: a domain spelt as the PROJ string of
        EPSG:25833 is in EPSG:25833 (`crs.same_crs`), so no point moves, the
        polygon is exactly the given one, and it is labelled as `dst`."""
        given = domain.DomainPolygon(polygon=Polygon(OUTER, [HOLE]), crs=proj4_of(25833))
        out = given.to_crs("EPSG:25833")
        assert out.crs == "EPSG:25833"
        assert out.polygon.equals_exact(given.polygon, tolerance=0.0)

    def test_a_vertex_with_no_image_is_refused_naming_both_crss(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        """Latitude 95 has no image in UTM 33 (pyproj answers inf). Refused,
        at reading or at the transform, never passed on as a coordinate."""
        ring = [(14.9, 59.4), (15.1, 59.4), (15.1, 95.0), (14.9, 59.6)]
        path = write(tmp_path, geojson(dict(mapping(Polygon(ring))), crs=None))
        with pytest.raises(domain.DomainError) as info:
            domain.read_domain(path).to_crs("EPSG:25833")
        for word in ("4326", "25833"):
            assert word in str(info.value), info.value

    def test_utm_numbers_in_a_file_without_a_crs_member_are_refused(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        """The commonest wrong file: UTM coordinates, no `crs` member, so WGS 84
        by the standard. Its "longitudes" are 500 000."""
        path = write(tmp_path, geojson(square(), crs=None))
        with pytest.raises(domain.DomainError) as info:
            domain.read_domain(path).to_crs("EPSG:25833")
        for word in ("4326", "25833"):
            assert word in str(info.value), info.value


class TestCheckExtent:
    """U4 (a), now after the transform, against the mosaic's node rectangle."""

    def test_a_vertex_on_the_border_of_the_node_rectangle_is_inside(
        self, domain: ModuleType, tmp_path: Path
    ) -> None:
        corners = [
            (500_000.0, 6_599_970.0),
            (500_080.0, 6_599_970.0),
            (500_080.0, 6_600_000.0),
            (500_000.0, 6_600_000.0),
        ]
        out = domain.read_domain(write(tmp_path, geojson(dict(mapping(Polygon(corners))))))
        assert domain.check_extent(out, META) is None

    @pytest.mark.parametrize("dx", [1.0, 1e-6])
    def test_a_vertex_outside_the_node_rectangle(
        self, domain: ModuleType, tmp_path: Path, dx: float
    ) -> None:
        """Refused, naming the vertex and saying it is outside."""
        outer = [*OUTER[:1], (500_080.0 + dx, 6_599_975.0), *OUTER[2:]]
        out = domain.read_domain(write(tmp_path, geojson(dict(mapping(Polygon(outer))))))
        with pytest.raises(domain.DomainError) as info:
            domain.check_extent(out, META)
        assert "outside" in str(info.value)
        assert "500080" in str(info.value).replace(",", "").replace(" ", "")

    def test_a_hole_vertex_is_checked_too(self, domain: ModuleType, tmp_path: Path) -> None:
        corners = [
            (500_000.0, 6_599_970.0),
            (500_080.0, 6_599_970.0),
            (500_080.0, 6_600_000.0),
            (500_000.0, 6_600_000.0),
        ]
        out = domain.read_domain(write(tmp_path, geojson(dict(mapping(Polygon(corners, [HOLE]))))))
        assert domain.check_extent(out, META) is None


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
        """The flag never overrides the file's own CRS: a disagreement is refused
        (16, unchanged by 15b: it is the file against the flag, not the DEM).

        The third row is a file with no ``crs`` member, which RFC 7946 makes
        WGS 84; ``--domain-crs`` does not supply one."""
        refused(domain, write(tmp_path, geojson(square(), crs=member)), *says, crs=flag)

    def test_wkt_without_a_crs(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, Polygon(OUTER).wkt, name="d.wkt"))

    @pytest.mark.parametrize("bad", ["EPSG:999999", "not a crs"])
    def test_an_unknown_domain_crs_is_named(
        self, domain: ModuleType, tmp_path: Path, bad: str
    ) -> None:
        refused(domain, write(tmp_path, Polygon(OUTER).wkt, name="d.wkt"), bad, crs=bad)

    def test_an_unknown_crs_member_is_named(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, geojson(square(), crs="EPSG:999999")), "EPSG:999999")

    def test_an_unknown_suffix(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, geojson(square()), name="d.shp"), ".shp")

    def test_malformed_json_is_a_domain_error(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, '{"type": "Polygon", '))

    def test_malformed_wkt_is_a_domain_error(self, domain: ModuleType, tmp_path: Path) -> None:
        refused(domain, write(tmp_path, "POLYGON ((1 2, 3", name="d.wkt"), crs="EPSG:25833")
