"""`rasputin network-catchment` and its MCP tool, on networks written by the test.

`docs/increments/85-network-catchment.md`, sections 3 and 7.2. The network
is test 2's confluence plus a unit beside it, as 1 km squares in EPSG:3006
(x from 500 000, y from 7 000 000, 1 km2 each):

    SE200 (Bravo)                 north of SE100, drains to SE100
    SE100 (Alpha) SE300 (Charlie) SE400 (Delta)
    SE500 (Echo)                  south of SE100, drains to the end

SE300 drains to SE100, SE400 to SE300, SE100 and SE500 to the end, -9999.
Above SE100: four units, 4 km2. Each unit's outlet point is the midpoint of
its border with the unit it drains to; SE100's is (500 500, 7 000 000).

Also here, the outlet radius comparison of section 7.1's mutation targets
(`<=` against `<`), since section 5 gives that comparison no pure function:
a seed at 199.9 m and 200.0 m is accepted at the default radius, 200.1 m is
refused. Distances in EPSG:3006 are planar, so these are exact to the
spacing of a double at 7e6, about 1e-9 m.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import pytest
import shapely
from pyproj import Transformer
from shapely.geometry import Point, Polygon, box, mapping

from catchment_fixtures import bowl, lake_box, lat, tile_of
from cli_helpers import invoke, plain, ran
from gpkg_fixtures import HAS_RTREE, Layer, Row, write_gpkg
from mosaic_fixtures import quadrants
from test_cli_mesh_mosaic import write_tiles
from tin_engine.crs import parse_crs
from tin_engine.io.domain_file import read_domain

COMMAND = "network-catchment"
CRS = "EPSG:3006"
X0 = 500_000.0
Y0 = 7_000_000.0
KM = 1000.0
KM2 = 1e6
END = "-9999"
ID, DOWN, AREA, NAME = "unit_code", "down", "up_area", "name"


def square(i: float, j: float) -> Polygon:
    x, y = X0 + i * KM, Y0 + j * KM
    return box(x, y, x + KM, y + KM)


def at(dx: float, dy: float) -> tuple[float, float]:
    return (X0 + dx, Y0 + dy)


#: (id, down, stated upstream area in m2, name, geometry)
Unit = tuple[Any, Any, Any, Any, Any]

NETWORK: list[Unit] = [
    ("SE100", END, 4 * KM2, "Alpha", square(0, 0)),
    ("SE200", "SE100", 1 * KM2, "Bravo", square(0, 1)),
    ("SE300", "SE100", 2 * KM2, "Charlie", square(1, 0)),
    ("SE400", "SE300", 1 * KM2, "Delta", square(2, 0)),
    ("SE500", END, 1 * KM2, "Echo", square(0, -1)),
]

#: (id, x, y): each unit's outlet point.
OUTLETS: list[tuple[Any, float, float]] = [
    ("SE100", *at(500, 0)),
    ("SE200", *at(500, KM)),
    ("SE300", *at(KM, 500)),
    ("SE400", *at(2 * KM, 500)),
    ("SE500", *at(500, -KM)),
]

ABOVE_SE100 = ["SE100", "SE200", "SE300", "SE400"]


def feature(geometry: Any, properties: dict[str, Any]) -> dict[str, Any]:
    return {"type": "Feature", "properties": properties, "geometry": mapping(geometry)}


def collection(path: Path, features: list[dict[str, Any]], crs: str | None = CRS) -> Path:
    doc: dict[str, Any] = {"type": "FeatureCollection", "features": features}
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(doc))
    return path


def network_geojson(path: Path, units: Sequence[Unit] = NETWORK, crs: str | None = CRS) -> Path:
    return collection(
        path,
        [
            feature(g, {ID: own, DOWN: down, AREA: area, NAME: name})
            for own, down, area, name, g in units
        ],
        crs,
    )


def network_gpkg(path: Path, units: Sequence[Unit] = NETWORK) -> Path:
    rows = [
        Row(pk, g, {ID: own, DOWN: down, AREA: area, NAME: name})
        for pk, (own, down, area, name, g) in enumerate(units, 1)
    ]
    layer = Layer(table="units", srs_id=3006, rows=rows, columns=(ID, DOWN, AREA, NAME))
    return write_gpkg(path, [layer])


def outlets_geojson(
    path: Path, outlets: Sequence[tuple[Any, float, float]] = OUTLETS, key: str = ID
) -> Path:
    return collection(path, [feature(Point(x, y), {key: own}) for own, x, y in outlets])


def lakes_geojson(path: Path, lake: Polygon, crs: str = CRS) -> Path:
    return collection(path, [feature(lake, {"name": "the lake"})], crs)


ATTRIBUTES = ("--id-attribute", ID, "--down-attribute", DOWN)
STATED = ("--upstream-area-attribute", AREA, "--area-unit", "m2")


def utm(x: float, y: float) -> tuple[str, ...]:
    return ("--seed", repr(x), repr(y), "--seed-crs", CRS)


@pytest.fixture(params=["geojson", "gpkg"])
def network(request: pytest.FixtureRequest, tmp_path: Path) -> tuple[str, ...]:
    """The network's file options, as GeoJSON and as a GeoPackage with its R-tree."""
    if request.param == "gpkg":
        if not HAS_RTREE:
            pytest.skip("this SQLite has no R-tree module")
        return ("--network", str(network_gpkg(tmp_path / "units.gpkg")), "--network-layer", "units")
    return ("--network", str(network_geojson(tmp_path / "units.geojson")))


@pytest.fixture
def geojson_network(tmp_path: Path) -> tuple[str, ...]:
    return ("--network", str(network_geojson(tmp_path / "units.geojson")))


@pytest.fixture
def outlets(tmp_path: Path) -> tuple[str, ...]:
    return ("--outlets", str(outlets_geojson(tmp_path / "outlets.geojson")))


def run(tmp_path: Path, *args: str, out: str = "n.geojson") -> tuple[dict[str, Any], str]:
    """Run the command, required to succeed; the written properties and the output."""
    target = tmp_path / out
    output = ran(COMMAND, *args, "--out", str(target))
    doc = json.loads(target.read_text())
    (only,) = doc["features"]
    return only["properties"], output


def refuse(tmp_path: Path, *args: str, says: tuple[str, ...]) -> str:
    """Refused with exit 2, saying each of `says`, writing nothing."""
    target = tmp_path / "refused.geojson"
    code, output = invoke(COMMAND, *args, "--out", str(target))
    assert code == 2, output
    assert "Traceback" not in output
    assert "No such option" not in output and "No such command" not in output, output
    for word in says:
        assert word in output, f"{word!r} not in {output!r}"
    assert not target.exists()
    return output


# ---------------------------------------------------------------------------
# The happy path: --unit
# ---------------------------------------------------------------------------


class TestUnitSeed:
    def test_the_unit_and_everything_upstream(
        self, tmp_path: Path, network: tuple[str, ...]
    ) -> None:
        props, _ = run(tmp_path, *network, *ATTRIBUTES, *STATED, "--unit", "SE100")
        assert props["units"] == ABOVE_SE100
        assert props["unit_count"] == 4
        assert props["area_km2"] == pytest.approx(4.0, rel=1e-9)
        assert props["stated_upstream_area_km2"] == pytest.approx(4.0, rel=1e-9)

    def test_a_unit_higher_up(self, tmp_path: Path, network: tuple[str, ...]) -> None:
        props, _ = run(tmp_path, *network, *ATTRIBUTES, *STATED, "--unit", "SE300")
        assert props["units"] == ["SE300", "SE400"]
        assert props["area_km2"] == pytest.approx(2.0, rel=1e-9)

    def test_the_domain_reader_reads_it_in_the_networks_crs(
        self, tmp_path: Path, network: tuple[str, ...]
    ) -> None:
        run(tmp_path, *network, *ATTRIBUTES, "--unit", "SE100")
        domain = read_domain(tmp_path / "n.geojson")
        assert parse_crs(domain.crs) == parse_crs(CRS)
        assert domain.polygon.is_valid
        assert domain.polygon.area == pytest.approx(4 * KM2, rel=1e-9)
        expected = shapely.union_all([g for *_, g in NETWORK[:4]])
        assert domain.polygon.symmetric_difference(expected).area < 1e-6  # m2, at 4e6 m2

    def test_the_path_is_printed(self, tmp_path: Path, geojson_network: tuple[str, ...]) -> None:
        _, output = run(tmp_path, *geojson_network, *ATTRIBUTES, "--unit", "SE100")
        assert str(tmp_path / "n.geojson") in output

    def test_geopackage_and_geojson_give_the_same_file(self, tmp_path: Path) -> None:
        if not HAS_RTREE:
            pytest.skip("this SQLite has no R-tree module")
        common = (*ATTRIBUTES, *STATED, "--name-attribute", NAME, "--unit", "SE100")
        gpkg = ("--network", str(network_gpkg(tmp_path / "u.gpkg")), "--network-layer", "units")
        run(tmp_path, *gpkg, *common, out="a.geojson")
        run(
            tmp_path,
            "--network",
            str(network_geojson(tmp_path / "u.geojson")),
            *common,
            out="b.geojson",
        )
        a = json.loads((tmp_path / "a.geojson").read_text())
        b = json.loads((tmp_path / "b.geojson").read_text())
        assert a == b

    def test_no_stated_area_attribute_is_null_and_unchecked(
        self, tmp_path: Path, geojson_network: tuple[str, ...]
    ) -> None:
        props, _ = run(tmp_path, *geojson_network, *ATTRIBUTES, "--unit", "SE100")
        assert props["stated_upstream_area_km2"] is None
        assert props["area_km2"] == pytest.approx(4.0, rel=1e-9)

    def test_a_stated_area_in_km2(self, tmp_path: Path) -> None:
        units = [(o, d, a / KM2, n, g) for o, d, a, n, g in NETWORK]
        path = network_geojson(tmp_path / "km2.geojson", units)
        props, _ = run(
            tmp_path,
            "--network",
            str(path),
            *ATTRIBUTES,
            "--upstream-area-attribute",
            AREA,
            "--area-unit",
            "km2",
            "--unit",
            "SE100",
        )
        assert props["stated_upstream_area_km2"] == pytest.approx(4.0, rel=1e-9)

    def test_numeric_ids_compare_as_text(self, tmp_path: Path) -> None:
        """A HydroBASINS-style network: integer IDs, written as floats, the end 0."""
        units = [(7.0, 0, None, None, square(0, 0)), (8, 7.0, None, None, square(1, 0))]
        path = network_geojson(tmp_path / "numeric.geojson", units)
        props, _ = run(tmp_path, "--network", str(path), *ATTRIBUTES, "--unit", "7")
        assert props["units"] == ["7", "8"]

    def test_holes_are_subtracted_from_the_area_and_filled_in_the_file(
        self, tmp_path: Path
    ) -> None:
        """Eight units round a missing ninth, chained to one outlet unit; the
        stated area is the eight units' own (a hole belongs to no unit), so the
        area check passes on the area with the hole subtracted."""
        cells = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
        units: list[Unit] = []
        for k, (i, j) in enumerate(cells):
            down = END if k == 0 else f"R{k - 1}"
            units.append((f"R{k}", down, (len(cells) - k) * KM2, None, square(i, j)))
        path = network_geojson(tmp_path / "ring.geojson", units)
        props, _ = run(tmp_path, "--network", str(path), *ATTRIBUTES, *STATED, "--unit", "R0")
        assert props["unit_count"] == 8
        assert props["area_km2"] == pytest.approx(8.0, rel=1e-9)
        domain = read_domain(tmp_path / "n.geojson")
        assert len(domain.polygon.interiors) == 0
        assert domain.polygon.area == pytest.approx(9 * KM2, rel=1e-9)

    def test_a_network_without_a_crs_member_is_in_degrees(self, tmp_path: Path) -> None:
        to = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
        units = [
            (o, d, a, n, shapely.transform(g, lambda x, y: to.transform(x, y), interleaved=False))
            for o, d, a, n, g in NETWORK
        ]
        path = network_geojson(tmp_path / "lonlat.geojson", units, crs=None)
        props, _ = run(tmp_path, "--network", str(path), *ATTRIBUTES, *STATED, "--unit", "SE100")
        # 0.1 %: SWEREF99 TM's scale on its central meridian is 0.9996 (0.08 % in area).
        assert props["area_km2"] == pytest.approx(4.0, rel=1e-3)
        assert parse_crs(read_domain(tmp_path / "n.geojson").crs) == parse_crs("EPSG:4326")


# ---------------------------------------------------------------------------
# The happy path: --seed at the home unit's outlet
# ---------------------------------------------------------------------------


class TestPointSeed:
    @pytest.mark.parametrize("metres", [150.0, 199.9, 200.0])
    def test_a_point_within_the_radius_of_its_units_outlet(
        self, tmp_path: Path, network: tuple[str, ...], outlets: tuple[str, ...], metres: float
    ) -> None:
        props, _ = run(tmp_path, *network, *ATTRIBUTES, *outlets, *utm(*at(500, metres)))
        assert props["units"] == ABOVE_SE100

    def test_a_point_just_beyond_the_radius_is_refused(
        self, tmp_path: Path, network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        refuse(
            tmp_path,
            *network,
            *ATTRIBUTES,
            *outlets,
            *utm(*at(500, 200.1)),
            says=("Invalid value for --seed", "SE100"),
        )

    def test_the_radius_is_an_option(
        self, tmp_path: Path, geojson_network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        seed = utm(*at(500, 150))
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            *outlets,
            *seed,
            "--outlet-radius",
            "100",
            says=("Invalid value for --seed",),
        )
        props, _ = run(
            tmp_path, *geojson_network, *ATTRIBUTES, *outlets, *seed, "--outlet-radius", "150"
        )
        assert props["units"] == ABOVE_SE100

    def test_a_seed_in_another_crs(
        self, tmp_path: Path, geojson_network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        lon, lat_ = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True).transform(*at(500, 150))
        props, _ = run(
            tmp_path, *geojson_network, *ATTRIBUTES, *outlets, "--seed", repr(lon), repr(lat_)
        )
        assert props["units"] == ABOVE_SE100

    def test_the_lake_under_the_point_brings_it_within_the_radius(
        self, tmp_path: Path, network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        """The point is 250 m from SE100's outlet; the lake under it, 85 m."""
        seed = utm(*at(500, 250))
        refuse(tmp_path, *network, *ATTRIBUTES, *outlets, *seed, says=("Invalid value for --seed",))
        lake = lakes_geojson(tmp_path / "lakes.geojson", box(*at(400, 85), *at(600, 300)))
        props, _ = run(tmp_path, *network, *ATTRIBUTES, *outlets, *seed, "--lakes", str(lake))
        assert props["units"] == ABOVE_SE100

    def test_the_refusal_names_the_unit_and_the_two_faithful_choices(
        self, tmp_path: Path, geojson_network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        said = refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            *outlets,
            "--name-attribute",
            NAME,
            *utm(*at(500, 500)),
            says=("Invalid value for --seed", "SE100", "Alpha", "500", "--unit SE100", "terrain"),
        )
        assert "km2" in said, said  # the home unit's own area

    def test_a_point_seed_needs_an_outlet_file(
        self, tmp_path: Path, geojson_network: tuple[str, ...]
    ) -> None:
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            *utm(*at(500, 50)),
            says=("Invalid value for --seed", "--outlets", "--unit"),
        )

    def test_no_outlet_point_for_the_home_unit(
        self, tmp_path: Path, geojson_network: tuple[str, ...]
    ) -> None:
        path = outlets_geojson(tmp_path / "o.geojson", OUTLETS[1:])
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            "--outlets",
            str(path),
            *utm(*at(500, 50)),
            says=("Invalid value for --seed", "SE100"),
        )

    def test_two_outlet_points_for_one_unit(
        self, tmp_path: Path, geojson_network: tuple[str, ...]
    ) -> None:
        path = outlets_geojson(tmp_path / "o.geojson", [*OUTLETS, ("SE100", *at(500, 10))])
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            "--outlets",
            str(path),
            *utm(*at(500, 50)),
            says=("Invalid value for --outlets", "SE100"),
        )


# ---------------------------------------------------------------------------
# Section 3.4: one test per refusal
# ---------------------------------------------------------------------------


class TestOptionRefusals:
    def test_both_unit_and_seed(
        self, tmp_path: Path, geojson_network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            *outlets,
            "--unit",
            "SE100",
            *utm(*at(500, 50)),
            says=("Invalid value for --unit", "--seed"),
        )

    def test_neither_unit_nor_seed(self, tmp_path: Path, geojson_network: tuple[str, ...]) -> None:
        refuse(tmp_path, *geojson_network, *ATTRIBUTES, says=("Invalid value for --unit", "--seed"))

    @pytest.mark.parametrize("option", ["--lakes", "--outlets"])
    def test_a_seed_file_without_seed(
        self, tmp_path: Path, geojson_network: tuple[str, ...], option: str
    ) -> None:
        path = outlets_geojson(tmp_path / "o.geojson")
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            "--unit",
            "SE100",
            option,
            str(path),
            says=(f"Invalid value for {option}", "--seed"),
        )

    @pytest.mark.parametrize(
        ("option", "value"), [("--outlet-radius", "300"), ("--check-dem", "glo30")]
    )
    def test_a_seed_option_without_seed(
        self, tmp_path: Path, geojson_network: tuple[str, ...], option: str, value: str
    ) -> None:
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            "--unit",
            "SE100",
            option,
            value,
            says=(f"Invalid value for {option}", "--seed"),
        )

    def test_area_unit_without_the_area_attribute(
        self, tmp_path: Path, geojson_network: tuple[str, ...]
    ) -> None:
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            "--area-unit",
            "m2",
            "--unit",
            "SE100",
            says=("Invalid value for --area-unit", "--upstream-area-attribute"),
        )

    def test_the_area_attribute_without_its_unit(
        self, tmp_path: Path, geojson_network: tuple[str, ...]
    ) -> None:
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            "--upstream-area-attribute",
            AREA,
            "--unit",
            "SE100",
            says=("Invalid value for --upstream-area-attribute", "--area-unit"),
        )

    @pytest.mark.parametrize(
        ("layer", "file"), [("--lakes-layer", "--lakes"), ("--outlets-layer", "--outlets")]
    )
    def test_a_layer_without_its_file(
        self,
        tmp_path: Path,
        geojson_network: tuple[str, ...],
        outlets: tuple[str, ...],
        layer: str,
        file: str,
    ) -> None:
        given = () if file == "--outlets" else outlets
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            *given,
            *utm(*at(500, 50)),
            layer,
            "t",
            says=(f"Invalid value for {layer}", file),
        )

    def test_out_must_be_geojson(self, tmp_path: Path, geojson_network: tuple[str, ...]) -> None:
        code, output = invoke(
            COMMAND,
            *geojson_network,
            *ATTRIBUTES,
            "--unit",
            "SE100",
            "--out",
            str(tmp_path / "n.vtk"),
        )
        assert code == 2, output
        assert "Invalid value for --out" in output, output
        assert not (tmp_path / "n.vtk").exists()


class TestFileRefusals:
    @pytest.mark.parametrize(
        "attributes",
        [
            ("--id-attribute", "nope_id", "--down-attribute", DOWN),
            ("--id-attribute", ID, "--down-attribute", "nope_down"),
            (*ATTRIBUTES, "--upstream-area-attribute", "nope_area", "--area-unit", "m2"),
            (*ATTRIBUTES, "--name-attribute", "nope_name"),
        ],
        ids=["id", "down", "area", "name"],
    )
    def test_a_missing_attribute_in_the_network(
        self, tmp_path: Path, network: tuple[str, ...], attributes: tuple[str, ...]
    ) -> None:
        missing = next(a for a in attributes if a.startswith("nope_"))
        refuse(
            tmp_path,
            *network,
            *attributes,
            "--unit",
            "SE100",
            says=("Invalid value for --network", missing),
        )

    def test_a_missing_attribute_in_the_outlets(
        self, tmp_path: Path, geojson_network: tuple[str, ...]
    ) -> None:
        path = outlets_geojson(tmp_path / "o.geojson", key="other_key")
        refuse(
            tmp_path,
            *geojson_network,
            *ATTRIBUTES,
            "--outlets",
            str(path),
            *utm(*at(500, 50)),
            says=("Invalid value for --outlets", ID),
        )

    @pytest.mark.parametrize("own", [None, ""], ids=["null", "empty"])
    def test_a_null_or_empty_own_id(self, tmp_path: Path, own: Any) -> None:
        units = [*NETWORK, (own, "SE100", 1 * KM2, "Nameless", square(-1, 0))]
        path = network_geojson(tmp_path / "bad.geojson", units)
        refuse(
            tmp_path,
            "--network",
            str(path),
            *ATTRIBUTES,
            "--unit",
            "SE100",
            says=("Invalid value for --network", ID),
        )

    def test_a_duplicate_own_id(self, tmp_path: Path, network_factory: Any) -> None:
        units = [*NETWORK, ("SE300", END, 1 * KM2, "Again", square(-1, 0))]
        refuse(
            tmp_path,
            *network_factory(units),
            *ATTRIBUTES,
            "--unit",
            "SE100",
            says=("Invalid value for --network", "SE300"),
        )

    def test_unit_names_no_unit(self, tmp_path: Path, network: tuple[str, ...]) -> None:
        refuse(
            tmp_path,
            *network,
            *ATTRIBUTES,
            "--unit",
            "SE999",
            says=("Invalid value for --unit", "SE999"),
        )

    def test_the_seed_in_no_unit(
        self, tmp_path: Path, network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        refuse(
            tmp_path,
            *network,
            *ATTRIBUTES,
            *outlets,
            *utm(*at(-500, 500)),
            says=("Invalid value for --seed", "no unit"),
        )

    def test_the_seed_in_the_notch_of_an_l_shaped_unit_is_in_no_unit(
        self, tmp_path: Path, network_factory: Any, outlets: tuple[str, ...]
    ) -> None:
        """The point is inside the unit's bounding box (an R-tree finds it) but
        not inside the unit."""
        ell = shapely.union_all([square(0, 0), square(0, 1), square(1, 1)])
        units: list[Unit] = [("SE100", END, None, None, ell)]
        refuse(
            tmp_path,
            *network_factory(units),
            *ATTRIBUTES,
            *outlets,
            *utm(*at(1500, 500)),
            says=("Invalid value for --seed", "no unit"),
        )

    def test_the_seed_on_the_border_of_two_units(
        self, tmp_path: Path, network: tuple[str, ...], outlets: tuple[str, ...]
    ) -> None:
        refuse(
            tmp_path,
            *network,
            *ATTRIBUTES,
            *outlets,
            *utm(*at(KM, 500)),
            says=("Invalid value for --seed", "SE100", "SE300"),
        )

    def test_the_seed_unit_on_a_cycle(self, tmp_path: Path, network_factory: Any) -> None:
        units = [
            ("SE100", "SE300", None, None, square(0, 0)),
            ("SE300", "SE100", None, None, square(1, 0)),
            ("SE500", "SE100", None, None, square(0, -1)),
        ]
        said = refuse(
            tmp_path,
            *network_factory(units),
            *ATTRIBUTES,
            "--unit",
            "SE100",
            says=("Invalid value for --network", "SE100", "SE300"),
        )
        assert "SE500" not in said, said

    def test_a_multipart_join(self, tmp_path: Path, network_factory: Any) -> None:
        apart = shapely.affinity.translate(square(1, 0), xoff=1.0)
        units = [("SE100", END, None, None, square(0, 0)), ("SE300", "SE100", None, None, apart)]
        refuse(
            tmp_path,
            *network_factory(units),
            *ATTRIBUTES,
            "--unit",
            "SE100",
            says=("Invalid value for --network", "km2"),
        )

    @pytest.mark.parametrize(("stated_km2", "accepted"), [(4.03, True), (4.1, False), (3.9, False)])
    def test_the_area_check(
        self, tmp_path: Path, network_factory: Any, stated_km2: float, accepted: bool
    ) -> None:
        units = [(NETWORK[0][0], END, stated_km2 * KM2, "Alpha", NETWORK[0][4]), *NETWORK[1:]]
        args = (*network_factory(units), *ATTRIBUTES, *STATED, "--unit", "SE100")
        if accepted:
            props, _ = run(tmp_path, *args)
            assert props["stated_upstream_area_km2"] == pytest.approx(stated_km2, rel=1e-9)
        else:
            refuse(tmp_path, *args, says=("Invalid value for --network", "lacks units"))


@pytest.fixture(params=["geojson", "gpkg"])
def network_factory(
    request: pytest.FixtureRequest, tmp_path: Path
) -> Callable[[Sequence[Unit]], tuple[str, ...]]:
    """Write `units` as a network file of each kind; its file options."""
    if request.param == "gpkg" and not HAS_RTREE:
        pytest.skip("this SQLite has no R-tree module")
    count = iter(range(1000))

    def written(units: Sequence[Unit]) -> tuple[str, ...]:
        k = next(count)
        if request.param == "gpkg":
            path = network_gpkg(tmp_path / f"net{k}.gpkg", units)
            return ("--network", str(path), "--network-layer", "units")
        return ("--network", str(network_geojson(tmp_path / f"net{k}.geojson", units)))

    return written


# ---------------------------------------------------------------------------
# The cross-check against the terrain route (section 3.5)
# ---------------------------------------------------------------------------

#: The bowl of `catchment_fixtures` (EPSG:25833, 100 m lattice): its lake
#: round the node (100, 200), its outlet channel down column 100. The home
#: unit holds the lake, its outlet on the lake's south side; the unit above
#: drains to it.
BOWL_CRS = "EPSG:25833"
LAKE_SOUTH = lat(100, 203.5)
BOWL_UNITS: list[Unit] = [
    ("H1", END, None, "Home", box(lat(90, 0)[0], LAKE_SOUTH[1], lat(110, 0)[0], lat(0, 150)[1])),
    ("U1", "H1", None, "Above", box(lat(90, 0)[0], lat(0, 150)[1], lat(110, 0)[0], lat(0, 100)[1])),
]
CHECK_KEYS = {
    "terrain_area_km2",
    "shared_km2",
    "network_only_km2",
    "terrain_only_km2",
    "units_mostly_outside",
}


def check_block(props: dict[str, Any]) -> dict[str, Any]:
    (block,) = [v for v in props.values() if isinstance(v, dict) and set(v) >= CHECK_KEYS]
    return block


@pytest.fixture
def bowl_network(tmp_path: Path) -> tuple[str, ...]:
    net = network_geojson(tmp_path / "bowl-units.geojson", BOWL_UNITS, crs=BOWL_CRS)
    outlet = collection(
        tmp_path / "bowl-outlets.geojson",
        [feature(Point(*LAKE_SOUTH), {ID: "H1"}), feature(Point(lat(100, 150)), {ID: "U1"})],
        BOWL_CRS,
    )
    lakes = lakes_geojson(tmp_path / "bowl-lakes.geojson", lake_box(), BOWL_CRS)
    seed = lat(100, 200)
    return (
        "--network",
        str(net),
        *ATTRIBUTES,
        "--outlets",
        str(outlet),
        "--lakes",
        str(lakes),
        "--seed",
        repr(seed[0]),
        repr(seed[1]),
        "--seed-crs",
        BOWL_CRS,
    )


@pytest.fixture
def dem_dir(tmp_path: Path) -> Path:
    """The bowl as four overlapping GeoTIFF tiles, as the `catchment` suite writes it."""
    write_tiles(tmp_path / "tiles", quadrants(tile_of(bowl()), row_cut=200, col_cut=100, overlap=1))
    return tmp_path / "tiles"


class TestCheckDem:
    def test_the_check_reports_the_four_areas_and_the_units_mostly_outside(
        self, tmp_path: Path, dem_dir: Path, bowl_network: tuple[str, ...]
    ) -> None:
        props, _ = run(tmp_path, *bowl_network, "--check-dem", str(dem_dir))
        assert props["units"] == ["H1", "U1"]
        block = check_block(props)
        assert block["terrain_area_km2"] > 0.0
        # The two partitions of each polygon, to 1e-6 relative at about 20 km2.
        assert block["shared_km2"] + block["network_only_km2"] == pytest.approx(
            props["area_km2"], rel=1e-6
        )
        assert block["shared_km2"] + block["terrain_only_km2"] == pytest.approx(
            block["terrain_area_km2"], rel=1e-6
        )

    def test_a_failed_terrain_route_still_writes_the_network(
        self, tmp_path: Path, bowl_network: tuple[str, ...]
    ) -> None:
        empty = tmp_path / "empty-cache"
        empty.mkdir()
        props, output = run(tmp_path, *bowl_network, "--check-dem", "glo30", "--cache", str(empty))
        assert props["units"] == ["H1", "U1"]
        assert "not in the cache" in output, output


# ---------------------------------------------------------------------------
# The MCP tool and the discovery text (sections 3.1 and 3.6)
# ---------------------------------------------------------------------------


class TestTheTool:
    def test_the_tool_requires_the_network_and_its_two_attributes(self) -> None:
        from tin_engine.cli import app
        from tin_engine.mcp.schema import derive

        tools = {d.tool_name: d for d in derive(app)}
        schema = tools["network_catchment"].tool().input_schema
        assert set(schema["required"]) == {"out", "network", "id_attribute", "down_attribute"}

    def test_the_server_instructions_name_the_tool_and_the_check(self) -> None:
        pytest.importorskip("mcp")
        from tin_engine.mcp import server

        text = server.INSTRUCTIONS
        assert "network_catchment" in text, text
        assert "check_dem" in text, text
        assert "downstream unit" in text, text

    def test_the_help_says_what_the_command_does(self) -> None:
        code, text = invoke(COMMAND, "--help")
        assert code == 0, text
        assert "downstream unit" in plain(text), text
