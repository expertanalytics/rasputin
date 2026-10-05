"""`io/station_set.py`: stations and reference polygons read back (29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "The station set" and "The
red suites", PR 3's `test_station_set.py`. Hand-written files only. Each
refusal is a `ValueError` (a Pydantic `ValidationError` is one) whose message
names what is wrong.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from pydantic import ValidationError
from shapely.geometry import MultiPolygon, Polygon, box

from nve_fixtures import CRS, collection, crs_member, feature, line, point, write


@pytest.fixture
def station_set() -> ModuleType:
    import tin_engine.io.station_set as module

    return module


def nve_station(
    number: str, x: float = 300_000.0, y: float = 6_600_000.0, **extra: Any
) -> dict[str, Any]:
    """A station as `fetch-stations` writes it."""
    props = {
        "station": number,
        "name": "Narsjø",
        "series": ["1001.0"],
        "nve_area_km2": 119.6,
        "hrd_start_daily": 1931,
        "watercourse": "002.DC",
        "river": "Glåma",
    }
    return point(x, y, **(props | extra))


class TestReadStations:
    def test_the_fields_and_the_crs(self, station_set: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "s.geojson", collection([nve_station("2.11.0", 1.5, 2.5)]))
        stations, crs = station_set.read_stations(path)
        assert crs == CRS
        (s,) = stations
        assert isinstance(s, station_set.Station)
        assert (s.station, s.name, s.x, s.y) == ("2.11.0", "Narsjø", 1.5, 2.5)
        assert s.series == ("1001.0",)
        assert s.nve_area_km2 == 119.6
        assert (s.watercourse, s.river) == ("002.DC", "Glåma")

    def test_file_order_is_kept(self, station_set: ModuleType, tmp_path: Path) -> None:
        numbers = ["313.10.0", "2.11.0", "19.79.0"]
        path = write(tmp_path / "s.geojson", collection([nve_station(n) for n in numbers]))
        stations, _ = station_set.read_stations(path)
        assert [s.station for s in stations] == numbers

    def test_the_files_own_crs_is_returned(self, station_set: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "s.geojson", collection([nve_station("2.11.0")], crs="EPSG:32633"))
        assert station_set.read_stations(path)[1] == "EPSG:32633"

    def test_a_users_own_points_file_reads(self, station_set: ModuleType, tmp_path: Path) -> None:
        """ "Any GeoJSON of points with a `station` property works": no
        watercourse, no river, no NVE area."""
        path = write(tmp_path / "mine.geojson", collection([point(10.0, 20.0, station="7.1.0")]))
        (s,), crs = station_set.read_stations(path)
        assert crs == CRS
        assert (s.station, s.x, s.y) == ("7.1.0", 10.0, 20.0)
        assert (s.watercourse, s.river, s.nve_area_km2) == (None, None, None)

    def test_a_station_is_frozen(self, station_set: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "s.geojson", collection([nve_station("2.11.0")]))
        (s,), _ = station_set.read_stations(path)
        with pytest.raises(ValidationError, match="frozen"):
            s.x = 0.0

    def test_a_file_without_crs_is_refused(self, station_set: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "s.geojson", collection([nve_station("2.11.0")], crs=None))
        with pytest.raises(ValueError, match=r"(?i)\bcrs\b"):
            station_set.read_stations(path)

    def test_a_duplicate_station_number_is_refused(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        doc = collection([nve_station("2.11.0"), nve_station("2.32.0"), nve_station("2.11.0", 5.0)])
        with pytest.raises(ValueError, match=r"2\.11\.0"):
            station_set.read_stations(write(tmp_path / "s.geojson", doc))

    def test_a_line_is_refused(self, station_set: ModuleType, tmp_path: Path) -> None:
        doc = collection([nve_station("2.11.0"), line([(0, 0), (1, 1)], station="2.32.0")])
        with pytest.raises(ValueError, match="LineString"):
            station_set.read_stations(write(tmp_path / "s.geojson", doc))

    @pytest.mark.parametrize("number", ["2.11", "2.11.0a", "2-11-0", "", " 2.11.0"])
    def test_a_bad_station_number_is_refused(
        self, station_set: ModuleType, tmp_path: Path, number: str
    ) -> None:
        """The pattern is `^\\d+\\.\\d+\\.\\d+$`: regine, main number, point number."""
        path = write(tmp_path / "s.geojson", collection([nve_station(number)]))
        with pytest.raises(ValueError, match="station"):
            station_set.read_stations(path)


class TestReadReferences:
    def test_polygons_and_multipolygons_by_station(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        two = MultiPolygon([box(0, 0, 10, 10), box(20, 0, 30, 10)])
        doc = collection(
            [
                feature(box(0, 0, 100, 100), station="2.11.0", reference_area_km2=0.01),
                feature(two, station="19.79.0", reference_area_km2=0.0002),
            ]
        )
        references, crs = station_set.read_references(write(tmp_path / "r.geojson", doc))
        assert crs == CRS
        assert set(references) == {"2.11.0", "19.79.0"}
        assert isinstance(references["2.11.0"], Polygon)
        assert references["2.11.0"].area == 10_000.0
        assert isinstance(references["19.79.0"], MultiPolygon)
        assert len(references["19.79.0"].geoms) == 2

    def test_a_file_without_crs_is_refused(self, station_set: ModuleType, tmp_path: Path) -> None:
        doc = collection([feature(box(0, 0, 1, 1), station="2.11.0")], crs=None)
        with pytest.raises(ValueError, match=r"(?i)\bcrs\b"):
            station_set.read_references(write(tmp_path / "r.geojson", doc))

    def test_a_duplicate_station_number_is_refused(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        doc = collection(
            [feature(box(0, 0, 1, 1), station="2.11.0"), feature(box(5, 5, 6, 6), station="2.11.0")]
        )
        with pytest.raises(ValueError, match=r"2\.11\.0"):
            station_set.read_references(write(tmp_path / "r.geojson", doc))

    def test_a_point_is_refused(self, station_set: ModuleType, tmp_path: Path) -> None:
        doc = collection([point(0.0, 0.0, station="2.11.0")])
        with pytest.raises(ValueError, match="Point"):
            station_set.read_references(write(tmp_path / "r.geojson", doc))


class TestAMalformedFile:
    """Code review round 1: a user's malformed station or reference file is
    refused with a `ValueError` naming what is wrong, never a `KeyError` or
    `TypeError` escaping from the reader."""

    def test_a_point_without_coordinates_names_its_station(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        broken = nve_station("2.32.0")
        broken["geometry"]["coordinates"] = None
        doc = collection([nve_station("2.11.0"), broken])
        with pytest.raises(ValueError, match=r"2\.32\.0"):
            station_set.read_stations(write(tmp_path / "s.geojson", doc))

    def test_a_stations_file_without_features(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        doc = collection([nve_station("2.11.0")])
        del doc["features"]
        with pytest.raises(ValueError, match="features"):
            station_set.read_stations(write(tmp_path / "s.geojson", doc))

    def test_a_crs_member_without_a_name(self, station_set: ModuleType, tmp_path: Path) -> None:
        doc = collection([nve_station("2.11.0")])
        doc["crs"] = {"type": "name"}
        with pytest.raises(ValueError, match=r"(?i)\bcrs\b"):
            station_set.read_stations(write(tmp_path / "s.geojson", doc))

    def test_a_reference_without_a_station(self, station_set: ModuleType, tmp_path: Path) -> None:
        doc = collection(
            [feature(box(0, 0, 1, 1), station="2.11.0"), feature(box(5, 5, 6, 6), area=1.0)]
        )
        with pytest.raises(ValueError, match="station"):
            station_set.read_references(write(tmp_path / "r.geojson", doc))

    def test_a_reference_with_null_properties(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        bare = feature(box(5, 5, 6, 6))
        bare["properties"] = None
        doc = collection([feature(box(0, 0, 1, 1), station="2.11.0"), bare])
        with pytest.raises(ValueError, match=r"station|properties"):
            station_set.read_references(write(tmp_path / "r.geojson", doc))


class TestAReferenceWithNoArea:
    """Change (b) of "PR 4's green step": `agreement` divides by NVE's polygon's
    area and perimeter, so a reference with no area would stop a batch at its
    station (a zero-area ring raises `ZeroDivisionError`, an empty polygon a
    `ValueError` from its NaN bounds). `read_references` refuses either when
    the file is read, naming the station. The bad feature comes second, after
    a good one, so the message names the right station."""

    @pytest.mark.parametrize(
        "geometry",
        [
            {"type": "Polygon", "coordinates": [[[0, 0], [10, 0], [20, 0], [0, 0]]]},
            {"type": "Polygon", "coordinates": []},
            {"type": "MultiPolygon", "coordinates": []},
        ],
        ids=["zero_area_ring", "empty_polygon", "empty_multipolygon"],
    )
    def test_it_is_refused_naming_the_station(
        self, station_set: ModuleType, tmp_path: Path, geometry: dict[str, Any]
    ) -> None:
        doc = collection(
            [feature(box(0, 0, 100, 100), station="2.11.0"), feature(geometry, station="2.32.0")]
        )
        with pytest.raises(ValueError, match=r"2\.32\.0"):
            station_set.read_references(write(tmp_path / "r.geojson", doc))


# ---------------------------------------------------------------------------
# PR 4, lake gauges: `read_nve_lakes`, `read_lakes` until audit PR C
# ("The station set", "Lake gauges")
# ---------------------------------------------------------------------------


def nve_lake(geometry: Any, objectid: int | None = 4_100_001, **extra: Any) -> dict[str, Any]:
    """A lake as `fetch-stations` writes it: `objectid`, `vatnlnr`, `navn`,
    `areal_km2`; `objectid=None` leaves the property out."""
    props: dict[str, Any] = {"vatnlnr": 495, "navn": "Narsjøen", "areal_km2": 0.01} | extra
    if objectid is not None:
        props["objectid"] = objectid
    return feature(geometry, **props)


class TestReadLakes:
    """`read_nve_lakes(path) -> (tuple[Lake, ...], crs)`. Before increment 29's
    PR 4, neither `read_lakes` nor `Lake` existed; audit PR C renamed the
    reader by what it reads (`docs/increments/python-audit.md`, section 10)."""

    def test_the_name_says_what_it_reads(self, station_set: ModuleType) -> None:
        """Audit PR C, red test 5: no alias keeps the old name."""
        assert callable(getattr(station_set, "read_nve_lakes", None))
        assert not hasattr(station_set, "read_lakes")
        assert "read_nve_lakes" in station_set.__all__
        assert "read_lakes" not in station_set.__all__

    def test_a_polygon_gives_one_lake_with_its_number_and_name(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        doc = collection([nve_lake(box(0, 0, 100, 100))])
        lakes, crs = station_set.read_nve_lakes(write(tmp_path / "l.geojson", doc))
        assert crs == CRS
        (lake,) = lakes
        assert isinstance(lake, station_set.Lake)
        assert (lake.number, lake.name) == (495, "Narsjøen")
        assert isinstance(lake.polygon, Polygon) and lake.polygon.area == 10_000.0

    def test_a_multipolygon_gives_one_lake_per_part_with_the_same_number_and_name(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        two = MultiPolygon([box(0, 0, 10, 10), box(20, 0, 30, 10)])
        doc = collection([nve_lake(two, vatnlnr=12, navn="Tvillingvatna")])
        lakes, _ = station_set.read_nve_lakes(write(tmp_path / "l.geojson", doc))
        assert len(lakes) == 2
        assert {(lk.number, lk.name) for lk in lakes} == {(12, "Tvillingvatna")}
        assert all(isinstance(lk.polygon, Polygon) for lk in lakes)
        assert sorted(lk.polygon.bounds for lk in lakes) == [(0, 0, 10, 10), (20, 0, 30, 10)]

    def test_file_order_is_kept(self, station_set: ModuleType, tmp_path: Path) -> None:
        doc = collection(
            [
                nve_lake(box(0, 0, 1, 1), objectid=3, vatnlnr=30),
                nve_lake(box(5, 5, 6, 6), objectid=1, vatnlnr=10),
                nve_lake(box(9, 9, 10, 10), objectid=2, vatnlnr=20),
            ]
        )
        lakes, _ = station_set.read_nve_lakes(write(tmp_path / "l.geojson", doc))
        assert [lk.number for lk in lakes] == [30, 10, 20]

    @pytest.mark.parametrize("vatnlnr", [0, None], ids=["zero", "null"])
    def test_lake_number_0_or_null_is_no_number(
        self, station_set: ModuleType, tmp_path: Path, vatnlnr: int | None
    ) -> None:
        doc = collection([nve_lake(box(0, 0, 1, 1), vatnlnr=vatnlnr, navn=None)])
        (lake,), _ = station_set.read_nve_lakes(write(tmp_path / "l.geojson", doc))
        assert lake.number is None and lake.name is None

    def test_a_file_without_crs_is_refused(self, station_set: ModuleType, tmp_path: Path) -> None:
        doc = collection([nve_lake(box(0, 0, 1, 1))], crs=None)
        with pytest.raises(ValueError, match=r"(?i)\bcrs\b"):
            station_set.read_nve_lakes(write(tmp_path / "l.geojson", doc))

    @pytest.mark.parametrize(
        "geometry",
        [
            {"type": "LineString", "coordinates": [[0, 0], [10, 0]]},
            {"type": "Point", "coordinates": [0, 0]},
            {"type": "Polygon", "coordinates": []},
            {"type": "MultiPolygon", "coordinates": []},
            {"type": "Polygon", "coordinates": [[[0, 0], [10, 0], [20, 0], [0, 0]]]},
        ],
        ids=["line_string", "point", "empty_polygon", "empty_multipolygon", "zero_area_ring"],
    )
    def test_a_bad_geometry_is_refused_naming_its_objectid(
        self, station_set: ModuleType, tmp_path: Path, geometry: dict[str, Any]
    ) -> None:
        """After a good lake, so the message names the right one."""
        doc = collection([nve_lake(box(0, 0, 1, 1)), nve_lake(geometry, objectid=7_654_321)])
        with pytest.raises(ValueError, match="7654321"):
            station_set.read_nve_lakes(write(tmp_path / "l.geojson", doc))

    def test_a_bad_geometry_without_objectid_is_refused_naming_its_index(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        """The fourth feature (index 3) has no `objectid`."""
        good = [nve_lake(box(k, 0, k + 1, 1), objectid=900 + k) for k in range(3)]
        bad = nve_lake({"type": "LineString", "coordinates": [[0, 0], [10, 0]]}, objectid=None)
        doc = collection([*good, bad])
        with pytest.raises(ValueError, match=r"(?i)(feature|index)\D{0,4}3\b"):
            station_set.read_nve_lakes(write(tmp_path / "l.geojson", doc))


# ---------------------------------------------------------------------------
# Audit PR C: the one GeoJSON reading path, with no default CRS
# ---------------------------------------------------------------------------


class TestOneGeojsonRule:
    """`docs/increments/python-audit.md`, section 10, red test 3: `features_of`
    reads through `io.geojson.read_collection(..., default_crs=None)`, and
    every refusal, the decoder's among them, starts with the file's name."""

    def test_a_single_feature_file_is_one_station(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        """Before: refused as not a FeatureCollection."""
        doc = nve_station("2.11.0", 1.5, 2.5) | {"crs": crs_member()}
        (s,), crs = station_set.read_stations(write(tmp_path / "s.geojson", doc))
        assert crs == CRS
        assert (s.station, s.x, s.y) == ("2.11.0", 1.5, 2.5)

    def test_a_feature_that_is_not_an_object_is_refused_naming_the_file(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        """Before: an `AttributeError` from `_geometry_type`, which
        `station-catchments` does not catch."""
        path = write(tmp_path / "s.geojson", collection([]) | {"features": [7]})
        with pytest.raises(ValueError) as info:
            station_set.read_stations(path)
        assert str(info.value) == (
            "s.geojson: no features list; the file is not a FeatureCollection"
        )

    def test_a_null_member_is_refused_as_null(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        """Before: refused as `no crs member`."""
        path = write(tmp_path / "s.geojson", collection([nve_station("2.11.0")]) | {"crs": None})
        with pytest.raises(ValueError) as info:
            station_set.read_stations(path)
        assert str(info.value) == "s.geojson: the crs member is null; the file must name its CRS"

    def test_an_unreadable_name_is_refused_naming_the_file(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        """Before: `cannot read the CRS ...`, without the file."""
        path = write(tmp_path / "s.geojson", collection([nve_station("2.11.0")], crs="EPSG:999999"))
        with pytest.raises(ValueError) as info:
            station_set.read_stations(path)
        assert str(info.value).startswith("s.geojson: cannot read the CRS 'EPSG:999999'")

    def test_a_file_that_is_not_json_is_refused_naming_the_file(
        self, station_set: ModuleType, tmp_path: Path
    ) -> None:
        """Before: the decoder's bare text."""
        path = tmp_path / "s.geojson"
        path.write_text('{"type": "FeatureCollection", "features": [', encoding="utf-8")
        with pytest.raises(ValueError) as info:
            station_set.read_stations(path)
        assert str(info.value).startswith("s.geojson: ")


# The four readers that go through `features_of`, by module and name.
FEATURES_OF_READERS = [
    ("tin_engine.io.station_set", "read_stations"),
    ("tin_engine.io.station_set", "read_references"),
    ("tin_engine.io.station_set", "read_nve_lakes"),
    ("tin_engine.io.rivers", "read_segments"),
]


class TestAFileWithNoGeometry:
    """`docs/increments/python-audit.md`, section 10's wording table, the row
    "a `crs` member and either a `Feature` without `geometry` or an object
    with neither `type` nor `features`". The shape rule reads either as one
    feature; the reader then refused it as `None is a None, not a Point` (or
    the Polygon, lake and LineString forms). It is refused in plain words
    naming the file and the feature: `has no geometry`.
    The `"geometry": 7` case is a known gap, kept out."""

    @pytest.mark.parametrize(
        "doc",
        [
            pytest.param({"type": "Feature", "crs": crs_member(), "properties": {}}, id="feature"),
            pytest.param({"crs": crs_member(), "foo": 1}, id="typeless"),
        ],
    )
    @pytest.mark.parametrize(
        ("module", "reader"), FEATURES_OF_READERS, ids=[r for _, r in FEATURES_OF_READERS]
    )
    def test_it_is_refused_saying_the_feature_has_no_geometry(
        self, module: str, reader: str, doc: dict[str, Any], tmp_path: Path
    ) -> None:
        read = getattr(importlib.import_module(module), reader)
        with pytest.raises(ValueError) as info:
            read(write(tmp_path / "nogeom.geojson", doc))
        message = str(info.value)
        assert "nogeom.geojson" in message
        assert "feature 0" in message
        assert "has no geometry" in message
