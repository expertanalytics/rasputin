"""`io/station_set.py`: stations and reference polygons read back (29, PR 3).

`docs/increments/29-nve-reference-catchments.md`, "The station set" and "The
red suites", PR 3's `test_station_set.py`. Hand-written files only. Each
refusal is a `ValueError` (a Pydantic `ValidationError` is one) whose message
names what is wrong.

HOW THIS FILE GOES RED: there is no `tin_engine/io/station_set.py`
(`ModuleNotFoundError` at collection of each test's import).
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from pydantic import ValidationError
from shapely.geometry import MultiPolygon, Polygon, box

from nve_fixtures import CRS, collection, feature, line, point, write


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
