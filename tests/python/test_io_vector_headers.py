"""What the listing reads from a vector file's header: GeoPackage and GeoJSON.

`docs/increments/57-vector-sources.md`, section 4.3, test T3, on fixture V
(`vector_fixtures`): `geopackage.layer_summary` and `distinct_values` over an
open connection, `geojson.collection_summary` over a parsed document. Each
reads only through the existing modules, which open nothing themselves.

The extent is in the file's CRS (metres here), compared exactly: the
R-tree's and `gpkg_contents`' figures are stored doubles, and the GeoJSON
bounds are the fixture's own coordinates.
"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterator
from contextlib import closing
from pathlib import Path

import pytest

from catchment_fixtures import lake_box
from dataset_fixtures import lake_file
from gpkg_fixtures import needs_rtree
from tin_engine.io import geojson, geopackage
from tin_engine.io.repository import open_geopackage
from vector_fixtures import (
    CODE,
    LAND_BOX,
    LAND_CODES,
    LANDCOVER,
    OTHER,
    landcover_gpkg,
    set_contents_bounds,
)

#: Stated bounds for the layer with no R-tree; deliberately not the rows' own,
#: so the test tells which source the extent came from.
STATED = (499_000.0, 6_559_000.0, 521_000.0, 6_601_000.0)


def opened(path: Path) -> Iterator[sqlite3.Connection]:
    with closing(open_geopackage(path)) as conn:
        yield conn


@pytest.fixture
def cover(tmp_path: Path) -> Iterator[sqlite3.Connection]:
    yield from opened(landcover_gpkg(tmp_path / "landcover.gpkg"))


def summary(conn: sqlite3.Connection) -> geopackage.LayerSummary:
    return geopackage.layer_summary(conn, geopackage.layer_info(conn, LANDCOVER))


class TestT3GeoPackage:
    @needs_rtree
    def test_the_layer_summary_from_the_rtree(self, cover: sqlite3.Connection) -> None:
        got = summary(cover)
        assert got.geometry == "GEOMETRY"
        assert got.features == 2
        assert got.attributes == (CODE,)  # less OBJECTID and Shape
        assert got.extent == LAND_BOX

    def test_without_an_rtree_the_contents_bounds(self, tmp_path: Path) -> None:
        path = landcover_gpkg(tmp_path / "landcover.gpkg", rtree=False)
        set_contents_bounds(path, LANDCOVER, STATED)
        with closing(open_geopackage(path)) as conn:
            assert summary(conn).extent == STATED

    def test_with_neither_no_extent(self, tmp_path: Path) -> None:
        path = landcover_gpkg(tmp_path / "landcover.gpkg", rtree=False)
        with closing(open_geopackage(path)) as conn:
            got = summary(conn)
        assert got.extent is None
        assert got.features == 2

    def test_partly_set_contents_bounds_are_no_extent(self, tmp_path: Path) -> None:
        """All four bounds are needed: one unset is no box."""
        path = landcover_gpkg(tmp_path / "landcover.gpkg", rtree=False)
        set_contents_bounds(path, LANDCOVER, STATED)
        with closing(sqlite3.connect(path)) as con:
            con.execute("UPDATE gpkg_contents SET max_y = NULL WHERE table_name = ?", (LANDCOVER,))
            con.commit()
        with closing(open_geopackage(path)) as conn:
            assert summary(conn).extent is None

    def test_the_distinct_values_and_their_count(self, cover: sqlite3.Connection) -> None:
        layer = geopackage.layer_info(cover, LANDCOVER)
        assert geopackage.distinct_values(cover, layer, CODE, 100) == (LAND_CODES, 2)

    def test_a_limit_cuts_the_values_not_the_count(self, cover: sqlite3.Connection) -> None:
        layer = geopackage.layer_info(cover, LANDCOVER)
        assert geopackage.distinct_values(cover, layer, CODE, 1) == (LAND_CODES[:1], 2)

    def test_two_features_tables_and_no_layer_is_the_readers_refusal(
        self, cover: sqlite3.Connection
    ) -> None:
        with pytest.raises(geopackage.GeoPackageError) as caught:
            geopackage.layer_info(cover, None)
        assert "name one layer of the features tables" in str(caught.value)
        assert LANDCOVER in str(caught.value) and OTHER in str(caught.value)


def lake_doc(tmp_path: Path) -> dict[str, object]:
    loaded: dict[str, object] = json.loads(lake_file(tmp_path / "lakes.geojson").read_text())
    return loaded


class TestT3GeoJson:
    def test_the_lake_file(self, tmp_path: Path) -> None:
        got = geojson.collection_summary(lake_doc(tmp_path))
        assert set(got.geometry) == {"Polygon"}
        assert got.features == 1
        assert got.attributes == ("name",)
        assert got.extent == lake_box().bounds
        assert got.crs_stated is True

    def test_no_crs_member(self, tmp_path: Path) -> None:
        doc = lake_doc(tmp_path)
        del doc["crs"]
        assert geojson.collection_summary(doc).crs_stated is False

    def test_types_and_attributes_in_first_seen_order(self) -> None:
        doc = {
            "type": "FeatureCollection",
            "features": [
                {"type": "Feature", "properties": {"objectid": 1, "elvid": "a"},
                 "geometry": {"type": "LineString", "coordinates": [[0, 0], [1, 1]]}},
                {"type": "Feature", "properties": {"zeta": 2, "objectid": 2},
                 "geometry": {"type": "Point", "coordinates": [5, -3]}},
            ],
        }  # fmt: skip
        got = geojson.collection_summary(doc)
        assert set(got.geometry) == {"LineString", "Point"}
        assert got.features == 2
        assert got.attributes == ("objectid", "elvid", "zeta")
        assert got.extent == (0.0, -3.0, 5.0, 1.0)
        assert got.crs_stated is False
