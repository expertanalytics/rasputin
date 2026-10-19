"""`tin_engine.gml_import.gml_to_geopackage`: a GML 3.2 document into one GeoPackage table.

`docs/increments/73-gml-vector-source.md`, sections 4.2 and 9, tests 5 and 6:
fixture A converts into a GeoPackage that the existing readers
(`layer_info`, `layer_summary`, `distinct_values`, `feature_input.read_source`)
read unchanged; the converter's rules on columns, their types, the one feature
type and the one CRS are each refused naming what they name.

The converter is imported inside each test, so each fails on its own while
`tin_engine.gml_import` does not exist. A refusal is asserted as a
`ValueError` (a `GmlError` is one) with its words.
"""

from __future__ import annotations

import importlib
import io
import sqlite3
from collections.abc import Iterator
from contextlib import closing
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from gml32_fixtures import (
    AR50_FIGURES,
    AR50_SAMPLE,
    AR50_VALUES,
    URI_4326,
    URI_25832,
    document,
    feature,
    polygon,
    values,
)
from tin_engine.feature_input import read_source
from tin_engine.io.geopackage import distinct_values, layer_info, layer_summary
from tin_engine.io.repository import create_sqlite, open_geopackage

#: A box in UTM 32 holding every fixture-A feature (eastings near 621 km,
#: northings near 6,769 km) and every hand-made square near the origin.
EVERYTHING = (-1.0e7, -1.0e7, 1.0e7, 1.0e8)


@pytest.fixture(scope="module")
def gml_import() -> ModuleType:
    return importlib.import_module("tin_engine.gml_import")


def converted(gml_import: ModuleType, stream: Any, path: Path, **kw: Any) -> Any:
    """``stream`` converted into a new GeoPackage at ``path``, the connection closed."""
    with closing(create_sqlite(path)) as conn:
        return gml_import.gml_to_geopackage(stream, conn, **kw)


@pytest.fixture
def opened(tmp_path: Path) -> Iterator[Any]:
    """Open ``tmp_path / "out.gpkg"`` read-only, closing it after the test."""
    conns: list[sqlite3.Connection] = []

    def open_() -> sqlite3.Connection:
        conns.append(open_geopackage(tmp_path / "out.gpkg"))
        return conns[-1]

    yield open_
    for c in conns:
        c.close()


def refused(gml_import: ModuleType, tmp_path: Path, stream: Any, *words: str, **kw: Any) -> str:
    with pytest.raises(ValueError) as caught:
        converted(gml_import, stream, tmp_path / "out.gpkg", **kw)
    message = str(caught.value)
    for word in words:
        assert word in message, f"{word!r} not in {message!r}"
    return message


# ---------------------------------------------------------------- test 5


class TestRoundTrip:
    @pytest.fixture
    def done(self, gml_import: ModuleType, tmp_path: Path) -> Any:
        with AR50_SAMPLE.open("rb") as stream:
            return converted(gml_import, stream, tmp_path / "out.gpkg")

    def test_the_returned_import(self, done: Any) -> None:
        assert done.layer == "ArealressursFlate"
        assert done.crs == URI_25832  # as stored: the srsName text
        assert (done.features, done.without_surface) == (6, 0)
        assert done.columns == ("gml_id", *AR50_VALUES)

    def test_layer_info(self, done: Any, opened: Any) -> None:
        info = layer_info(opened(), None)
        assert (info.table, info.crs, info.srs_id) == ("ArealressursFlate", "EPSG:25832", 25832)
        assert info.rtree is not None

    def test_layer_summary(self, done: Any, opened: Any) -> None:
        conn = opened()
        summary = layer_summary(conn, layer_info(conn, None))
        assert (summary.features, summary.geometry) == (6, "MULTIPOLYGON")
        assert summary.attributes == ("gml_id", *AR50_VALUES)
        assert summary.extent is not None

    def test_feature_input_reads_it(self, done: Any, tmp_path: Path) -> None:
        read = read_source(tmp_path / "out.gpkg", None, "arealtype", lambda _: EVERYTHING)
        assert (read.crs, read.layer, read.scanned) == ("EPSG:25832", "ArealressursFlate", False)
        assert [fid for fid, _, _ in read.rows] == [1, 2, 3, 4, 5, 6]
        # Areas to 0.1 m2 as measured; the largest is 353,418.5 m2.
        assert [round(g.area, 1) for _, g, _ in read.rows] == [f.area_m2 for f in AR50_FIGURES]
        assert [v for _, _, v in read.rows] == [60, 60, 60, 60, 60, 30]

    def test_holes_survive(self, done: Any, tmp_path: Path) -> None:
        read = read_source(tmp_path / "out.gpkg", None, None, lambda _: EVERYTHING)
        holes = [sum(len(p.interiors) for p in g.geoms) for _, g, _ in read.rows]
        assert holes == [f.holes for f in AR50_FIGURES]

    def test_distinct_values(self, done: Any, opened: Any) -> None:
        conn = opened()
        assert distinct_values(conn, layer_info(conn, None), "arealtype", 100) == (("30", "60"), 2)

    def test_gml_ids_in_document_order(self, done: Any, opened: Any) -> None:
        rows = opened().execute('SELECT gml_id FROM "ArealressursFlate" ORDER BY fid').fetchall()
        assert [r for (r,) in rows] == [f.gml_id for f in AR50_FIGURES]


# ---------------------------------------------------------------- test 6


def column_types(conn: sqlite3.Connection, table: str = "F") -> dict[str, str]:
    return {row[1]: row[2] for row in conn.execute(f'PRAGMA table_info("{table}")')}


class TestColumns:
    def test_ar50_types(self, gml_import: ModuleType, tmp_path: Path, opened: Any) -> None:
        with AR50_SAMPLE.open("rb") as stream:
            converted(gml_import, stream, tmp_path / "out.gpkg")
        types = column_types(opened(), "ArealressursFlate")
        assert types["arealtype"] == "INTEGER"
        assert types["informasjon"] == "TEXT"
        assert types["gml_id"] == "TEXT"
        assert list(types)[:3] == ["fid", "geom", "gml_id"]

    @pytest.mark.parametrize(
        ("text", "kind"),
        [("60", "INTEGER"), ("-5", "INTEGER"), ("+7", "INTEGER"), ("007", "INTEGER"),
         ("1.5", "TEXT"), ("1e3", "TEXT"), ("AR50", "TEXT"), ("6 0", "TEXT")],
    )  # fmt: skip
    def test_a_column_is_integer_only_for_a_decimal_integer(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any, text: str, kind: str
    ) -> None:
        converted(
            gml_import,
            document(feature("f1", polygon(), props=values(a=text))),
            tmp_path / "out.gpkg",
        )
        assert column_types(opened())["a"] == kind

    def test_a_later_non_integer_value_is_kept_as_text(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any
    ) -> None:
        stream = document(
            feature("f1", polygon(), props=values(a="60")),
            feature("f2", polygon(), props=values(a="x7")),
        )
        converted(gml_import, stream, tmp_path / "out.gpkg")
        rows = opened().execute('SELECT a, typeof(a) FROM "F" ORDER BY fid').fetchall()
        assert rows == [(60, "integer"), ("x7", "text")]

    def test_a_missing_and_a_nil_value_are_null(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any
    ) -> None:
        stream = document(
            feature("f1", polygon(), props=values(a="60", b="x")),
            feature("f2", polygon(), props=values(a="61")),
            feature("f3", polygon(), props=values(a=None, b="y")),
        )
        done = converted(gml_import, stream, tmp_path / "out.gpkg")
        rows = opened().execute('SELECT gml_id, a, b FROM "F" ORDER BY fid').fetchall()
        assert rows == [("f1", 60, "x"), ("f2", 61, None), ("f3", None, "y")]
        assert done.columns == ("gml_id", "a", "b")

    def test_a_later_value_the_first_feature_lacked(
        self, gml_import: ModuleType, tmp_path: Path
    ) -> None:
        stream = document(
            feature("f1", polygon(), props=values(a="60")),
            feature("f.extra", polygon(), props=values(a="60", extra="1")),
        )
        refused(gml_import, tmp_path, stream, "f.extra", "extra")

    def test_a_second_feature_type(self, gml_import: ModuleType, tmp_path: Path) -> None:
        stream = document(
            feature("f1", polygon()), feature("f.other", polygon(), type_name="Other")
        )
        refused(gml_import, tmp_path, stream, "f.other", "Other")

    @pytest.mark.parametrize("name", ["gml_id", "fid", "geom"])
    def test_a_value_with_a_reserved_name(
        self, gml_import: ModuleType, tmp_path: Path, name: str
    ) -> None:
        stream = document(feature("f.reserved", polygon(), props=values(**{name: "1"})))
        refused(gml_import, tmp_path, stream, "f.reserved", name)

    def test_many_chunks_keep_every_feature_in_order(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any
    ) -> None:
        n = 2 * gml_import.CHUNK + 1
        stream = document(*(feature(f"f{k}", polygon(), props=values(a=str(k))) for k in range(n)))
        done = converted(gml_import, stream, tmp_path / "out.gpkg")
        assert done.features == n
        rows = opened().execute('SELECT fid, gml_id, a FROM "F" ORDER BY fid').fetchall()
        assert rows == [(k + 1, f"f{k}", k) for k in range(n)]


class TestCrs:
    def test_two_srs_names(self, gml_import: ModuleType, tmp_path: Path) -> None:
        other = "http://www.opengis.net/def/crs/EPSG/0/25833"
        stream = document(feature("f1", polygon()), feature("f.second", polygon(srs_name=other)))
        refused(gml_import, tmp_path, stream, "f.second", URI_25832, other)

    def test_no_srs_name_and_no_crs(self, gml_import: ModuleType, tmp_path: Path) -> None:
        stream = document(feature("f.bare", polygon(srs_name=None)))
        refused(gml_import, tmp_path, stream, "feature f.bare names no srsName; give --crs")

    def test_a_crs_that_differs_from_the_srs_name(
        self, gml_import: ModuleType, tmp_path: Path
    ) -> None:
        stream = document(feature("f.both", polygon()))
        refused(gml_import, tmp_path, stream, "EPSG:25833", crs="EPSG:25833")

    def test_a_crs_the_srs_name_equals_is_taken(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any
    ) -> None:
        converted(
            gml_import, document(feature("f1", polygon())), tmp_path / "out.gpkg", crs="EPSG:25832"
        )
        assert layer_info(opened(), None).srs_id == 25832

    def test_a_latitude_first_crs(self, gml_import: ModuleType, tmp_path: Path) -> None:
        stream = document(feature("f.lat", polygon(srs_name=URI_4326)))
        refused(gml_import, tmp_path, stream, URI_4326)

    def test_a_crs_that_does_not_parse(self, gml_import: ModuleType, tmp_path: Path) -> None:
        stream = document(feature("f.junk", polygon(srs_name="urn:not:a:crs")))
        refused(gml_import, tmp_path, stream, "urn:not:a:crs")

    def test_the_given_crs_when_no_srs_name(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any
    ) -> None:
        stream = document(feature("f1", polygon(srs_name=None)))
        done = converted(gml_import, stream, tmp_path / "out.gpkg", crs="EPSG:25833")
        assert done.crs == "EPSG:25833"
        info = layer_info(opened(), None)
        assert (info.srs_id, info.crs) == (25833, "EPSG:25833")

    def test_a_later_geometry_without_srs_name_takes_the_fixed_one(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any
    ) -> None:
        stream = document(feature("f1", polygon()), feature("f2", polygon(srs_name=None)))
        assert converted(gml_import, stream, tmp_path / "out.gpkg").features == 2
        assert layer_info(opened(), None).srs_id == 25832


class TestWithoutSurface:
    def test_a_feature_with_no_geometry_is_skipped_and_counted(
        self, gml_import: ModuleType, tmp_path: Path, opened: Any
    ) -> None:
        stream = document(
            feature("f0"), feature("f1", polygon()), feature("f2"), feature("f3", polygon())
        )
        done = converted(gml_import, stream, tmp_path / "out.gpkg")
        assert (done.features, done.without_surface) == (2, 2)
        rows = opened().execute('SELECT fid, gml_id FROM "F" ORDER BY fid').fetchall()
        assert rows == [(1, "f1"), (2, "f3")]

    def test_a_document_with_no_surface_at_all(
        self, gml_import: ModuleType, tmp_path: Path
    ) -> None:
        refused(gml_import, tmp_path, document(feature("f0"), feature("f1")), "surface")

    def test_the_converter_opens_nothing(self, gml_import: ModuleType) -> None:
        """A pure function of a stream and a connection: it takes no path."""
        import inspect

        assert list(inspect.signature(gml_import.gml_to_geopackage).parameters) == [
            "stream",
            "conn",
            "crs",
        ]


def test_a_decoder_refusal_passes_through(gml_import: ModuleType, tmp_path: Path) -> None:
    stream = io.BytesIO(document(feature("f.bad", polygon("0 0 10 0 0 0"))).getvalue())
    refused(gml_import, tmp_path, stream, "f.bad")
