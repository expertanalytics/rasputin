"""Feature layers read and written through the plugin interface.

`docs/increments/116-feature-table.md`, *3. Tests*: P1 to P5, plus one test
for each type rule and refusal its design names. Every fixture is built here
and its figures are the ones it is written with. P6 is the existing lists:
`test_plugin_api.py` (each name the core's own object, with its signature),
`test_cli_help.py` T6 (each name on `docs/plugins.md`) and `test_layering.py`.
"""

from __future__ import annotations

import json
import re
import sqlite3
from collections.abc import Sequence
from contextlib import closing
from pathlib import Path
from typing import Any

import pytest
from shapely.geometry import MultiPolygon, Polygon, box

from tin_engine import cli, sources
from tin_engine.crs import srs_entry
from tin_engine.feature_input import FeatureError
from tin_engine.io import geopackage_write
from tin_engine.io.feature_table import FeatureTable, read_feature_table
from tin_engine.sources import RemoteSource, features_file

UTM33 = "EPSG:25833"
#: Off every round number, UTM-shaped.
X0, Y0 = 500_000.3, 6_600_000.7
TABLE = "parcels"
COLUMNS = (
    ("code", "esriFieldTypeString"),
    ("n", "esriFieldTypeInteger"),
    ("x", "esriFieldTypeDouble"),
)


def square(x: float, y: float, side: float = 10.0) -> MultiPolygon:
    return MultiPolygon([box(x, y, x + side, y + side)])


#: Rows 7 and 9 sit 100 m apart; row 12 sits 10 km east. A 10 m square's
#: R-tree box is widened by its width plus height, 20 m, in the query.
ROWS: tuple[tuple[int, MultiPolygon, tuple[Any, ...]], ...] = (
    (7, square(X0, Y0), ("a7", 3, 0.125)),
    (9, square(X0 + 100.0, Y0), ("b9", -4, 2.5e-7)),
    (12, square(X0 + 10_000.0, Y0), ("c12", 0, 1234567.890625)),
)


def write_layer(
    path: Path,
    columns: Sequence[tuple[str, str]] = COLUMNS,
    rows: Sequence[tuple[int, MultiPolygon, tuple[Any, ...]]] = ROWS,
    extra_sql: Sequence[str] = (),
) -> Path:
    """A GeoPackage in EPSG:25833 with one features table `TABLE`."""
    srs = srs_entry(UTM33)
    with closing(sqlite3.connect(path)) as conn:
        geopackage_write.create_geopackage(conn, srs)
        geopackage_write.create_features_table(
            conn, TABLE, srs.srs_id, "MULTIPOLYGON", "fid", columns
        )
        geopackage_write.insert_features(conn, TABLE, srs.srs_id, rows)
        for sql in extra_sql:
            conn.execute(sql)
        conn.commit()
    return path


def write_geojson(
    path: Path, properties: Sequence[dict[str, Any]], ids: Sequence[int] = ()
) -> Path:
    features = []
    for i, props in enumerate(properties):
        ring = [[10.0 + i, 60.0], [10.1 + i, 60.0], [10.1 + i, 60.1], [10.0 + i, 60.0]]
        feature: dict[str, Any] = {
            "type": "Feature",
            "properties": props,
            "geometry": {"type": "Polygon", "coordinates": [ring]},
        }
        if ids:
            feature["id"] = ids[i]
        features.append(feature)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}), "utf-8")
    return path


# ---------------------------------------------------------------- P1


class TestAGeoPackageLayerReadsWhole:
    """P1: columns and types, fids, values exactly, geometries exactly, CRS text."""

    def read(self, tmp_path: Path) -> Any:
        return read_feature_table(write_layer(tmp_path / "f.gpkg"), TABLE)

    def test_columns_and_their_field_types_in_file_order(self, tmp_path: Path) -> None:
        assert self.read(tmp_path).columns == COLUMNS

    def test_fids_and_values_exactly(self, tmp_path: Path) -> None:
        got = [(fid, values) for fid, _, values in self.read(tmp_path).rows]
        assert got == [(fid, values) for fid, _, values in ROWS]

    def test_each_geometry_equals_the_written_one_exactly(self, tmp_path: Path) -> None:
        for (_, got, _), (fid, want, _) in zip(self.read(tmp_path).rows, ROWS, strict=True):
            assert got is not None, fid
            assert got.equals_exact(want, 0.0), fid

    def test_the_crs_text_and_the_layer(self, tmp_path: Path) -> None:
        table = self.read(tmp_path)
        assert (table.crs, table.layer) == (UTM33, TABLE)

    def test_the_only_layer_is_read_when_none_is_named(self, tmp_path: Path) -> None:
        table = read_feature_table(write_layer(tmp_path / "f.gpkg"), None)
        assert [fid for fid, _, _ in table.rows] == [7, 9, 12]

    def test_the_table_is_the_promised_dataclass(self, tmp_path: Path) -> None:
        assert isinstance(self.read(tmp_path), FeatureTable)


# ---------------------------------------------------------------- P2


class TestABoxKeepsTheRowsNearIt:
    """P2: a box meeting rows 7 and 9 only gives those two rows."""

    def test_two_rows(self, tmp_path: Path) -> None:
        asked: list[str] = []

        def box_for(crs: str) -> tuple[float, float, float, float]:
            asked.append(crs)
            return (X0 - 1.0, Y0 - 1.0, X0 + 111.0, Y0 + 11.0)

        table = read_feature_table(write_layer(tmp_path / "f.gpkg"), TABLE, box_for)
        assert [fid for fid, _, _ in table.rows] == [7, 9]
        assert [values for _, _, values in table.rows] == [ROWS[0][2], ROWS[1][2]]
        assert asked == [UTM33]


# ---------------------------------------------------------------- GeoPackage types


@pytest.mark.parametrize(
    ("declared", "field"),
    [
        ("INTEGER", "esriFieldTypeInteger"),
        ("INT", "esriFieldTypeInteger"),
        ("MEDIUMINT", "esriFieldTypeInteger"),
        ("SMALLINT", "esriFieldTypeInteger"),
        ("TINYINT", "esriFieldTypeInteger"),
        ("BOOLEAN", "esriFieldTypeInteger"),
        ("REAL", "esriFieldTypeDouble"),
        ("DOUBLE", "esriFieldTypeDouble"),
        ("FLOAT", "esriFieldTypeDouble"),
        ("TEXT", "esriFieldTypeString"),
        ("TEXT(12)", "esriFieldTypeString"),
        ("DATE", "esriFieldTypeString"),
        ("DATETIME", "esriFieldTypeString"),
    ],
)
def test_each_declared_column_type_has_its_field_type(
    tmp_path: Path, declared: str, field: str
) -> None:
    path = write_layer(
        tmp_path / "f.gpkg", extra_sql=[f'ALTER TABLE {TABLE} ADD COLUMN "v" {declared}']
    )
    table = read_feature_table(path, TABLE)
    assert table.columns == (*COLUMNS, ("v", field))


# ---------------------------------------------------------------- P4


def test_p4_a_blob_column_is_refused_by_name(tmp_path: Path) -> None:
    path = write_layer(
        tmp_path / "f.gpkg", extra_sql=[f'ALTER TABLE {TABLE} ADD COLUMN "raw" BLOB']
    )
    with pytest.raises(FeatureError, match=r"\braw\b"):
        read_feature_table(path, TABLE)


def test_an_undeclared_column_type_is_refused_by_name(tmp_path: Path) -> None:
    path = write_layer(
        tmp_path / "f.gpkg", extra_sql=[f'ALTER TABLE {TABLE} ADD COLUMN "odd" GEOMETRYX']
    )
    with pytest.raises(FeatureError, match=r"\bodd\b"):
        read_feature_table(path, TABLE)


# ---------------------------------------------------------------- P3


class TestAGeoJsonLayer:
    """P3: keys in order of first appearance, types from the values, fid the
    1-based position, a feature's own `id` not kept."""

    def test_columns_values_and_fids(self, tmp_path: Path) -> None:
        path = write_geojson(
            tmp_path / "f.geojson", [{"kind": "a"}, {"kind": "b", "n": 2}], ids=[50, 60]
        )
        table = read_feature_table(path, None)
        assert table.columns == (("kind", "esriFieldTypeString"), ("n", "esriFieldTypeInteger"))
        assert [(fid, values) for fid, _, values in table.rows] == [(1, ("a", None)), (2, ("b", 2))]
        assert table.layer is None

    def test_the_geometries_are_read(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.json", [{"kind": "a"}])
        ((_, geometry, _),) = read_feature_table(path, None).rows
        ring = [(10.0, 60.0), (10.1, 60.0), (10.1, 60.1), (10.0, 60.0)]
        assert geometry is not None and geometry.equals_exact(Polygon(ring), 0.0)

    def test_a_key_both_text_and_integer_is_refused_by_name(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.geojson", [{"code": "1"}, {"code": 1}])
        with pytest.raises(FeatureError, match=r"\bcode\b"):
            read_feature_table(path, None)

    def test_integers_and_floats_make_a_double(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.geojson", [{"x": 2}, {"x": 0.5}])
        table = read_feature_table(path, None)
        assert table.columns == (("x", "esriFieldTypeDouble"),)
        assert [values for _, _, values in table.rows] == [(2,), (0.5,)]

    def test_booleans_make_an_integer(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.geojson", [{"ok": True}, {"ok": False}])
        table = read_feature_table(path, None)
        assert table.columns == (("ok", "esriFieldTypeInteger"),)
        assert [values for _, _, values in table.rows] == [(1,), (0,)]

    def test_a_nested_value_is_refused_by_name(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.geojson", [{"deep": {"a": 1}}])
        with pytest.raises(FeatureError, match=r"\bdeep\b"):
            read_feature_table(path, None)

    def test_the_box_is_not_used(self, tmp_path: Path) -> None:
        def box_for(crs: str) -> tuple[float, float, float, float]:
            raise AssertionError(f"box_for called with {crs}")

        path = write_geojson(tmp_path / "f.geojson", [{"kind": "a"}, {"kind": "b"}])
        assert len(read_feature_table(path, None, box_for).rows) == 2


def test_an_unknown_suffix_is_refused_by_name(tmp_path: Path) -> None:
    path = tmp_path / "f.shp"
    path.write_bytes(b"")
    with pytest.raises(FeatureError, match=r"f\.shp"):
        read_feature_table(path, None)


# ---------------------------------------------------------------- refusals and no columns


def test_a_missing_geopackage_is_refused_by_name(tmp_path: Path) -> None:
    path = tmp_path / "gone.gpkg"
    with pytest.raises(FeatureError, match=r"^gone\.gpkg: "):
        read_feature_table(path, TABLE)
    assert not path.exists(), "the read-only open created the file"


def test_a_file_that_is_not_a_geopackage_is_refused_by_name(tmp_path: Path) -> None:
    path = tmp_path / "text.gpkg"
    path.write_text("not a database, only text\n" * 8, "utf-8")
    with pytest.raises(FeatureError, match=r"^text\.gpkg: not a readable GeoPackage"):
        read_feature_table(path, TABLE)


@pytest.mark.parametrize("name", ["f.geojson", "f.json"])
def test_a_layer_given_with_geojson_is_refused(tmp_path: Path, name: str) -> None:
    path = write_geojson(tmp_path / name, [{"kind": "a"}])
    with pytest.raises(FeatureError) as caught:
        read_feature_table(path, TABLE)
    assert str(caught.value) == f"{name}: a layer applies to a GeoPackage only"


def test_a_geopackage_layer_with_no_attribute_columns(tmp_path: Path) -> None:
    rows = tuple((fid, geometry, ()) for fid, geometry, _ in ROWS)
    table = read_feature_table(write_layer(tmp_path / "f.gpkg", columns=(), rows=rows), TABLE)
    assert table.columns == ()
    assert [(fid, values) for fid, _, values in table.rows] == [(7, ()), (9, ()), (12, ())]


def test_a_geojson_layer_with_no_properties(tmp_path: Path) -> None:
    table = read_feature_table(write_geojson(tmp_path / "f.geojson", [{}, {}]), None)
    assert table.columns == ()
    assert [(fid, values) for fid, _, values in table.rows] == [(1, ()), (2, ())]


# ---------------------------------------------------------------- P5


def local_entry(path: Path) -> RemoteSource:
    return RemoteSource.model_validate(
        {"id": "loc", "kind": "local", "path": str(path), "layer": TABLE,
         "credit": "Test layer", "licence_note": "Test licence"}
    )  # fmt: skip


def remote_entry() -> RemoteSource:
    return RemoteSource.model_validate(
        {"id": "far", "kind": "arcgis-features", "crs": UTM33, "fields": ["ID"],
         "url": "https://maps.example.org/arcgis/rest/services/T/MapServer/0",
         "credit": "Test layer, {url} (accessed on {date})", "licence_note": "Test licence"}
    )  # fmt: skip


class TestFeaturesFile:
    """P5: a features source's file, the way `mesh` finds it."""

    def test_a_local_entry_is_its_own_path(self, tmp_path: Path) -> None:
        path = write_layer(tmp_path / "f.gpkg")
        assert features_file(local_entry(path), None) == path

    def test_a_remote_entry_is_its_copy_in_the_cache(self, tmp_path: Path) -> None:
        copy = tmp_path / "cache" / "far" / "features.gpkg"
        copy.parent.mkdir(parents=True)
        write_layer(copy)
        assert features_file(remote_entry(), tmp_path / "cache") == copy

    def test_a_missing_copy_is_refused_by_name(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match=re.escape("far/features.gpkg is not in the cache")):
            features_file(remote_entry(), tmp_path / "cache")

    def test_a_missing_cache_is_refused(self) -> None:
        with pytest.raises(ValueError, match="no cache directory"):
            features_file(remote_entry(), None)

    def test_the_cache_root_moves_to_sources(self) -> None:
        assert cli.cache_root is sources.cache_root
