"""`tin_engine.feature_input`: sources, class maps, CRS, region and clip (16b-1).

`docs/increments/16b-terrain-polygons.md` R1, R4-R7, "Degeneracy policy" and
"Tests for @tester" (class maps, CRS, region and pre-clip, clip degeneracies),
with Ola's Q6 (b) (the legacy GML, read with the `clc18_kode` map).

Pinned by this suite (see "Pinned by the red suite (16b-1/2)"):

- `ClassMap(name, attribute, classes, otherwise, notice)`, frozen, extra
  fields forbidden. `CLASS_MAPS` maps `property`, `corine`, `corine-water`
  and `clc18_kode` to the built-in maps. `clc18_kode` is `corine` keyed on
  the Norwegian GML's attribute. Every CORINE map's `notice` names the
  Copernicus Land Monitoring Service and says the data were modified;
  `property`'s is empty.
- `FeatureSource(path, class_map, layer=None, crs=None)` and
  `FeatureRequest(sources, vocabulary=DEFAULT_VOCABULARY)`, frozen.
- `open_features(request, domain, dem_crs) -> FeatureSet`, `domain` already in
  the DEM's CRS. `FeatureSet.features` is a tuple of `TerrainFeature` (`fid`,
  `mask`, `lines`) in source order; `FeatureSet.outside` counts features read
  and then wholly clipped away; `FeatureSet.empty` counts empty geometries
  skipped.
- A GeoJSON feature's `fid` is its `id` member, else its position in the file
  (from 0). A GeoPackage row's is its primary key; a GML feature's its `fid`
  attribute.
- By suffix: `.geojson` and `.json`, `.gpkg`, `.gml`; anything else refused.
- `FeatureError`, a `ValueError`, for every refusal; it names the file, and
  for a feature its `fid`.
- The lines of a polygon feature: its exterior's pieces, then each hole's, part
  by part for a `MultiPolygon`.

Committed red at `e99c8ea` (amended at `3990449` and `972312c`):
`tin_engine.feature_input` did not exist yet, and because it is imported lazily
each test failed on its own. It landed in `5079da8` and the suite has been
green since.
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import closing
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely
from pyproj import CRS, Geod, Transformer
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
    mapping,
)

import feature_fixtures as ff
from crs_fixtures import axes_swapped, refuse_point_moves
from feature_fixtures import (
    CROSSING,
    UTM33,
    Feat,
    all_vertices,
    at,
    domain_of,
    moved,
    on_boundary,
    open_one,
    square,
    within,
    write_geojson,
)
from gpkg_fixtures import EXTRACT, LEGACY_GML, Layer, Row, needs_rtree, write_gpkg
from test_cli_mesh_domain import quarter_circle
from test_io_gml import document, member
from tin_engine.features import DEFAULT_VOCABULARY, EdgeProperty, EdgeVocabulary

V = DEFAULT_VOCABULARY
LAEA = "EPSG:3035"
WATER_CODES = ("511", "512", "521", "522", "523")
BOX = domain_of(square(0, 0, 300, 300))
INNER = square(100, 100, 200, 200)
#: `INNER` as a GeoJSON geometry, and a `crs` member naming the DEM's CRS.
GEOMETRY = json.loads(json.dumps(mapping(INNER)))
MEMBER = {"type": "name", "properties": {"name": UTM33}}


@pytest.fixture(scope="module")
def fi() -> ModuleType:
    return ff.feature_input()


def one(tmp_path: Path, geometry: Any, props: dict[str, Any] | None = None, **kw: Any) -> Any:
    path = write_geojson(
        tmp_path / "f.geojson", [Feat("a", geometry, props or {"property": "road"})]
    )
    return open_one(path, kw.pop("domain", BOX), kw.pop("map_name", "property"), **kw)


def lines_of(fs: Any) -> list[LineString]:
    return [line for f in fs.features for line in f.lines]


# --------------------------------------------------------------- class maps


class TestClassMaps:
    def test_the_built_in_maps(self, fi: ModuleType) -> None:
        assert set(fi.CLASS_MAPS) == {"property", "corine", "corine-water", "clc18_kode"}
        assert {name: m.attribute for name, m in fi.CLASS_MAPS.items()} == {
            "property": "property",
            "corine": "Code_18",
            "corine-water": "Code_18",
            "clc18_kode": "clc18_kode",
        }
        assert all(m.name == name for name, m in fi.CLASS_MAPS.items())

    def test_the_corine_maps_carry_the_attribution(self, fi: ModuleType) -> None:
        """R10 and "Test data": a mesh built from CORINE says so, and says the
        data were modified."""
        for name in ("corine", "corine-water", "clc18_kode"):
            notice = fi.CLASS_MAPS[name].notice
            assert "Copernicus Land Monitoring Service" in notice, name
            assert "modified" in notice.lower(), name
        assert fi.CLASS_MAPS["property"].notice == ""

    def test_a_class_map_is_frozen_and_strict(self, fi: ModuleType) -> None:
        corine = fi.CLASS_MAPS["corine"]
        with pytest.raises((TypeError, ValueError)):
            corine.attribute = "other"
        with pytest.raises(ValueError):
            fi.ClassMap(name="x", attribute="a", classes={}, colour="red")

    def test_property_one_name(self, tmp_path: Path) -> None:
        fs = one(tmp_path, INNER, {"property": "road"})
        assert [f.mask for f in fs.features] == [V.mask("road")]

    def test_property_a_list_of_names_is_their_union(self, tmp_path: Path) -> None:
        fs = one(tmp_path, INNER, {"property": ["road", "water"]})
        assert [f.mask for f in fs.features] == [V.mask("road", "water")]

    def test_property_every_vocabulary_name_including_the_new_bits(self, tmp_path: Path) -> None:
        for prop in V.properties:
            fs = one(tmp_path, INNER, {"property": prop.name})
            assert [f.mask for f in fs.features] == [1 << prop.bit], prop.name

    def test_property_an_unknown_name_is_refused_naming_feature_and_value(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = write_geojson(tmp_path / "f.geojson", [Feat("f-17", INNER, {"property": "glacier"})])
        with pytest.raises(fi.FeatureError) as caught:
            open_one(path, BOX)
        assert "f-17" in str(caught.value) and "glacier" in str(caught.value)

    def test_property_missing_is_refused_naming_the_feature(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = write_geojson(tmp_path / "f.geojson", [Feat("f-18", INNER, {"name": "x"})])
        with pytest.raises(fi.FeatureError, match="f-18"):
            open_one(path, BOX)

    @pytest.mark.parametrize("code", WATER_CODES)
    def test_corine_water_codes_are_land_cover_and_water(self, tmp_path: Path, code: str) -> None:
        """MUTANT 3 (in `test_feature_chains_bits`), unit by unit."""
        fs = one(tmp_path, INNER, {"Code_18": code}, map_name="corine")
        assert [f.mask for f in fs.features] == [V.mask("land_cover", "water")]

    @pytest.mark.parametrize("code", ["111", "311", "412", "423", "999", "5", "51", "513", "530"])
    def test_corine_every_other_code_is_land_cover(self, tmp_path: Path, code: str) -> None:
        """`otherwise = ("land_cover",)`: an unlisted code, a malformed one
        and the 5xx codes the nomenclature does not have included."""
        fs = one(tmp_path, INNER, {"Code_18": code}, map_name="corine")
        assert [f.mask for f in fs.features] == [V.mask("land_cover")]

    @pytest.mark.parametrize("code", WATER_CODES)
    def test_corine_water_keeps_water_only(self, tmp_path: Path, code: str) -> None:
        fs = one(tmp_path, INNER, {"Code_18": code}, map_name="corine-water")
        assert [f.mask for f in fs.features] == [V.mask("water")]

    def test_corine_water_drops_land(self, tmp_path: Path) -> None:
        features = [
            Feat(1, INNER, {"Code_18": "311"}),
            Feat(2, square(10, 10, 50, 50), {"Code_18": "512"}),
        ]
        path = write_geojson(tmp_path / "f.geojson", features)
        fs = open_one(path, BOX, "corine-water")
        assert [f.fid for f in fs.features] == [2]

    @pytest.mark.parametrize(
        ("code", "names"), [("512", ("land_cover", "water")), ("322", ("land_cover",))]
    )
    def test_clc18_kode(self, tmp_path: Path, code: str, names: tuple[str, ...]) -> None:
        fs = one(tmp_path, INNER, {"clc18_kode": code}, map_name="clc18_kode")
        assert [f.mask for f in fs.features] == [V.mask(*names)]

    def test_a_map_naming_a_name_the_vocabulary_lacks_is_refused_at_the_start(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        """R4: checked when the run starts, even with no feature to map."""
        bad = fi.ClassMap(
            name="bad", attribute="kind", classes={"x": ("glacier",)}, otherwise="drop"
        )
        path = write_geojson(tmp_path / "f.geojson", [])
        request = fi.FeatureRequest(sources=(fi.FeatureSource(path=path, class_map=bad),))
        with pytest.raises(fi.FeatureError, match="glacier"):
            fi.open_features(request, BOX, UTM33)

    def test_an_otherwise_name_the_vocabulary_lacks_is_refused_too(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        bad = fi.ClassMap(name="bad", attribute="kind", classes={}, otherwise=("glacier",))
        path = write_geojson(tmp_path / "f.geojson", [])
        request = fi.FeatureRequest(sources=(fi.FeatureSource(path=path, class_map=bad),))
        with pytest.raises(fi.FeatureError, match="glacier"):
            fi.open_features(request, BOX, UTM33)

    def test_the_request_vocabulary_is_the_one_masks_come_from(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        mine = EdgeVocabulary(properties=(EdgeProperty(name="road", bit=9),))
        path = write_geojson(tmp_path / "f.geojson", [Feat("a", INNER, {"property": "road"})])
        own_map = fi.ClassMap(name="mine", attribute="property", classes={"road": ("road",)})
        request = fi.FeatureRequest(
            sources=(fi.FeatureSource(path=path, class_map=own_map),), vocabulary=mine
        )
        fs = fi.open_features(request, BOX, UTM33)
        assert [f.mask for f in fs.features] == [1 << 9]


# ---------------------------------------------------------------- sources


class TestSources:
    def test_geojson_order_and_ids(self, tmp_path: Path) -> None:
        features = [
            Feat("z", square(10, 10, 20, 20), {"property": "road"}),
            Feat(None, square(30, 10, 40, 20), {"property": "road"}),
            Feat(7, square(50, 10, 60, 20), {"property": "road"}),
        ]
        fs = open_one(write_geojson(tmp_path / "f.geojson", features), BOX)
        assert [f.fid for f in fs.features] == ["z", 1, 7]

    def test_json_suffix_is_geojson(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.json", [Feat("a", INNER, {"property": "road"})])
        assert len(open_one(path, BOX).features) == 1

    def test_an_unknown_suffix_is_refused_naming_it(self, tmp_path: Path, fi: ModuleType) -> None:
        path = tmp_path / "f.shp"
        path.write_bytes(b"")
        with pytest.raises(fi.FeatureError, match=r"\.shp"):
            open_one(path, BOX)

    def test_a_layer_on_geojson_is_refused(self, tmp_path: Path, fi: ModuleType) -> None:
        path = write_geojson(tmp_path / "f.geojson", [])
        with pytest.raises(fi.FeatureError):
            open_one(path, BOX, layer="x")

    def test_unparseable_geojson_is_refused_naming_the_file(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = tmp_path / "broken.geojson"
        path.write_text('{"type": "FeatureCollection", "features": [')
        with pytest.raises(fi.FeatureError, match=r"broken\.geojson"):
            open_one(path, BOX)

    def test_a_missing_file_is_refused_naming_it(self, tmp_path: Path, fi: ModuleType) -> None:
        for name in ("absent.geojson", "absent.gpkg", "absent.gml"):
            with pytest.raises(fi.FeatureError, match=name.replace(".", r"\.")):
                open_one(tmp_path / name, BOX)

    def test_the_error_is_a_value_error(self, fi: ModuleType) -> None:
        assert issubclass(fi.FeatureError, ValueError)

    @needs_rtree
    def test_geopackage_rows_in_key_order_with_their_values(self, tmp_path: Path) -> None:
        polys = {9: square(10, 10, 20, 20), 2: square(30, 10, 40, 20), 5: square(50, 10, 60, 20)}
        rows = [Row(pk, g, {"Code_18": "512" if pk == 5 else "311"}) for pk, g in polys.items()]
        path = write_gpkg(tmp_path / "f.gpkg", [Layer("clc", 25833, rows)])
        fs = open_one(path, BOX, "corine")
        assert [f.fid for f in fs.features] == [2, 5, 9]
        assert [f.mask for f in fs.features] == [
            V.mask("land_cover"),
            V.mask("land_cover", "water"),
            V.mask("land_cover"),
        ]

    def test_geopackage_several_layers_need_one_named(self, tmp_path: Path, fi: ModuleType) -> None:
        rows = [Row(1, INNER, {"Code_18": "311"})]
        path = write_gpkg(
            tmp_path / "f.gpkg",
            [Layer("a", 25833, rows, rtree=False), Layer("b", 25833, rows, rtree=False)],
        )
        with pytest.raises(fi.FeatureError, match=r"f\.gpkg"):
            open_one(path, BOX, "corine")
        assert len(open_one(path, BOX, "corine", layer="b").features) == 1

    def test_geopackage_a_missing_attribute_column_is_refused(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = write_gpkg(
            tmp_path / "f.gpkg", [Layer("a", 25833, [Row(1, INNER)], columns=("kind",))]
        )
        with pytest.raises(fi.FeatureError, match="Code_18"):
            open_one(path, BOX, "corine")

    def test_geopackage_refusals_become_feature_errors_naming_the_file(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = write_gpkg(
            tmp_path / "f.gpkg", [Layer("a", 25833, [Row(1, INNER)])], application_id=1
        )
        with pytest.raises(fi.FeatureError, match=r"f\.gpkg"):
            open_one(path, BOX, "corine")

    def test_a_file_that_is_not_sqlite_is_refused(self, tmp_path: Path, fi: ModuleType) -> None:
        path = tmp_path / "f.gpkg"
        path.write_bytes(b"this is not a database" * 100)
        with pytest.raises(fi.FeatureError, match=r"f\.gpkg"):
            open_one(path, BOX, "corine")

    def test_the_geopackage_is_left_unchanged(self, tmp_path: Path) -> None:
        path = write_gpkg(
            tmp_path / "f.gpkg", [Layer("a", 25833, [Row(1, INNER, {"Code_18": "311"})])]
        )
        before = path.read_bytes()
        open_one(path, BOX, "corine")
        assert path.read_bytes() == before
        with closing(sqlite3.connect(path)) as con:
            assert con.execute("PRAGMA integrity_check").fetchone()[0] == "ok"

    def test_gml_with_the_clc18_kode_map(self, tmp_path: Path) -> None:
        lonlat = moved(INNER, UTM33, "EPSG:4326")
        path = tmp_path / "f.gml"
        path.write_bytes(document(member("sql_statement.1", lonlat, "512")).getvalue())
        fs = open_one(path, BOX, "clc18_kode")
        assert [f.fid for f in fs.features] == ["sql_statement.1"]
        assert [f.mask for f in fs.features] == [V.mask("land_cover", "water")]

    def test_gml_without_a_crs_needs_one_given(self, tmp_path: Path, fi: ModuleType) -> None:
        path = tmp_path / "f.gml"
        path.write_bytes(document(member("s.1", INNER, "311", srs=None)).getvalue())
        with pytest.raises(fi.FeatureError, match=r"f\.gml"):
            open_one(path, BOX, "clc18_kode")
        assert len(open_one(path, BOX, "clc18_kode", crs=UTM33).features) == 1

    def test_gml_refusals_become_feature_errors_naming_the_file(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = tmp_path / "f.gml"
        path.write_bytes(b"<not gml")
        with pytest.raises(fi.FeatureError, match=r"f\.gml"):
            open_one(path, BOX, "clc18_kode")

    def test_the_legacy_gml_over_the_benchmark_tile(self) -> None:
        """Q6 (b): the committed legacy file, through the `clc18_kode` map,
        clipped to the quarter circle. Every feature's bits are `corine`'s for
        its code, and the sea polygon (`sql_statement.11029`, 523) is among
        them."""
        domain = domain_of(Polygon(quarter_circle()))
        fs = open_one(LEGACY_GML, domain, "clc18_kode")
        fids = {f.fid for f in fs.features}
        assert "sql_statement.11029" in fids
        assert len(fs.features) >= 10
        assert {f.mask for f in fs.features} == {
            V.mask("land_cover"),
            V.mask("land_cover", "water"),
        }
        assert within(domain.polygon, all_vertices(lines_of(fs))).max() <= CROSSING


# ---------------------------------------------------------------- geometry


class TestGeometry:
    @pytest.mark.parametrize(
        "geometry",
        [
            pytest.param(Point(*at(150, 150)), id="point"),
            pytest.param(MultiPoint([at(150, 150), at(160, 160)]), id="multipoint"),
            pytest.param(GeometryCollection([INNER]), id="geometrycollection"),
        ],
    )
    def test_loose_points_and_collections_are_refused_naming_the_feature(
        self, tmp_path: Path, fi: ModuleType, geometry: Any
    ) -> None:
        path = write_geojson(tmp_path / "f.geojson", [Feat("pt-3", geometry, {"property": "road"})])
        with pytest.raises(fi.FeatureError, match="pt-3"):
            open_one(path, BOX)

    def test_a_point_in_gml_is_refused_too(self, tmp_path: Path, fi: ModuleType) -> None:
        path = tmp_path / "f.gml"
        path.write_bytes(document(member("s.9", Point(15.0, 59.5), "311")).getvalue())
        with pytest.raises(fi.FeatureError, match=r"s\.9"):
            open_one(path, BOX, "clc18_kode")

    def test_z_is_dropped(self, tmp_path: Path) -> None:
        ring = [
            (*at(100, 100), 5.0),
            (*at(200, 100), 6.0),
            (*at(200, 200), 7.0),
            (*at(100, 100), 5.0),
        ]
        fs = one(tmp_path, Polygon(ring))
        (line,) = lines_of(fs)
        assert not line.has_z
        assert set(line.coords) == {p[:2] for p in ring}

    @pytest.mark.parametrize(
        "geometry",
        [
            pytest.param({"type": "Polygon", "coordinates": []}, id="empty-polygon"),
            pytest.param({"type": "LineString", "coordinates": []}, id="empty-line"),
            pytest.param(None, id="null"),
        ],
    )
    def test_an_empty_geometry_is_skipped_and_counted(self, tmp_path: Path, geometry: Any) -> None:
        features = [
            Feat("e", geometry, {"property": "road"}),
            Feat("f", INNER, {"property": "road"}),
        ]
        fs = open_one(write_geojson(tmp_path / "f.geojson", features), BOX)
        assert [f.fid for f in fs.features] == ["f"]
        assert fs.empty == 1

    def test_an_empty_flag_row_is_skipped_and_counted(self, tmp_path: Path) -> None:
        from gpkg_fixtures import blob

        rows = [
            Row(1, INNER, {"Code_18": "311"}),
            Row(2, None, {"Code_18": "311"}, raw=blob(None, 25833, empty=True)),
        ]
        fs = open_one(
            write_gpkg(tmp_path / "f.gpkg", [Layer("a", 25833, rows, rtree=False)]), BOX, "corine"
        )
        assert [f.fid for f in fs.features] == [1]
        assert fs.empty == 1

    def test_a_self_intersecting_ring_is_linework_not_refused(self, tmp_path: Path) -> None:
        """R1: feature polygons are not validated. Its source vertices are all
        kept and its length is unchanged; a vertex GEOS may add at the
        self-crossing (150, 150) is not pinned either way."""
        bowtie = Polygon([at(100, 100), at(200, 200), at(200, 100), at(100, 200)])
        assert not bowtie.is_valid
        fs = one(tmp_path, bowtie)
        assert len(fs.features) == 1
        got = {tuple(p) for p in all_vertices(lines_of(fs))}
        assert got - set(bowtie.exterior.coords) <= {at(150, 150)}
        assert got >= set(bowtie.exterior.coords)
        assert sum(line.length for line in lines_of(fs)) == pytest.approx(bowtie.length, abs=1e-6)


# ------------------------------------------------------------------- CRS


class TestCrs:
    def test_geojson_without_a_crs_member_is_wgs84(self, tmp_path: Path) -> None:
        lonlat = moved(INNER, UTM33, "EPSG:4326")
        path = write_geojson(
            tmp_path / "f.geojson", [Feat("a", lonlat, {"property": "road"})], crs=None
        )
        (line,) = lines_of(open_one(path, BOX))
        assert set(line.coords) == set(moved(lonlat, "EPSG:4326", UTM33).exterior.coords)

    def test_geojson_in_the_dems_crs_is_not_moved(self, tmp_path: Path) -> None:
        (line,) = lines_of(one(tmp_path, INNER))
        assert set(line.coords) == set(INNER.exterior.coords)

    def test_the_dems_epsg_code_by_definition_is_not_moved(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Audit PR B, red test 9: a source spelt as the WKT of the DEM's
        EPSG:25833, without its ID and axes swapped, is in the DEM's CRS (`crs.same_crs`), so its
        geometries are those of the source spelt `EPSG:25833`, and no point
        moves, neither for the source region nor for the features."""
        leaving = LineString([at(150, 150), at(400, 170)])  # clipped at the domain
        features = [Feat(fid, g, {"property": "road"}) for fid, g in (("a", INNER), ("b", leaving))]
        by_code = open_one(write_geojson(tmp_path / "code.geojson", features), BOX)
        spelt = write_geojson(tmp_path / "proj.geojson", features, crs=axes_swapped(25833))
        refuse_point_moves(monkeypatch)
        by_proj = open_one(spelt, BOX)
        assert [f.fid for f in by_proj.features] == [f.fid for f in by_code.features]
        assert [line.coords[:] for line in lines_of(by_proj)] == [
            line.coords[:] for line in lines_of(by_code)
        ]

    def test_a_given_crs_agreeing_with_the_files_is_accepted(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.geojson", [Feat("a", INNER, {"property": "road"})])
        assert len(open_one(path, BOX, crs="EPSG:25833").features) == 1
        assert len(open_one(path, BOX, crs=CRS.from_epsg(25833).to_wkt()).features) == 1

    def test_a_given_crs_disagreeing_with_the_files_is_refused(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = write_geojson(tmp_path / "f.geojson", [Feat("a", INNER, {"property": "road"})])
        with pytest.raises(fi.FeatureError, match=r"f\.geojson"):
            open_one(path, BOX, crs="EPSG:3035")

    def test_a_given_crs_disagreeing_with_rfc_7946_is_refused(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        """domain.py's rule exactly: no `crs` member means EPSG:4326."""
        path = write_geojson(
            tmp_path / "f.geojson", [Feat("a", INNER, {"property": "road"})], crs=None
        )
        with pytest.raises(fi.FeatureError):
            open_one(path, BOX, crs="EPSG:25833")

    def test_an_unreadable_crs_is_refused(self, tmp_path: Path, fi: ModuleType) -> None:
        path = write_geojson(
            tmp_path / "f.geojson", [Feat("a", INNER, {"property": "road"})], crs="EPSG:0"
        )
        with pytest.raises(fi.FeatureError):
            open_one(path, BOX)

    def test_a_geopackage_in_3035(self, tmp_path: Path) -> None:
        laea = moved(INNER, UTM33, "EPSG:3035")
        path = write_gpkg(
            tmp_path / "f.gpkg", [Layer("a", 3035, [Row(1, laea, {"Code_18": "311"})])]
        )
        (line,) = lines_of(open_one(path, BOX, "corine"))
        assert set(line.coords) == set(moved(laea, "EPSG:3035", UTM33).exterior.coords)

    def test_a_geopackage_given_another_crs_is_refused(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = write_gpkg(
            tmp_path / "f.gpkg", [Layer("a", 3035, [Row(1, INNER, {"Code_18": "311"})])]
        )
        with pytest.raises(fi.FeatureError):
            open_one(path, BOX, "corine", crs="EPSG:25833")

    def test_a_vertex_with_no_image_is_refused_naming_the_feature(
        self, tmp_path: Path, fi: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """R5. Forced: the transform into the DEM's CRS gives `inf` for one
        vertex of feature `lost`, as pyproj does for a point it cannot map."""
        import tin_engine.crs as crs_module

        real = crs_module.reprojector
        laea = moved(Polygon([at(100, 100), at(200, 100), at(150, 150), at(100, 200)]), UTM33, LAEA)
        marker = laea.exterior.coords[2]  # in the source CRS

        def fake(src: Any, dst: Any) -> Any:
            move = real(src, dst)
            if CRS.from_user_input(dst) != CRS.from_user_input(UTM33):
                return move

            def apply(xy: Any) -> np.ndarray:
                out = move(xy)
                points = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
                out[np.all(np.isclose(points, marker, rtol=0, atol=1e-6), axis=1)] = np.inf
                return out

            return apply

        monkeypatch.setattr(crs_module, "reprojector", fake)
        monkeypatch.setattr(fi, "reprojector", fake, raising=False)
        path = write_geojson(
            tmp_path / "f.geojson", [Feat("lost", laea, {"property": "road"})], crs=LAEA
        )
        with pytest.raises(fi.FeatureError, match="lost"):
            open_one(path, BOX)


# --------------------------------------------------------- region and clip


class TestOneGeojsonRule:
    """Audit PR C (`docs/increments/python-audit.md`, section 11, red test 2):
    `read_source`'s GeoJSON branch reads through `io.geojson.read_collection`
    with RFC 7946's default, and a refusal reaches the existing handler as
    `<file name>: <words>`. Each row here differs from the code before it."""

    @staticmethod
    def rows(fi: ModuleType, path: Path) -> Any:
        return fi.read_source(path, None, "property", lambda _crs: (0.0, 0.0, 0.0, 0.0))

    @staticmethod
    def written(tmp_path: Path, doc: object) -> Path:
        path = tmp_path / "f.geojson"
        path.write_text(json.dumps(doc), encoding="utf-8")
        return path

    @staticmethod
    def collection(**members: Any) -> dict[str, Any]:
        feature = {"type": "Feature", "properties": {"property": "road"}, "geometry": GEOMETRY}
        return {"type": "FeatureCollection", "features": [feature]} | members

    @pytest.mark.parametrize(
        ("member", "says"),
        [
            (None, "the crs member is null; the file must name its CRS"),
            ({}, "the crs member has no name; it must name the CRS"),
        ],
        ids=["null", "empty_object"],
    )
    def test_a_member_naming_nothing_is_refused_naming_the_file(
        self, tmp_path: Path, fi: ModuleType, member: Any, says: str
    ) -> None:
        """Before: both read as EPSG:4326, as if the member were absent."""
        path = self.written(tmp_path, self.collection(crs=member))
        with pytest.raises(fi.FeatureError) as info:
            self.rows(fi, path)
        assert str(info.value) == f"f.geojson: {says}"

    def test_a_feature_file_reads(self, tmp_path: Path, fi: ModuleType) -> None:
        """Before: refused as not a FeatureCollection. Its fid is its position."""
        doc = {"type": "Feature", "properties": {"property": "road"}, "geometry": GEOMETRY}
        read = self.rows(fi, self.written(tmp_path, doc | {"crs": MEMBER}))
        ((fid, geometry, value),) = read.rows
        assert (fid, value, read.crs) == (0, "road", UTM33)
        assert geometry.equals(INNER)

    def test_a_bare_geometry_file_reads(self, tmp_path: Path, fi: ModuleType) -> None:
        """Before: refused. One feature with no properties, so no class value."""
        read = self.rows(fi, self.written(tmp_path, GEOMETRY | {"crs": MEMBER}))
        ((fid, geometry, value),) = read.rows
        assert (fid, value, read.crs) == (0, None, UTM33)
        assert geometry.equals(INNER)

    def test_features_that_are_not_a_list_are_refused(self, tmp_path: Path, fi: ModuleType) -> None:
        """Before: `"features": {}` read as no features."""
        path = self.written(tmp_path, self.collection(crs=MEMBER, features={}))
        with pytest.raises(fi.FeatureError) as info:
            self.rows(fi, path)
        assert str(info.value) == (
            "f.geojson: no features list; the file is not a FeatureCollection"
        )

    def test_a_numeric_name_is_read_as_its_digits(self, tmp_path: Path, fi: ModuleType) -> None:
        """Before: the int 4326 reached `FeatureSet(crs=...)`, whose field is
        `tuple[str, ...]`, and pydantic refused it."""
        member = {"type": "name", "properties": {"name": 4326}}
        read = self.rows(fi, self.written(tmp_path, self.collection(crs=member)))
        assert read.crs == "4326"
        assert isinstance(read.crs, str)

    def test_a_numeric_name_opens_as_a_feature_set(self, tmp_path: Path) -> None:
        """The same file through `open_features`: its `crs` is the text."""
        lonlat = moved(INNER, UTM33, "EPSG:4326")
        feature = {"type": "Feature", "properties": {"property": "road"}}
        doc = {
            "type": "FeatureCollection",
            "crs": {"type": "name", "properties": {"name": 4326}},
            "features": [feature | {"geometry": mapping(lonlat)}],
        }
        fs = open_one(self.written(tmp_path, doc), BOX)
        assert fs.crs == ("4326",)
        assert len(fs.features) == 1


#: Every empty JSON value but `null` (`docs/increments/python-audit.md`,
#: section 11, the ruling after code review round 2).
EMPTY_NOT_NULL = [
    pytest.param("", id="empty_string"),
    pytest.param(0, id="zero"),
    pytest.param(False, id="false"),
    pytest.param([], id="empty_list"),
    pytest.param({}, id="empty_object"),
]


class TestAnEmptyGeometryIsNull:
    """Audit PR C, section 11's ruling after code review round 2: `read_source`
    treats every empty `geometry` as `null`. Before, `""`, `0`, `false`, `[]`
    and `{}` crashed `--features` and `catchment --lakes` with an
    `AttributeError` (`'str' object has no attribute 'is_empty'`). Each test
    compares with the same file whose `geometry` is `null`, so it pins no
    outcome of its own: only that the two read alike."""

    @staticmethod
    def written(tmp_path: Path, geometry: Any, shape: str) -> Path:
        """A top-level `Feature`, or a collection whose first feature has
        `geometry` and whose second is `INNER`; both with a `crs` member."""
        feature = {"type": "Feature", "properties": {"property": "water"}, "geometry": geometry}
        if shape == "feature":
            doc = feature | {"crs": MEMBER}
        else:
            kept = feature | {"geometry": GEOMETRY}
            doc = {"type": "FeatureCollection", "crs": MEMBER, "features": [feature, kept]}
        path = tmp_path / f"{shape}-{json.dumps(geometry)}.geojson"
        path.write_text(json.dumps(doc), encoding="utf-8")
        return path

    @staticmethod
    def rows(fi: ModuleType, path: Path) -> Any:
        read = fi.read_source(path, None, "property", lambda _crs: (0.0, 0.0, 0.0, 0.0))
        return read.crs, read.rows

    @staticmethod
    def opened(path: Path) -> Any:
        fs = open_one(path, BOX)
        return ([f.fid for f in fs.features], fs.outside, fs.clipped, fs.empty, fs.crs, fs.counts)

    @staticmethod
    def lakes(fi: ModuleType, path: Path) -> Any:
        lakes, crs = fi.read_lake_polygons(path, None, at(150, 150), UTM33)
        return lakes, crs

    @pytest.mark.parametrize("shape", ["feature", "collection"])
    @pytest.mark.parametrize("empty", EMPTY_NOT_NULL)
    def test_read_source_reads_it_as_null(
        self, tmp_path: Path, fi: ModuleType, empty: Any, shape: str
    ) -> None:
        null = self.rows(fi, self.written(tmp_path, None, shape))
        assert self.rows(fi, self.written(tmp_path, empty, shape)) == null

    @pytest.mark.parametrize("shape", ["feature", "collection"])
    @pytest.mark.parametrize("empty", EMPTY_NOT_NULL)
    def test_features_reads_it_as_null(self, tmp_path: Path, empty: Any, shape: str) -> None:
        """`--features`, through `open_features`."""
        null = self.opened(self.written(tmp_path, None, shape))
        assert self.opened(self.written(tmp_path, empty, shape)) == null

    @pytest.mark.parametrize("shape", ["feature", "collection"])
    @pytest.mark.parametrize("empty", EMPTY_NOT_NULL)
    def test_catchment_lakes_reads_it_as_null(
        self, tmp_path: Path, fi: ModuleType, empty: Any, shape: str
    ) -> None:
        """`catchment --lakes`, through `read_lake_polygons`."""
        null = self.lakes(fi, self.written(tmp_path, None, shape))
        assert self.lakes(fi, self.written(tmp_path, empty, shape)) == null


class TestRenamedAndMoved:
    """Audit PR C, red tests 5 and 7: the any-source lake reader is named by
    what it reads, and `TerrainFeature` lives in `features` (layer 0), which
    `feature_input` imports it back from."""

    def test_read_lake_polygons_replaces_read_lakes(self, fi: ModuleType) -> None:
        assert callable(getattr(fi, "read_lake_polygons", None))
        assert not hasattr(fi, "read_lakes")

    def test_read_lake_polygons_skips_what_is_not_a_polygon(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        """Increment 22's behaviour, unchanged by the rename."""
        features = [
            Feat("lake", INNER, {"property": "water"}),
            Feat("shore", LineString([at(0, 0), at(300, 0)]), {"property": "road"}),
        ]
        path = write_geojson(tmp_path / "l.geojson", features)
        lakes, crs = fi.read_lake_polygons(path, None, at(150, 150), UTM33)
        assert crs == UTM33
        (lake,) = lakes
        assert lake.equals(INNER)

    def test_terrain_feature_is_one_class(self, fi: ModuleType) -> None:
        import tin_engine.features as features

        assert features.TerrainFeature is fi.TerrainFeature
        assert fi.TerrainFeature.__module__ == "tin_engine.features"


class TestClip:
    def test_a_ring_wholly_inside_is_one_closed_line(self, tmp_path: Path) -> None:
        (line,) = lines_of(one(tmp_path, INNER))
        assert line.coords[0] == line.coords[-1]
        assert set(line.coords) == set(INNER.exterior.coords)

    def test_a_polygon_gives_its_exterior_then_its_holes(self, tmp_path: Path) -> None:
        hole = square(130, 130, 170, 170)
        poly = Polygon(INNER.exterior.coords, [hole.exterior.coords])
        lines = lines_of(one(tmp_path, poly))
        assert [set(line.coords) for line in lines] == [
            set(INNER.exterior.coords),
            set(hole.exterior.coords),
        ]

    def test_a_multipolygon_gives_its_parts_in_order(self, tmp_path: Path) -> None:
        a, b = square(200, 200, 250, 250), square(20, 20, 60, 60)
        lines = lines_of(one(tmp_path, MultiPolygon([a, b])))
        assert [set(line.coords) for line in lines] == [
            set(a.exterior.coords),
            set(b.exterior.coords),
        ]

    def test_a_multilinestring_gives_its_parts(self, tmp_path: Path) -> None:
        a = LineString([at(10, 10), at(90, 90)])
        b = LineString([at(10, 290), at(90, 210)])
        lines = lines_of(one(tmp_path, MultiLineString([a, b])))
        assert [list(line.coords) for line in lines] == [list(a.coords), list(b.coords)]

    def test_a_ring_crossing_the_boundary_is_open_pieces(self, tmp_path: Path) -> None:
        """R6: the part outside is dropped, and no edge is added along the
        boundary. The pieces' length is the ring's length inside the domain,
        and each piece ends on the boundary or at a source vertex (GEOS may
        split a piece at the ring's own start, which the noder rejoins)."""
        crossing = square(250, 100, 350, 200)
        lines = lines_of(one(tmp_path, crossing))
        assert lines and all(line.coords[0] != line.coords[-1] for line in lines)
        ends = [p for line in lines for p in (line.coords[0], line.coords[-1])]
        loose = [p for p in ends if p not in set(crossing.exterior.coords)]
        assert loose and on_boundary(BOX.polygon, np.array(loose)).max() <= CROSSING
        inside = shapely.intersection(LineString(crossing.exterior.coords), BOX.polygon).length
        assert sum(line.length for line in lines) == pytest.approx(inside, abs=1e-6)
        assert sum(line.length for line in lines) == pytest.approx(200.0, abs=1e-6)

    def test_a_line_leaving_the_domain_loses_its_exterior_part(self, tmp_path: Path) -> None:
        """Increment 8's `wall-leaves-domain`, before the engine."""
        wall = LineString([at(150, 150), at(700, 150)])
        (line,) = lines_of(one(tmp_path, wall, {"property": "wall"}))
        assert line.coords[0] == at(150, 150)
        assert line.coords[-1] == pytest.approx(at(300, 150), abs=CROSSING)
        assert within(BOX.polygon, np.asarray(line.coords)).max() <= CROSSING

    def test_a_ring_touching_the_boundary_at_one_point_keeps_no_point_piece(
        self, tmp_path: Path
    ) -> None:
        diamond = Polygon([at(300, 150), at(250, 200), at(200, 150), at(250, 100)])
        lines = lines_of(one(tmp_path, diamond))
        assert all(isinstance(line, LineString) and len(line.coords) >= 2 for line in lines)
        assert sum(line.length for line in lines) == pytest.approx(diamond.length, abs=1e-6)

    def test_a_ring_outside_touching_at_one_point_is_dropped_and_counted(
        self, tmp_path: Path
    ) -> None:
        """The degenerate point piece is all that is left: no line survives."""
        outside = Polygon([at(300, 150), at(350, 200), at(400, 150), at(350, 100)])
        fs = one(tmp_path, outside)
        assert fs.features == () and fs.outside == 1

    def test_a_feature_edge_along_the_boundary_is_kept(self, tmp_path: Path) -> None:
        along = square(200, 0, 300, 100)  # its south and east edges are the domain's
        lines = lines_of(one(tmp_path, along))
        total = sum(line.length for line in lines)
        assert total == pytest.approx(along.length, abs=1e-6)

    def test_a_feature_in_a_domain_hole_is_dropped_and_counted(self, tmp_path: Path) -> None:
        holed = domain_of(Polygon(square(0, 0, 300, 300).exterior.coords, [INNER.exterior.coords]))
        fs = one(tmp_path, square(120, 120, 180, 180), domain=holed)
        assert fs.features == () and fs.outside == 1

    def test_a_feature_just_outside_is_dropped_and_counted(self, tmp_path: Path) -> None:
        fs = one(tmp_path, square(320, 100, 360, 200))
        assert fs.features == () and fs.outside == 1

    def test_a_feature_far_away_is_not_a_feature(self, tmp_path: Path) -> None:
        fs = one(tmp_path, square(50_000, 50_000, 50_100, 50_100))
        assert fs.features == ()

    def test_no_feature_left_is_legitimate(self, tmp_path: Path) -> None:
        fs = open_one(write_geojson(tmp_path / "f.geojson", []), BOX)
        assert fs.features == () and fs.outside == 0 and fs.empty == 0

    def test_a_domain_hole_cuts_rings_as_lines(self, tmp_path: Path) -> None:
        holed = domain_of(Polygon(square(0, 0, 300, 300).exterior.coords, [INNER.exterior.coords]))
        crossing = square(150, 150, 250, 250)
        lines = lines_of(one(tmp_path, crossing, domain=holed))
        assert lines and all(line.coords[0] != line.coords[-1] for line in lines)
        assert sum(line.length for line in lines) == pytest.approx(300.0, abs=1e-6)
        assert within(holed.polygon, all_vertices(lines)).max() <= CROSSING


# -------------------------------------------------------- the pre-clip (I4)


class TestPreClip:
    def test_i4_the_extract_gives_the_same_noded_graph_as_no_pre_clip(self) -> None:
        """I4 on the extract, relationally: `open_features` + `start_chains`
        through the engine give the same noded constraint edges, with the same
        masks, as this test's own pipeline with no region and no pre-clip
        (every row read whole, moved by pyproj, clipped as lines, the `corine`
        map written from the design), handed to `start_chains` as
        `TerrainFeature(fid=, mask=, lines=)`. Compared after the noder, where input
        order and piece splitting no longer matter."""
        domain = domain_of(Polygon(quarter_circle()))
        fs = open_one(EXTRACT, domain, "corine")
        produced = ff.run_engine(ff.start(domain, fs.features))

        rows = ff.extract_features(EXTRACT, "U2018_CLC2018_V2020_20u1")
        fi = ff.feature_input()
        oracle_features = []
        for pk, geometry, code in rows:
            lines = tuple(
                part
                for line in ff.boundary_lines(moved(geometry, LAEA))
                for part in shapely.get_parts(shapely.intersection(line, domain.polygon))
                if isinstance(part, LineString) and part.length > 0
            )
            if lines:
                oracle_features.append(
                    fi.TerrainFeature(fid=pk, mask=ff.corine_mask(code), lines=lines)
                )
        expected = ff.run_engine(ff.start(domain, oracle_features))
        assert produced.masks == expected.masks


# ------------------------------------------------ the region, in metres (R5)

WGS84 = "EPSG:4326"


class TestRegion:
    """R5 as revised on 2026-09-28: the domain is buffered by 100 m in the
    DEM's CRS, densified, moved into the source CRS, and the region is the
    convex hull. Pinned: `source_region(domain, dem_crs, source_crs)`, a
    shapely `Polygon` in the source CRS. The first design buffered in the
    source CRS, which for EPSG:4326 is 100 degrees."""

    @pytest.mark.parametrize("source_crs", [WGS84, LAEA, UTM33])
    def test_the_region_is_the_domain_plus_about_100_metres(
        self, fi: ModuleType, source_crs: str
    ) -> None:
        """Moved back into the DEM's CRS, the region holds the domain with
        at least 90 m to spare and reaches no more than 150 m from it (a
        mitred corner of a 100 m buffer is 141 m out)."""
        region = fi.source_region(BOX, UTM33, source_crs)
        assert isinstance(region, Polygon)
        back = moved(region, source_crs, UTM33)
        assert back.contains(BOX.polygon.buffer(90.0))
        assert BOX.polygon.buffer(150.0).contains(back)

    def test_a_wgs84_region_is_hundredths_of_a_degree_not_degrees(self, fi: ModuleType) -> None:
        """The case the first design got wrong: at 59.5° N, 100 m is about
        0.0009° of latitude and 0.0018° of longitude."""
        region = fi.source_region(BOX, UTM33, WGS84)
        image = moved(BOX.polygon, UTM33, WGS84).bounds
        grown = np.subtract(region.bounds, image) * [-1, -1, 1, 1]
        assert np.all(grown > 0.0)
        assert np.all(grown < 0.01)


# ------------------------------------ the pre-clip drops whole edges (R5)

#: A region in some source CRS; the pre-clip is pure geometry in that CRS.
REGION = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])


def chains(fi: ModuleType, geometry: Any) -> list[list[tuple[float, float]]]:
    return [[(float(x), float(y)) for x, y in c.coords] for c in fi.pre_clip(geometry, REGION)]


class TestPreClipKeepsWholeEdges:
    """R5 as revised on 2026-09-28. Pinned: `pre_clip(geometry, region)`, a
    tuple of `LineString`s in the source CRS: a polygon's exterior's chains,
    then each hole's; a closed ring repeats its first vertex. An edge is kept
    when the closed segment meets the closed region; each maximal run of
    edges that miss it is dropped; no vertex is added.

    Amended 2026-09-28 (R5, "Long edges"): `pre_clip(geometry, region,
    widening)`, where `widening(a, b)` takes an edge's two vertices as `(x,
    y)` tuples and gives `w(e)` in source units; an edge is kept when its
    distance to the region is `<= w(e)`. Left out, `w = 0`, which is the
    rule every other test in this class pins."""

    def test_a_ring_wholly_inside_stays_closed_with_its_start_vertex(self, fi: ModuleType) -> None:
        ring = Polygon([(6, 2), (8, 2), (8, 8), (2, 8)])
        assert chains(fi, ring) == [list(ring.exterior.coords)]

    def test_a_ring_poking_out_with_no_edge_wholly_outside_stays_closed(
        self, fi: ModuleType
    ) -> None:
        """Every edge meets the region (the long one only at its corner
        (10, 10)), so none is dropped, although two vertices lie outside. The
        exact clip in the DEM's CRS cuts it later."""
        ring = Polygon([(5, 5), (15, 5), (5, 15)])
        assert chains(fi, ring) == [list(ring.exterior.coords)]

    def test_a_ring_losing_one_run_is_one_open_chain_joined_across_its_start(
        self, fi: ModuleType
    ) -> None:
        """Edge (20, 5)-(20, 8) misses the region. What is left runs from its
        far end, through the ring's start (5, 5), to its near end."""
        ring = Polygon([(5, 5), (8, 5), (20, 5), (20, 8), (8, 8)])
        assert chains(fi, ring) == [[(20, 8), (8, 8), (5, 5), (8, 5), (20, 5)]]

    def test_a_ring_losing_two_runs_is_two_open_chains(self, fi: ModuleType) -> None:
        ring = Polygon([(5, 2), (20, 2), (20, 4), (5, 4), (5, 6), (20, 6), (20, 8), (2, 8)])
        got = sorted(chains(fi, ring))
        assert got == sorted(
            [
                [(20, 4), (5, 4), (5, 6), (20, 6)],
                [(20, 8), (2, 8), (5, 2), (20, 2)],
            ]
        )

    def test_a_run_of_several_edges_is_dropped_as_one(self, fi: ModuleType) -> None:
        ring = Polygon([(5, 5), (20, 5), (30, 5), (30, 8), (20, 8)])
        assert chains(fi, ring) == [[(20, 8), (5, 5), (20, 5)]]

    def test_a_line_losing_its_middle_is_two_chains(self, fi: ModuleType) -> None:
        line = LineString([(1, 5), (8, 5), (20, 5), (20, 7), (8, 7), (2, 7)])
        assert sorted(chains(fi, line)) == sorted(
            [[(1, 5), (8, 5), (20, 5)], [(20, 7), (8, 7), (2, 7)]]
        )

    def test_a_line_losing_its_ends_keeps_its_middle(self, fi: ModuleType) -> None:
        line = LineString([(-20, 5), (-10, 5), (5, 5), (20, 5), (30, 5)])
        assert chains(fi, line) == [[(-10, 5), (5, 5), (20, 5)]]

    def test_an_edge_touching_the_region_at_one_point_is_kept(self, fi: ModuleType) -> None:
        """Closed segment, closed region: an edge through the corner (10, 10)
        and nowhere else in the region is kept; the two edges beside it miss."""
        touching = LineString([(20, 20), (15, 5), (5, 15), (-5, 25)])
        assert chains(fi, touching) == [[(15, 5), (5, 15)]]

    def test_a_polygon_whose_exterior_is_dropped_keeps_its_hole_closed(
        self, fi: ModuleType
    ) -> None:
        hole = [(3, 3), (7, 3), (7, 7), (3, 7), (3, 3)]
        polygon = Polygon([(-50, -50), (50, -50), (50, 50), (-50, 50)], [hole])
        assert chains(fi, polygon) == [hole]

    def test_parts_in_order_exterior_before_holes(self, fi: ModuleType) -> None:
        a = Polygon([(1, 1), (4, 1), (4, 4)], [[(2, 1.5), (3.5, 1.5), (3.5, 3)]])
        b = Polygon([(6, 6), (9, 6), (9, 9)])
        got = chains(fi, MultiPolygon([a, b]))
        assert got == [
            list(a.exterior.coords),
            list(a.interiors[0].coords),
            list(b.exterior.coords),
        ]

    @pytest.mark.parametrize(
        "geometry",
        [
            pytest.param(Polygon([(5, 5), (8, 5), (20, 5), (20, 8), (8, 8)]), id="one-run"),
            pytest.param(
                Polygon([(5, 2), (20, 2), (20, 4), (5, 4), (5, 6), (20, 6), (20, 8), (2, 8)]),
                id="two-runs",
            ),
            pytest.param(LineString([(1, 5), (8, 5), (20, 5), (20, 7), (8, 7)]), id="line"),
            pytest.param(Polygon([(-5, 3), (15, 3.3), (15, 7.1), (-5, 6.9)]), id="crossing"),
            pytest.param(
                LineString([(-7.3, -2.1), (13.9, 12.7), (25.0, 3.3), (4.4, -9.1)]),
                id="slanted",
            ),
        ],
    )
    def test_no_vertex_is_added(self, fi: ModuleType, geometry: Any) -> None:
        """Intersecting with the region would add a vertex on its boundary
        wherever an edge crosses it; every case here has such crossings."""
        source = {tuple(p) for p in shapely.get_coordinates(geometry).tolist()}
        cut = shapely.intersection(shapely.MultiLineString(ff.boundary_lines(geometry)), REGION)
        assert {tuple(p) for p in shapely.get_coordinates(cut).tolist()} - source
        got = {p for c in chains(fi, geometry) for p in c}
        assert got and got <= source

    def test_nothing_kept_is_nothing(self, fi: ModuleType) -> None:
        around = Polygon([(-50, -50), (50, -50), (50, 50), (-50, 50)])
        assert fi.pre_clip(around, REGION) == ()

    def test_an_edge_within_its_widening_is_kept_and_one_beyond_it_is_not(
        self, fi: ModuleType
    ) -> None:
        """Edge (12, 5)-(12, 20) lies exactly 2 from the region. It is kept
        when its own widening is 2 (the distance test is closed) and dropped
        when it is 1.9. The widening is per edge: the edge (12, 20)-(30, 20),
        about 10.2 from the region, gets 0, so it is dropped either way and
        the kept chain ends at (12, 20)."""
        line = LineString([(5, 5), (12, 5), (12, 20), (30, 20)])
        tall = {(12.0, 5.0), (12.0, 20.0)}

        def widening(w: float) -> Any:
            def of(a: tuple[float, float], b: tuple[float, float]) -> float:
                return w if {tuple(map(float, a)), tuple(map(float, b))} == tall else 0.0

            return of

        kept = fi.pre_clip(line, REGION, widening(2.0))
        assert [list(c.coords) for c in kept] == [[(5, 5), (12, 5), (12, 20)]]
        dropped = fi.pre_clip(line, REGION, widening(1.9))
        assert [list(c.coords) for c in dropped] == [[(5, 5), (12, 5)]]


class TestPreClipThroughOpenFeatures:
    def test_a_polygon_around_the_domain_with_no_edge_kept_is_dropped_and_counted(
        self, tmp_path: Path
    ) -> None:
        """It meets the region (it holds it), so it is read; none of its edges
        do, so the pre-clip keeps nothing and it counts as outside."""
        fs = one(tmp_path, square(-50_000, -50_000, 50_000, 50_000))
        assert fs.features == () and fs.outside == 1


# ------------------------------- a 10 km WGS84 edge across the region (I4)

#: A 10 km domain, so that a 10 km edge can cross the region's boundary near
#: its middle, where its bend is largest, and still enter the domain.
WIDE = domain_of(square(0, 0, 10_000, 10_000))


def wgs84_east_west_edge() -> LineString:
    """An edge along one parallel, straight in EPSG:4326, from 5.1 km west of
    the domain to 4.9 km inside it: about 10 km, its middle near the region's
    west boundary (100 m west of the domain)."""
    to_wgs84 = Transformer.from_crs(UTM33, WGS84, always_xy=True)
    lon_a, lat = to_wgs84.transform(*at(-5_100, 5_000))
    lon_b, _ = to_wgs84.transform(*at(4_900, 5_000))
    return LineString([(lon_a, lat), (lon_b, lat)])


def clipped_in_the_dem(line: LineString) -> tuple[LineString, ...]:
    """Moved vertex by vertex into the DEM's CRS, clipped to `WIDE` there."""
    pieces = shapely.get_parts(shapely.intersection(moved(line, WGS84), WIDE.polygon))
    return tuple(p for p in pieces if isinstance(p, LineString) and p.length > 0)


class TestLongWgs84Edge:
    def test_an_edge_crossing_the_region_boundary_is_not_split_there(self, tmp_path: Path) -> None:
        """I4 for an EPSG:4326 source (R5, revised 2026-09-28). The produced
        noded edges equal those of this test's own pipeline with no pre-clip:
        the whole edge moved vertex by vertex and clipped in the DEM's CRS.

        First, that this can fail: a pre-clip that intersects the edge with
        the region (built here with shapely and pyproj as R5 describes it)
        splits it on the region's boundary, and its piece in the domain meets
        the domain's west side metres away from the unsplit edge's."""
        edge = wgs84_east_west_edge()
        (whole,) = clipped_in_the_dem(edge)
        region = shapely.convex_hull(
            moved(WIDE.polygon.buffer(100.0).segmentize(1_000.0), UTM33, WGS84)
        )
        (split,) = clipped_in_the_dem(shapely.intersection(edge, region))
        west = np.array([at(0, 0), at(0, 10_000)])
        hit_whole = shapely.intersection(whole, LineString(west))
        hit_split = shapely.intersection(split, LineString(west))
        assert hit_whole.distance(hit_split) > 1.0

        path = write_geojson(
            tmp_path / "f.geojson", [Feat("edge", edge, {"property": "road"})], crs=None
        )
        fs = open_one(path, WIDE)
        produced = ff.run_engine(ff.start(WIDE, fs.features))
        fi = ff.feature_input()
        oracle = fi.TerrainFeature(fid="edge", mask=V.mask("road"), lines=(whole,))
        expected = ff.run_engine(ff.start(WIDE, [oracle]))
        assert produced.masks == expected.masks


# -------------------- a 60 km WGS84 edge that misses the region (R5, long edges)

#: The edge's parallel, and its half-width in longitude: about 60 km at 70° N,
#: centred on UTM 33's central meridian, where its bend is largest.
LONG_LAT, LONG_HALF = 70.0, 0.788


def far_north_domain() -> tuple[LineString, Any]:
    """The edge in EPSG:4326, and a domain 2 km square in UTM 33 centred over
    the edge's midpoint, its southern side 140 m north of that midpoint."""
    edge = LineString([(15.0 - LONG_HALF, LONG_LAT), (15.0 + LONG_HALF, LONG_LAT)])
    mx, my = Transformer.from_crs(WGS84, UTM33, always_xy=True).transform(15.0, LONG_LAT)
    west, east, south, north = mx - 1_000, mx + 1_000, my + 140, my + 2_140
    return edge, domain_of(Polygon([(west, south), (east, south), (east, north), (west, north)]))


class TestLongEdgeWidening:
    def test_an_edge_outside_the_region_whose_dem_straight_line_enters_is_kept(
        self, tmp_path: Path
    ) -> None:
        """R5, "Long edges" (Ola, 2026-09-28: "yes, widen per edge"), and I4
        with no length limit. The edge is about 40 m outside the region in
        EPSG:4326, but the engine draws it straight in UTM 33, where the
        chord between its two images runs about 194 m north of its midpoint,
        through the domain. It is kept, and the produced noded edges equal
        those of this test's own pipeline with no pre-clip: the edge moved
        vertex by vertex and clipped to the domain in UTM 33.

        First, that this can fail: the edge misses the region (the domain
        buffered by 100 m in UTM 33, densified, moved into EPSG:4326, its
        convex hull, as R5 builds it), so a pre-clip with no widening drops
        it; yet its UTM-straight line crosses the domain."""
        edge, domain = far_north_domain()
        assert geodesic_length(edge) == pytest.approx(60_000.0, rel=0.01)
        region = shapely.convex_hull(
            moved(domain.polygon.buffer(100.0).segmentize(1_000.0), UTM33, WGS84)
        )
        assert not edge.intersects(region)
        assert edge.distance(region) == pytest.approx(0.00036, abs=0.00005)  # about 40 m
        dem_straight = moved(edge, WGS84)
        (whole,) = (
            p
            for p in shapely.get_parts(shapely.intersection(dem_straight, domain.polygon))
            if isinstance(p, LineString) and p.length > 0
        )
        mid_y = Transformer.from_crs(WGS84, UTM33, always_xy=True).transform(15.0, LONG_LAT)[1]
        assert whole.coords[0][1] - mid_y == pytest.approx(194.0, abs=2.0)

        path = write_geojson(
            tmp_path / "f.geojson", [Feat("edge", edge, {"property": "road"})], crs=None
        )
        fs = open_one(path, domain)
        assert [f.fid for f in fs.features] == ["edge"] and fs.outside == 0
        produced = ff.run_engine(ff.start(domain, fs.features))
        fi = ff.feature_input()
        oracle = fi.TerrainFeature(fid="edge", mask=V.mask("road"), lines=(whole,))
        expected = ff.run_engine(ff.start(domain, [oracle]))
        assert produced.masks == expected.masks


def geodesic_length(line: LineString) -> float:
    lon, lat = zip(*line.coords, strict=True)
    return float(Geod(ellps="WGS84").line_length(lon, lat))


# ---------------------------------------------- several sources (16e, R6/D2)


class TestManySources:
    """`docs/increments/16e-multi-features.md` R6/D2: `open_features` on a
    request of several sources returns per-source tuples. `crs` and `layers`
    are already per-source (this pins them against regression); `counts` is the
    new field D2 adds — features kept per source, in source order."""

    def two_source_request(self, tmp_path: Path, fi: ModuleType) -> Any:
        """A first source of two features inside `INNER`, a second of one."""
        big = square(120, 120, 180, 180)
        small = square(130, 130, 150, 150)
        one_more = square(160, 160, 175, 175)
        a = write_geojson(
            tmp_path / "a.geojson",
            [
                Feat("big", big, {"property": "land_cover"}),
                Feat("small", small, {"property": "water"}),
            ],
        )
        b = write_geojson(tmp_path / "b.geojson", [Feat("one", one_more, {"property": "wall"})])
        return fi.FeatureRequest(sources=(ff.source(a), ff.source(b)))

    def test_counts_are_per_source_in_order(self, tmp_path: Path, fi: ModuleType) -> None:
        fs = fi.open_features(self.two_source_request(tmp_path, fi), BOX, UTM33)
        assert fs.counts == (2, 1)
        assert sum(fs.counts) == len(fs.features)

    def test_crs_and_layers_are_per_source(self, tmp_path: Path, fi: ModuleType) -> None:
        fs = fi.open_features(self.two_source_request(tmp_path, fi), BOX, UTM33)
        assert fs.crs == (UTM33, UTM33)
        assert fs.layers == (None, None)
