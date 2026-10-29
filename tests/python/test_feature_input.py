"""`tin_engine.feature_input`: sources, legends, CRS, region and clip (16b-1).

`docs/increments/16b-terrain-polygons.md` R1, R4-R7, "Degeneracy policy" and
"Tests for @tester" (legends, CRS, region and pre-clip, clip degeneracies),
with Ola's Q6 (b) (the legacy GML, read with the `clc18_kode` map).

Pinned by this suite (`docs/increments/111-one-land-cover-map.md`, section 8):

- A class-polygon file is read by `read_class_polygons` and a line file by
  `read_linestrings`; `feature_fixtures.open_one` and `open_lines` run them
  through `land_cover_map.build` and `edges_of` as `rasputin mesh` does.
  The shipped land-cover class files `corine` and `clc18_kode` (`corine` keyed
  on the Norwegian GML's attribute) carry the Copernicus attribution.
- `edges_of(...).features` is a tuple of `TerrainFeature` (`fid`, `mask`,
  `lines`): the areas first, `fid` the index in the map's areas, then the
  lines, `fid` `"<key> <i>"`. `outside` counts features wholly clipped away;
  the readers' `empty` counts empty rows skipped.
- A GeoJSON feature's `fid`, named in a refusal, is its `id` member, else its
  position in the file (from 0). A GeoPackage row's is its primary key; a GML
  feature's its `fid` attribute.
- By suffix: `.geojson` and `.json`, `.gpkg`, `.gml`; anything else refused.
- `FeatureError`, a `ValueError`, for every refusal; it names the file, and
  for a feature its `fid`.
- The lines of a polygon feature: its exterior's pieces, then each hole's, part
  by part for a `MultiPolygon`.

Committed red at `e99c8ea` (amended at `3990449` and `972312c`):
`tin_engine.feature_input` did not exist yet, and because it is imported lazily
each test failed on its own. It landed in `5079da8` and the suite has been
green since.

Increment 30b (`docs/increments/30b-clip-speed.md`, section 7) added the
classes at the end: `TestUntouched` (R1, R2) for `_untouched`, committed red
before it existed; and `TestTheClipStaysTheSame` (pins P1, P2, P3, P5),
green before the change.
"""

from __future__ import annotations

import itertools
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
from shapely.geometry.base import BaseGeometry

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
    open_lines,
    open_one,
    square,
    within,
    write_geojson,
)
from gpkg_fixtures import EXTRACT, LEGACY_GML, Layer, Row, needs_rtree, write_gpkg
from test_cli_mesh_domain import quarter_circle
from test_io_gml import document, member
from tin_engine.features import DEFAULT_VOCABULARY

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


LINEAR = ("LineString", "MultiLineString")


def one(tmp_path: Path, geometry: Any, props: dict[str, Any] | None = None, **kw: Any) -> Any:
    """One feature: a polygon as a CORINE class polygon (forest 311 unless
    ``props`` says), a line as a `--lines` line of the key ``kw["key"]``
    (`road` unless given)."""
    domain = kw.pop("domain", BOX)
    key = kw.pop("key", "road")
    if isinstance(geometry, BaseGeometry) and geometry.geom_type in LINEAR:
        path = write_geojson(tmp_path / "f.geojson", [Feat("a", geometry, props or {})])
        return open_lines(path, domain, key, **kw)
    path = write_geojson(tmp_path / "f.geojson", [Feat("a", geometry, props or {"Code_18": "311"})])
    return open_one(path, domain, kw.pop("map_name", "corine"), **kw)


def lines_of(fs: Any) -> list[LineString]:
    return [line for f in fs.features for line in f.lines]


# --------------------------------------------------------------- legends


BUILT_IN = ("corine", "clc18_kode")


class TestLegends:
    def test_the_built_in_maps(self) -> None:
        maps = {name: ff.classes(name) for name in BUILT_IN}
        assert {name: m.attribute for name, m in maps.items()} == {
            "corine": "Code_18",
            "clc18_kode": "clc18_kode",
        }
        assert all(m.name == name for name, m in maps.items())

    def test_the_corine_maps_carry_the_attribution(self, fi: ModuleType) -> None:
        """R10 and "Test data": a mesh built from CORINE says so, and says the
        data were modified."""
        for name in BUILT_IN:
            notice = ff.classes(name).notice
            assert "Copernicus Land Monitoring Service" in notice, name
            assert "modified" in notice.lower(), name

    def test_a_legend_is_frozen_and_strict(self, fi: ModuleType) -> None:
        corine = ff.classes("corine")
        with pytest.raises((TypeError, ValueError)):
            corine.attribute = "other"
        with pytest.raises(ValueError):
            fi.LandCoverClasses(
                name="x",
                attribute="a",
                codes="c",
                legend=corine.legend,
                classes=frozenset(),
                colour="red",
            )

    @pytest.mark.parametrize("code", WATER_CODES)
    def test_corine_water_codes_are_land_cover_and_water(self, tmp_path: Path, code: str) -> None:
        """MUTANT 3 (in `test_feature_chains_bits`), unit by unit."""
        fs = one(tmp_path, INNER, {"Code_18": code}, map_name="corine")
        assert [f.mask for f in fs.features] == [V.mask("land_cover", "water")]

    @pytest.mark.parametrize("code", ["111", "311", "412", "423", "999", "5", "51", "513", "530"])
    def test_corine_every_other_code_is_land_cover(self, tmp_path: Path, code: str) -> None:
        """`unlisted = "keep"`: an unlisted code and the 5xx codes the
        nomenclature does not have included."""
        fs = one(tmp_path, INNER, {"Code_18": code}, map_name="corine")
        assert [f.mask for f in fs.features] == [V.mask("land_cover")]
        assert [f.code for f in fs.features] == [int(code)]

    @pytest.mark.parametrize(
        ("code", "names"), [("512", ("land_cover", "water")), ("322", ("land_cover",))]
    )
    def test_clc18_kode(self, tmp_path: Path, code: str, names: tuple[str, ...]) -> None:
        fs = one(tmp_path, INNER, {"clc18_kode": code}, map_name="clc18_kode")
        assert [f.mask for f in fs.features] == [V.mask(*names)]

    def test_only_the_chosen_classes_give_lines(self, tmp_path: Path, fi: ModuleType) -> None:
        """`EdgeChoice(codes=...)`: a forest beside a lake, the lake's border
        chosen; the forest gives no feature."""
        features = [
            Feat(1, INNER, {"Code_18": "311"}),
            Feat(2, square(10, 10, 50, 50), {"Code_18": "512"}),
        ]
        path = write_geojson(tmp_path / "f.geojson", features)
        fs = open_one(path, BOX, choice=fi.EdgeChoice(codes=frozenset({512})))
        assert [(f.fid, f.code, f.mask) for f in fs.features] == [
            (1, 512, V.mask("land_cover", "water"))
        ]


# ---------------------------------------------------------------- sources


class TestSources:
    def test_geojson_order_and_ids(self, tmp_path: Path, fi: ModuleType) -> None:
        """File order is kept; an area's fid is its index in the map. The
        source's own id (its `id` member, else its position) is what a
        refusal names."""
        features = [
            Feat("z", square(10, 10, 20, 20), {"Code_18": "311"}),
            Feat(None, square(30, 10, 40, 20), {"Code_18": "512"}),
            Feat(7, square(50, 10, 60, 20), {"Code_18": "312"}),
        ]
        fs = open_one(write_geojson(tmp_path / "f.geojson", features), BOX)
        assert [(f.fid, f.code) for f in fs.features] == [(0, 311), (1, 512), (2, 312)]
        ok = Feat("ok", square(30, 10, 40, 20), {"Code_18": "311"})
        for fid, named in (("z", "feature z"), (None, "feature 1"), (7, "feature 7")):
            bad = Feat(fid, square(10, 10, 20, 20), {"Code_18": "x"})
            path = write_geojson(tmp_path / "bad.geojson", [ok, bad])
            with pytest.raises(fi.FeatureError, match=rf"bad\.geojson: {named}\b"):
                open_one(path, BOX)

    def test_json_suffix_is_geojson(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.json", [Feat("a", INNER, {"Code_18": "311"})])
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
        assert [f.code for f in fs.features] == [311, 512, 311]
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
        assert [f.code for f in fs.features] == [512]
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
        them: a sea class (523) is among them."""
        domain = domain_of(Polygon(quarter_circle()))
        fs = open_one(LEGACY_GML, domain, "clc18_kode")
        assert 523 in {f.code for f in fs.features}
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
            pytest.param(LineString([at(10, 10), at(90, 90)]), id="line"),
        ],
    )
    def test_loose_points_collections_and_lines_are_refused_as_class_polygons(
        self, tmp_path: Path, fi: ModuleType, geometry: Any
    ) -> None:
        """A line goes to `--lines`, not `--land-cover`."""
        path = write_geojson(tmp_path / "f.geojson", [Feat("pt-3", geometry, {"Code_18": "311"})])
        with pytest.raises(fi.FeatureError, match="pt-3"):
            open_one(path, BOX)

    @pytest.mark.parametrize(
        "geometry",
        [
            pytest.param(Point(*at(150, 150)), id="point"),
            pytest.param(GeometryCollection([INNER]), id="geometrycollection"),
            pytest.param(INNER, id="polygon"),
            pytest.param(MultiPolygon([INNER]), id="multipolygon"),
        ],
    )
    def test_points_collections_and_polygons_are_refused_as_lines(
        self, tmp_path: Path, fi: ModuleType, geometry: Any
    ) -> None:
        """A polygon goes to `--prefer`; the message names it."""
        path = write_geojson(tmp_path / "f.geojson", [Feat("pt-4", geometry, {})])
        with pytest.raises(fi.FeatureError, match="pt-4") as caught:
            open_lines(path, BOX)
        if geometry.geom_type in ("Polygon", "MultiPolygon"):
            assert "--prefer" in str(caught.value)

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
            Feat("e", geometry, {"Code_18": "311"}),
            Feat("f", INNER, {"Code_18": "512"}),
        ]
        fs = open_one(write_geojson(tmp_path / "f.geojson", features), BOX)
        assert [f.code for f in fs.features] == [512]
        assert fs.empty == 1

    @pytest.mark.parametrize(
        "geometry",
        [
            pytest.param({"type": "LineString", "coordinates": []}, id="empty-line"),
            pytest.param(None, id="null"),
        ],
    )
    def test_an_empty_line_row_is_skipped_and_counted(self, tmp_path: Path, geometry: Any) -> None:
        line = LineString([at(10, 10), at(90, 90)])
        features = [Feat("e", geometry, {}), Feat("f", line, {})]
        fs = open_lines(write_geojson(tmp_path / "f.geojson", features), BOX)
        assert [f.fid for f in fs.features] == ["road 0"]
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
        assert [f.fid for f in fs.features] == [0]
        assert fs.empty == 1

    def test_a_self_intersecting_ring_is_linework_not_refused(self, tmp_path: Path) -> None:
        """R1: class polygons are read, not validated nor refused. Its source vertices are all
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
            tmp_path / "f.geojson", [Feat("a", lonlat, {"Code_18": "311"})], crs=None
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
        leaving = square(150, 150, 400, 170)  # clipped at the domain
        features = [Feat(fid, g, {"Code_18": "311"}) for fid, g in (("a", INNER), ("b", leaving))]
        by_code = open_one(write_geojson(tmp_path / "code.geojson", features), BOX)
        spelt = write_geojson(tmp_path / "proj.geojson", features, crs=axes_swapped(25833))
        refuse_point_moves(monkeypatch)
        by_proj = open_one(spelt, BOX)
        assert [f.fid for f in by_proj.features] == [f.fid for f in by_code.features]
        assert [line.coords[:] for line in lines_of(by_proj)] == [
            line.coords[:] for line in lines_of(by_code)
        ]

    def test_a_given_crs_agreeing_with_the_files_is_accepted(self, tmp_path: Path) -> None:
        path = write_geojson(tmp_path / "f.geojson", [Feat("a", INNER, {"Code_18": "311"})])
        assert len(open_one(path, BOX, crs="EPSG:25833").features) == 1
        assert len(open_one(path, BOX, crs=CRS.from_epsg(25833).to_wkt()).features) == 1

    def test_a_given_crs_disagreeing_with_the_files_is_refused(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        path = write_geojson(tmp_path / "f.geojson", [Feat("a", INNER, {"Code_18": "311"})])
        with pytest.raises(fi.FeatureError, match=r"f\.geojson"):
            open_one(path, BOX, crs="EPSG:3035")

    def test_a_given_crs_disagreeing_with_rfc_7946_is_refused(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        """domain.py's rule exactly: no `crs` member means EPSG:4326."""
        path = write_geojson(
            tmp_path / "f.geojson", [Feat("a", INNER, {"Code_18": "311"})], crs=None
        )
        with pytest.raises(fi.FeatureError):
            open_one(path, BOX, crs="EPSG:25833")

    def test_an_unreadable_crs_is_refused(self, tmp_path: Path, fi: ModuleType) -> None:
        path = write_geojson(
            tmp_path / "f.geojson", [Feat("a", INNER, {"Code_18": "311"})], crs="EPSG:0"
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

        def fake(src: Any, dst: Any, *, threads: int = 1) -> Any:
            move = real(src, dst, threads=threads)
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
            tmp_path / "f.geojson", [Feat("lost", laea, {"Code_18": "311"})], crs=LAEA
        )
        with pytest.raises(fi.FeatureError, match="lost"):
            open_one(path, BOX)


# --------------------------------------------------------- region and clip


class TestOneGeojsonRule:
    """Audit PR C (`docs/increments/python-audit.md`, section 12, red test 2):
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
        """The same file through `read_class_polygons`: its `crs` is the text."""
        lonlat = moved(INNER, UTM33, "EPSG:4326")
        feature = {"type": "Feature", "properties": {"Code_18": "311"}}
        doc = {
            "type": "FeatureCollection",
            "crs": {"type": "name", "properties": {"name": 4326}},
            "features": [feature | {"geometry": mapping(lonlat)}],
        }
        fs = open_one(self.written(tmp_path, doc), BOX)
        assert fs.crs == "4326"
        assert len(fs.features) == 1


#: Every empty JSON value but `null` (`docs/increments/python-audit.md`,
#: section 12, the ruling after code review round 2).
EMPTY_NOT_NULL = [
    pytest.param("", id="empty_string"),
    pytest.param(0, id="zero"),
    pytest.param(False, id="false"),
    pytest.param([], id="empty_list"),
    pytest.param({}, id="empty_object"),
]


class TestAnEmptyGeometryIsNull:
    """Audit PR C, section 12's ruling after code review round 2: `read_source`
    treats every empty `geometry` as `null`. Before, `""`, `0`, `false`, `[]`
    and `{}` crashed `--features` and `catchment --lakes` with an
    `AttributeError` (`'str' object has no attribute 'is_empty'`). Each test
    compares with the same file whose `geometry` is `null`, so it pins no
    outcome of its own: only that the two read alike."""

    @staticmethod
    def written(tmp_path: Path, geometry: Any, shape: str) -> Path:
        """A top-level `Feature`, or a collection whose first feature has
        `geometry` and whose second is `INNER`; both with a `crs` member."""
        feature = {"type": "Feature", "properties": {"Code_18": "512"}, "geometry": geometry}
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
        read = fi.read_source(path, None, "Code_18", lambda _crs: (0.0, 0.0, 0.0, 0.0))
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
        """`--land-cover`, through `read_class_polygons` and `edges_of`."""
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
        (line,) = lines_of(one(tmp_path, wall, key="wall"))
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
        """I4 on the extract, relationally: `open_one` + `start_chains`
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


class TestAPolygonAroundTheDomain:
    def test_a_class_polygon_around_the_domain_labels_it_with_no_line(self, tmp_path: Path) -> None:
        """No edge of it meets the domain, and a point of it is inside: it
        covers the domain, so it is kept with no line, for its class code,
        and not counted outside (as a coded polygon is today)."""
        fs = one(tmp_path, square(-50_000, -50_000, 50_000, 50_000))
        assert [(f.code, f.lines) for f in fs.features] == [(311, ())]
        assert fs.outside == 0


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

        path = write_geojson(tmp_path / "f.geojson", [Feat("edge", edge, {})], crs=None)
        fs = open_lines(path, WIDE)
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

        path = write_geojson(tmp_path / "f.geojson", [Feat("edge", edge, {})], crs=None)
        fs = open_lines(path, domain)
        assert [f.fid for f in fs.features] == ["road 0"] and fs.outside == 0
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
    """`edges_of`'s `counts`: one entry per source edges were taken from,
    the map first (when the choice names a class), then each chosen
    `--lines` key in the order its first line appears."""

    def opened(self, tmp_path: Path, fi: ModuleType) -> Any:
        """A map of two polygons inside `INNER`; a `wall` line, then two
        `road` lines."""
        big = square(120, 120, 180, 180)
        small = square(130, 130, 150, 150)
        cover = write_geojson(
            tmp_path / "a.geojson",
            [Feat("big", big, {"Code_18": "311"}), Feat("small", small, {"Code_18": "512"})],
        )
        known = ff.classes("corine")
        region = ff.region_of(BOX)
        read = fi.read_class_polygons(cover, None, known, UTM33, ff.box_of(BOX))
        built = ff.land_cover_map().build([read.areas], [], known, region)
        lines = [
            ("wall", LineString([at(160, 160), at(175, 175)])),
            ("road", LineString([at(10, 10), at(90, 10)])),
            ("road", LineString([at(10, 20), at(90, 20)])),
        ]
        choice = fi.EdgeChoice(every_class=True, lines=frozenset({"road", "wall"}))
        return fi.edges_of(built, lines, choice, ff.cleanup(), BOX, V, region=region), read

    def test_counts_are_per_source_in_order(self, tmp_path: Path, fi: ModuleType) -> None:
        fs, _ = self.opened(tmp_path, fi)
        assert fs.counts == (2, 1, 2)
        assert sum(fs.counts) == len(fs.features)
        assert [f.fid for f in fs.features] == [0, 1, "wall 0", "road 1", "road 2"]

    def test_the_reader_gives_the_crs_and_layer(self, tmp_path: Path, fi: ModuleType) -> None:
        _, read = self.opened(tmp_path, fi)
        assert (read.crs, read.layer, read.scanned, read.empty) == (UTM33, None, False, 0)


# ------------------------------------ 30b: the features clip made fast

#: R1 and R2's offset: `docs/increments/30b-probes/geos_cases.py`'s, whole
#: metres, so every coordinate below is exactly representable.
U0 = (500_000.0, 6_600_000.0)


def u(x: float, y: float) -> tuple[float, float]:
    return (U0[0] + x, U0[1] + y)


#: A square with an extra vertex at (100, 50) on its east side, and a square
#: with a hole (40..60)²: `geos_cases.py`'s two domains.
EXTRA = Polygon([u(0, 0), u(100, 0), u(100, 50), u(100, 100), u(0, 100)])
HOLED = Polygon(
    [u(0, 0), u(100, 0), u(100, 100), u(0, 100)],
    [[u(40, 40), u(60, 40), u(60, 60), u(40, 60)]],
)

#: R1: `geos_cases.py`'s 14 cases, then two the domain does not cover, then
#: increment 35's five (section 3.2, 5 test 3 and its review): lines whose
#: ends lie on the outline or a hole. `True` exactly where
#: `shapely.intersection(line, domain)` gives the line back unchanged
#: (shapely 2.1.2, GEOS 3.13.1, 30b; shapely 2.2.0, GEOS 3.14.1, 35).
UNTOUCHED_CASES = [
    pytest.param(EXTRA, [u(10, 10), u(20, 30), u(50, 20)], True, id="plain-inside"),
    pytest.param(EXTRA, [u(10, 10), u(30, 10), u(30, 30), u(10, 30), u(10, 10)], True,
                 id="closed-ring-inside"),
    pytest.param(EXTRA, [u(10, 10), u(20, 10), u(30, 10)], True, id="collinear-middle-vertex"),
    pytest.param(EXTRA, [u(10, 10), u(20, 30), u(20, 30), u(50, 20)], False,
                 id="repeated-vertex"),
    pytest.param(EXTRA, [u(0, 10), u(30, 10), u(30, 30), u(0, 30), u(0, 10)], False,
                 id="closed-ring-touching-the-outline-at-its-start"),
    pytest.param(EXTRA, [u(10, 10), u(30, 10), u(30, 30), u(0, 20), u(10, 10)], False,
                 id="closed-ring-touching-the-outline-mid-ring"),
    pytest.param(EXTRA, [u(50, 50), u(100, 40), u(100, 60), u(60, 70)], False,
                 id="domain-vertex-on-the-line"),
    pytest.param(EXTRA, [u(100, 10), u(100, 90)], False, id="along-the-outline"),
    pytest.param(EXTRA, [u(100, 10), u(100, 90), u(50, 90)], False, id="along-then-inside"),
    pytest.param(EXTRA, [u(10, 10), u(50, 50), u(50, 10), u(10, 50)], False, id="self-crossing"),
    pytest.param(EXTRA, [u(10, 10), u(50, 10), u(30, 30), u(10, 10), u(10, 50)], False,
                 id="self-touching-at-a-vertex"),
    pytest.param(EXTRA, [u(10, 10), u(50, 10), u(30, 10)], False, id="doubling-back"),
    pytest.param(HOLED, [u(10, 10), u(40, 50), u(10, 90)], False, id="touching-a-hole"),
    pytest.param(EXTRA, [u(10, 10), u(10, 10)], False, id="zero-length"),
    pytest.param(EXTRA, [u(50, 50), u(150, 50)], False, id="partly-outside"),
    pytest.param(EXTRA, [u(200, 10), u(250, 30)], False, id="wholly-outside"),
    # Increment 35: an end may lie on the outline; the interior may not touch it.
    pytest.param(EXTRA, [u(0, 10), u(30, 30), u(100, 20)], True, id="ends-on-the-outline"),
    pytest.param(EXTRA, [u(10, 10), u(30, 30), u(100, 20)], True, id="one-end-on-the-outline"),
    pytest.param(HOLED, [u(10, 10), u(30, 30), u(40, 50)], True, id="end-on-a-hole"),
    pytest.param(EXTRA, [u(10, 10), u(30, 30), u(100, 50)], True, id="end-on-a-domain-vertex"),
    pytest.param(EXTRA, [u(0, 10), u(100, 50), u(0, 90)], False,
                 id="ends-on-the-outline-through-a-domain-vertex"),
]  # fmt: skip


class TestUntouched:
    """30b, section 7, R1 and R2, widened by increment 35, section 3.2:
    `_untouched(line, domain) -> bool`, True exactly when the line has two
    vertices or more, repeats no vertex in a row, is simple, and its interior
    lies in the domain's interior with its ends in the domain's closure
    (`relate_pattern(domain, line, "T**F**FF*")`; 30b asked
    `contains_properly`, which refused an end on the outline). Where it is
    True, the clip keeps the chain as it is instead of intersecting it with
    the domain, so it must be True only where that intersection would give the
    chain back unchanged. Called, as `_take` calls it, with the domain
    prepared."""

    @pytest.mark.parametrize(("domain", "coords", "expected"), UNTOUCHED_CASES)
    def test_r1_case_by_case(
        self, fi: ModuleType, domain: Polygon, coords: list[tuple[float, float]], expected: bool
    ) -> None:
        shapely.prepare(domain)
        line = LineString(coords)
        before = (shapely.to_wkb(line), shapely.to_wkb(domain))
        assert fi._untouched(line, domain) is expected
        assert (shapely.to_wkb(line), shapely.to_wkb(domain)) == before  # pure (3.4)

    @pytest.mark.parametrize(
        ("domain", "coords"),
        [pytest.param(*p.values[:2], id=p.id) for p in UNTOUCHED_CASES if p.values[2]],
    )
    def test_r2_every_true_row_is_what_geos_gives_back(
        self, fi: ModuleType, domain: Polygon, coords: list[tuple[float, float]]
    ) -> None:
        """Increment 35, test 3's oracle, on every True row of the table (the
        35 rows included): GEOS's own `intersection(line, domain)` is one
        `LineString` with the line's bytes. This is what makes skipping the
        intersection safe, and it runs on every GEOS CI has (section 8)."""
        shapely.prepare(domain)
        line = LineString(coords)
        assert fi._untouched(line, domain)
        parts = shapely.get_parts(shapely.intersection(line, domain))
        assert len(parts) == 1 and isinstance(parts[0], LineString), [p.wkt for p in parts]
        assert shapely.to_wkb(parts[0]) == shapely.to_wkb(line)

    def test_r2_where_it_is_true_geos_gives_the_line_back_unchanged(self, fi: ModuleType) -> None:
        """The oracle is GEOS's own intersection, so this is the CI guard
        against a GEOS that reorients, restarts or re-vertexes a line lying
        in the interior (section 10). 400 seeded lines of 2 to 5 vertices
        near a 48-vertex domain (integer vertices) with a square hole; by
        index, a fifth each: plain, rounded to a 10 m lattice, a closed
        triangle (either orientation), one vertex repeated, one vertex moved
        onto a domain vertex. At this seed 91 are untouched: 55 open, 19
        counter-clockwise rings, 17 clockwise. The bounds below are a tenth
        of the set and one of each shape, so it cannot pass vacuously."""
        t = np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False)
        r = 400.0 + 60.0 * np.sin(5.0 * t)
        outer = np.round(np.column_stack([r * np.cos(t), r * np.sin(t)])) + U0
        hole = np.array([(-60.0, -60.0), (-60.0, 60.0), (60.0, 60.0), (60.0, -60.0)]) + U0
        domain = Polygon(outer, [hole])
        assert domain.is_valid
        shapely.prepare(domain)
        rng = np.random.default_rng(302)
        seen = {"open": 0, "ccw": 0, "cw": 0}
        for i in range(400):
            kind = i % 5
            n = 3 if kind == 2 else int(rng.integers(2, 6))
            xy = rng.uniform(-470.0, 470.0, 2) + np.cumsum(rng.normal(0.0, 40.0, (n, 2)), axis=0)
            if kind == 1:
                xy = np.round(xy / 10.0) * 10.0
            elif kind == 2:
                xy = np.vstack([xy, xy[:1]])
            elif kind == 3:
                j = int(rng.integers(len(xy)))
                xy = np.insert(xy, j, xy[j], axis=0)
            elif kind == 4:
                xy[int(rng.integers(len(xy)))] = outer[int(rng.integers(len(outer)))] - U0
            line = LineString(xy + U0)
            if not fi._untouched(line, domain):
                continue
            parts = shapely.get_parts(shapely.intersection(line, domain))
            assert len(parts) == 1 and isinstance(parts[0], LineString), (i, line.wkt)
            assert shapely.to_wkb(parts[0]) == shapely.to_wkb(line), (i, line.wkt)
            shape = "open" if not line.is_closed else "ccw" if shapely.is_ccw(line) else "cw"
            seen[shape] += 1
        assert sum(seen.values()) >= 40, seen
        assert min(seen.values()) >= 1, seen


#: P5's L-shaped shell: it misses the region (`BOX` grown by 100 m), but its
#: envelope overlaps the region's.
ELL = [at(500, -200), at(600, -200), at(600, 600), at(-200, 600), at(-200, 500), at(500, 500)]
#: P1: `BOX` with an extra vertex at (300, 150) on its east side.
EAST_VERTEX = domain_of(Polygon([at(0, 0), at(300, 0), at(300, 150), at(300, 300), at(0, 300)]))


def pieces_in(chain: LineString, domain: Polygon) -> list[list[tuple[float, float]]]:
    """The oracle: the lines of positive length `shapely.intersection` gives."""
    parts = shapely.get_parts(shapely.intersection(chain, domain))
    return [list(p.coords) for p in parts if isinstance(p, LineString) and p.length > 0]


class TestTheClipStaysTheSame:
    """30b, section 7, pins P1-P5: green before the change and after it. The
    change skips work (the whole feature, an edge's segment test, a chain's
    intersection); these pin that no skip changes what is kept."""

    @pytest.mark.parametrize(
        ("domain", "geometry"),
        [
            pytest.param(
                BOX,
                Polygon([at(250, 200), at(200, 150), at(250, 100), at(300, 150)]),
                id="ring-touching-the-outline-at-one-vertex",
            ),
            pytest.param(
                EAST_VERTEX,
                LineString([at(300, 50), at(300, 250)]),
                id="line-along-the-outline",
            ),
            pytest.param(
                BOX,
                Polygon([at(100, 100), at(200, 100), at(200, 100), at(200, 200), at(100, 200)]),
                id="ring-with-a-repeated-vertex",
            ),
        ],
    )
    def test_p1_a_chain_the_domain_covers_but_touches_is_cut_as_geos_cuts_it(
        self, tmp_path: Path, domain: Any, geometry: Any
    ) -> None:
        """The domain covers each chain, yet GEOS's intersection does not
        give it back whole; the feature's lines are those pieces, vertex by
        vertex. A rule that skipped the intersection whenever the domain
        covers the chain (section 3.5) fails here. GEOS gives a chain back
        whole when it touches the outline only at its own start vertex, or
        runs along a side with no domain vertex in between, so the touch is
        mid-ring and the line passes `EAST_VERTEX`'s (300, 150)."""
        (chain,) = ff.boundary_lines(geometry)
        expected = pieces_in(chain, domain.polygon)
        assert domain.polygon.covers(chain) and expected != [list(chain.coords)]
        fs = one(tmp_path, geometry, domain=domain)
        assert [list(line.coords) for line in lines_of(fs)] == expected

    def test_p1_a_self_crossing_ring_is_repaired_into_its_two_lobes(self, tmp_path: Path) -> None:
        """Re-pinned at increment 111 (``a4bcab1e``; Ola, 2026-10-29: "Please
        update the pin"): the land-cover map repairs a polygon whose cut to
        the read box is invalid (``land_cover_map._fixed``), so the bowtie's
        lines are its two closed lobes, not GEOS's cut of the raw ring.
        Both give the same segments; only where a ring starts moved."""
        bowtie = Polygon([at(100, 100), at(200, 200), at(200, 100), at(100, 200)])
        fs = one(tmp_path, bowtie)
        west = [at(100, 100), at(100, 200), at(150, 150), at(100, 100)]
        east = [at(200, 200), at(200, 100), at(150, 150), at(200, 200)]
        assert [list(line.coords) for line in lines_of(fs)] == [west, east]
        (chain,) = ff.boundary_lines(bowtie)

        def segments(lines: list[list[tuple[float, float]]]) -> set[frozenset[Any]]:
            return {frozenset(s) for c in lines for s in itertools.pairwise(c)}

        assert segments([west, east]) == segments(pieces_in(chain, BOX.polygon))

    @pytest.mark.parametrize(
        "where",
        [
            pytest.param(square(50_000, 50_000, 50_100, 50_100), id="far-outside-the-region"),
            pytest.param(INNER, id="inside"),
        ],
    )
    def test_p2_a_value_the_map_lacks_is_refused_wherever_the_feature_lies(
        self, tmp_path: Path, fi: ModuleType, where: Polygon
    ) -> None:
        """Section 3.2: the refusal checks come before the whole-feature
        test, so a feature missing the region is still refused."""
        from legend_fixtures import coded_map, write

        strict = write(tmp_path / "strict.toml", coded_map("strict", unlisted="refuse"))
        path = write_geojson(tmp_path / "f.geojson", [Feat("far-9", where, {"cover": 7})])
        with pytest.raises(fi.FeatureError, match=r"far-9.*7.*strict"):
            open_one(path, BOX, str(strict))

    def test_p3_features_missing_the_region_are_counted_outside(self, tmp_path: Path) -> None:
        """The region is `BOX` grown by 100 m, so 150 m east of it is out."""
        features = [
            Feat("east", square(450, 0, 500, 50), {"Code_18": "311"}),
            Feat("in", INNER, {"Code_18": "512"}),
            Feat("far", square(50_000, 50_000, 50_100, 50_100), {"Code_18": "311"}),
        ]
        fs = open_one(write_geojson(tmp_path / "f.geojson", features), BOX)
        assert [f.code for f in fs.features] == [512]
        assert (fs.outside, fs.clipped, fs.empty) == (2, 0, 0)
        west = [Feat("west", LineString([at(-2_000, 150), at(-1_500, 150)]), {})]
        fs = open_lines(write_geojson(tmp_path / "l.geojson", west), BOX)
        assert fs.features == ()
        assert (fs.outside, fs.clipped, fs.empty) == (1, 0, 0)

    @pytest.mark.parametrize(
        ("geometry", "expected"),
        [
            pytest.param(
                Polygon(
                    square(2_000, 2_000, 2_100, 2_100).exterior.coords, [INNER.exterior.coords]
                ),
                [list(INNER.exterior.coords)],
                id="hole-inside-shell-s-envelope-misses-the-region-s",
            ),
            pytest.param(
                Polygon(ELL, [INNER.exterior.coords]),
                [list(INNER.exterior.coords)],
                id="hole-inside-shell-s-envelope-overlaps-the-region-s",
            ),
            pytest.param(
                LineString([at(-500, 150), at(800, 150)]),
                [[at(0, 150), at(300, 150)]],
                id="line-with-both-ends-outside-the-region",
            ),
        ],
    )
    def test_p5_the_whole_feature_test_keeps_what_the_pre_clip_keeps(
        self, tmp_path: Path, geometry: Any, expected: list[list[tuple[float, float]]]
    ) -> None:
        """`clip_bytes.py`'s three `probe::` cases. Each shell misses the
        region and is dropped; its hole lies in the domain and is kept whole.
        The line's ends are outside the region; its piece in the domain is
        kept. A whole-feature test on the polygon (it reads only the shell's
        envelope) loses the first; one on every feature's boundary (a line's
        is its two ends) loses the third."""
        fs = one(tmp_path, geometry)
        assert len(fs.features) == 1 and fs.outside == 0
        assert [list(line.coords) for line in lines_of(fs)] == expected


# ------------------------------------ increment 35: the clean-up step, faster

#: The CLI's clean-up defaults (`test_cli_features_cleanup.py`): repair 1 m,
#: merge on, band 50 m, outline rule 5 m.
CLEANUP = {"repair_m": 1.0, "merge_same_class": True, "tolerance_m": 50.0, "outline_snap_m": 5.0}


def cleaned(path: Path, domain: Any = BOX) -> Any:
    """``path`` as `--land-cover` with `corine` and `--edges all`, cleaned
    with `CLEANUP`."""
    return open_one(path, domain, "corine", switches=CLEANUP)


def counting_reprojector(fi: ModuleType, monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """`fi.reprojector` wrapped around the real one; the list gets the point
    count of every `apply` call."""
    import tin_engine.crs as crs_module

    calls: list[int] = []
    real = crs_module.reprojector

    def fake(src: Any, dst: Any, *, threads: int = 1) -> Any:
        move = real(src, dst, threads=threads)

        def apply(xy: Any) -> np.ndarray:
            calls.append(len(np.asarray(xy, dtype=np.float64).reshape(-1, 2)))
            return move(xy)

        return apply

    monkeypatch.setattr(fi, "reprojector", fake)
    return calls


class TestOneTransformPerSource:
    """Increment 35, section 3.1 (a) and test 2: `read_class_polygons` moves
    a source's polygons in one `shapely.transform` over all of them,
    not one call per polygon, so the threads of 3.1 (b) have an array worth
    splitting; the per-polygon `_finite` refusal keeps its text and its fid."""

    #: Three squares tiling `BOX`, in EPSG:3035, five coordinates each.
    TILES = (square(0, 0, 150, 150), square(150, 0, 300, 150), square(0, 150, 300, 300))

    def source(self, tmp_path: Path, *extra: Feat) -> Path:
        feats = [
            Feat(i + 1, moved(p, UTM33, LAEA), {"Code_18": c})
            for i, (p, c) in enumerate(zip(self.TILES, ("311", "312", "512"), strict=True))
        ]
        return write_geojson(tmp_path / "cover.geojson", [*feats, *extra], crs=LAEA)

    def test_one_apply_call_carries_every_cover_vertex(
        self, tmp_path: Path, fi: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        calls = counting_reprojector(fi, monkeypatch)
        fs = cleaned(self.source(tmp_path))
        assert len(fs.features) == 3 and fs.outside == 0
        assert calls == [sum(shapely.get_num_coordinates(p) for p in self.TILES)]  # one call, 15

    def test_a_vertex_with_no_image_is_refused_naming_its_fid(
        self, tmp_path: Path, fi: ModuleType
    ) -> None:
        """x = 1e300 in EPSG:3035 has no image (PROJ gives `inf`); the batch
        must still name the one polygon, not the first or the source."""
        xy = list(moved(self.TILES[0], UTM33, LAEA).exterior.coords)[:-1]
        xy[1] = (1e300, xy[1][1])
        path = self.source(tmp_path, Feat("lost", Polygon(xy), {"Code_18": "311"}))
        with pytest.raises(
            fi.FeatureError, match=r"^cover\.geojson: feature lost has a vertex with no image"
        ):
            cleaned(path)


class TestClippedCount:
    """Increment 35, section 3.3 and test 4: `clipped` counts the kept
    features that crossed the domain's boundary; a chain the domain covers is
    not clipped, even one along the outline, which GEOS gives back whole."""

    INSIDE = LineString([at(50, 50), at(100, 80)])
    ALSO_INSIDE = LineString([at(10, 10), at(20, 20)])
    CROSSING = LineString([at(150, 150), at(400, 150)])  # `BOX` ends at x = 300
    ALONG = LineString([at(300, 50), at(300, 250)])  # on `BOX`'s east side, no vertex between

    def edges(self, fi: ModuleType, *lines: LineString) -> Any:
        keyed = [("road", line) for line in lines]
        choice = fi.EdgeChoice(lines=frozenset({"road"}))
        region = ff.region_of(BOX)
        return fi.edges_of(None, keyed, choice, ff.cleanup(), BOX, V, region=region)

    def test_one_untouched_and_one_crossing_line_is_one_clipped_feature(
        self, fi: ModuleType
    ) -> None:
        fs = self.edges(fi, self.INSIDE, self.CROSSING)
        assert fs.clipped == 1 and len(fs.features) == 2
        assert [[list(line.coords) for line in f.lines] for f in fs.features] == [
            [list(self.INSIDE.coords)],
            [[at(150, 150), at(300, 150)]],
        ]

    def test_two_untouched_lines_are_not_clipped(self, fi: ModuleType) -> None:
        fs = self.edges(fi, self.INSIDE, self.ALSO_INSIDE)
        assert fs.clipped == 0 and len(fs.features) == 2

    def test_a_covered_line_the_rule_refused_is_not_clipped(self, fi: ModuleType) -> None:
        """A line along the outline: `_untouched` is False (its interior
        touches the boundary), yet the domain covers it and GEOS gives it
        back whole, so it is kept and not clipped. This pins that `clipped`
        is still decided by `covers` where the rule did not pass, not by the
        rule's flag alone (mutant M3)."""
        assert not fi._untouched(self.ALONG, BOX.polygon)
        fs = self.edges(fi, self.ALONG)
        assert fs.clipped == 0 and len(fs.features) == 1
        (feature,) = fs.features
        assert [list(line.coords) for line in feature.lines] == [list(self.ALONG.coords)]


def outline_segments(outline: Polygon) -> np.ndarray:
    rims = [shapely.get_coordinates(r) for r in shapely.get_rings(outline)]
    return shapely.linestrings(np.concatenate([np.stack([r[:-1], r[1:]], axis=1) for r in rims]))


def cut_on_the_outline(fi: ModuleType, polygon: Any, outline: Polygon) -> list[LineString]:
    """32, section 15.2, without a tree: every ring of `polygon` cut where an
    edge has both ends within `IN_LINE` of one outline segment, by the full
    distance matrix of ring points against outline segments."""
    segments = outline_segments(outline)
    out: list[LineString] = []
    for ring in shapely.get_rings(shapely.get_parts(polygon)):
        xy = shapely.get_coordinates(ring)
        near = shapely.distance(shapely.points(xy)[:, None], segments[None, :]) <= fi.IN_LINE
        on = (near[:-1] & near[1:]).any(axis=1)
        out += fi._chains(xy[:-1], on.tolist())
    return out


class TestOneOutlineTreePerCleanup:
    """Increment 35, section 3.4 and test 5: the clean-up builds the
    outline's segment `STRtree` once and hands it to `_outline_lines` for
    every polygon, instead of each call building its own. The CORINE extract
    over the quarter-circle domain (536 outline vertices; `TestPreClip`'s
    pairing) reaches `_outline_lines` with 8 polygons under the CLI's
    defaults."""

    @pytest.fixture
    def run(
        self, fi: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> tuple[Any, Any, list[list[str]], list[tuple[Any, tuple[LineString, ...]]]]:
        """(the edges, domain, the `feature_input.py` functions on the stack
        of each `STRtree` built over the outline's segments, nearest first,
        every `_outline_lines` call's polygon and result)."""
        import inspect

        domain = domain_of(Polygon(quarter_circle()))
        segments = len(outline_segments(domain.polygon))
        built: list[list[str]] = []
        real_tree = shapely.STRtree

        class SpiedTree(real_tree):  # type: ignore[misc, valid-type]
            def __init__(self, geoms: Any, *args: Any, **kwargs: Any) -> None:
                if len(np.asarray(geoms, dtype=object).reshape(-1)) == segments:
                    built.append(
                        [
                            f.function
                            for f in inspect.stack()
                            if Path(f.filename).name == "feature_input.py"
                        ]
                    )
                super().__init__(geoms, *args, **kwargs)

        monkeypatch.setattr(shapely, "STRtree", SpiedTree)
        calls: list[tuple[Any, tuple[LineString, ...]]] = []
        real_lines = fi._outline_lines

        def spied_lines(*args: Any) -> tuple[LineString, ...]:
            result = real_lines(*args)
            calls.append((args[0], result))
            return result

        monkeypatch.setattr(fi, "_outline_lines", spied_lines)
        fs = cleaned(EXTRACT, domain)
        return fs, domain, built, calls

    def test_the_tree_is_built_once_under_clean_not_once_per_polygon(self, run: Any) -> None:
        """`built` holds, per tree built over the outline's segments, the
        `feature_input.py` functions on the stack, nearest first. The outline
        rule's own tree (`_Outline.__init__`, under `snap_to_outline`) is not
        this one."""
        _, _, built, calls = run
        assert len(calls) == 8, len(calls)
        ours = [frames for frames in built if "snap_to_outline" not in frames]
        assert len(ours) == 1, ours
        assert "_outline_lines" not in ours[0], ours

    def test_the_lines_are_those_of_a_fresh_cut(self, fi: ModuleType, run: Any) -> None:
        _, domain, _, calls = run
        assert sum(len(lines) for _, lines in calls) == 86  # not vacuous
        for polygon, lines in calls:
            expected = cut_on_the_outline(fi, polygon, domain.polygon)
            assert [shapely.to_wkb(g) for g in lines] == [shapely.to_wkb(g) for g in expected]


def untouched_30b(line: Any, domain: Polygon) -> bool:
    """30b's `_untouched`, verbatim
    (`src_python/tin_engine/feature_input.py@c4130194:244-254` at `c7dba32e`)."""
    xy = shapely.get_coordinates(line)
    return (
        len(xy) >= 2
        and not (xy[1:] == xy[:-1]).all(axis=1).any()
        and bool(shapely.contains_properly(domain, line))
        and bool(shapely.is_simple(line))
    )


def single_call_reprojector(src: Any, dst: Any, *, threads: int = 1) -> Any:
    """`crs.reprojector` before 35, verbatim: one `transform` call. `threads`
    is accepted as the call sites pass it (38) and unused: one call, no pool."""
    del threads
    import tin_engine.crs as crs_module

    transformer = Transformer.from_crs(
        crs_module.parse_crs(src), crs_module.parse_crs(dst), always_xy=True
    )

    def apply(xy: Any) -> np.ndarray:
        points = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        x, y = transformer.transform(points[:, 0], points[:, 1])
        return np.column_stack([np.asarray(x, np.float64), np.asarray(y, np.float64)])

    return apply


class TestTheSetIsUnchanged:
    """Increment 35, test 6: the four changes skip work; none changes what
    `edges_of` gives back. The CORINE extract over the quarter circle,
    once as built and once with 30b's `_untouched` and the one-call
    reprojector put back: the same features, line for line."""

    def test_as_built_equals_the_old_rules(
        self, fi: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        domain = domain_of(Polygon(quarter_circle()))
        new = cleaned(EXTRACT, domain)
        with monkeypatch.context() as m:
            m.setattr(fi, "_untouched", untouched_30b)
            m.setattr(fi, "reprojector", single_call_reprojector)
            old = cleaned(EXTRACT, domain)
        assert len(new.features) == len(old.features) == 8
        for a, b in zip(new.features, old.features, strict=True):
            assert (a.fid, a.mask) == (b.fid, b.mask)
            assert [shapely.to_wkb(g) for g in a.lines] == [shapely.to_wkb(g) for g in b.lines]
        assert (new.clipped, new.outside, new.empty) == (old.clipped, old.outside, old.empty)
        assert new.cover_vertices == old.cover_vertices


# ------------------------------------- splits by a second attribute (ar50)


#: `docs/increments/91-ar50-vegetation-cover.md`, section 8, test 4: one
#: feature per line, (properties, the code `ar50` gives it).
SPLIT_CASES: list[tuple[dict[str, Any], int]] = [
    ({"arealtype": 50, "vegetasjonsdekke": 52}, 5052),
    ({"arealtype": 50, "vegetasjonsdekke": 52.0}, 5052),
    ({"arealtype": 50, "vegetasjonsdekke": "52"}, 5052),
    ({"arealtype": 50, "vegetasjonsdekke": 52.5}, 50),
    ({"arealtype": 50, "vegetasjonsdekke": True}, 50),
    ({"arealtype": 50, "vegetasjonsdekke": 98}, 50),
    ({"arealtype": 50}, 50),
    ({"arealtype": 50, "treslag": 31}, 50),
    ({"arealtype": 30, "treslag": 31}, 3031),
    ({"arealtype": 30, "treslag": 39}, 30),
    ({"arealtype": 30, "vegetasjonsdekke": 52}, 30),
    ({"arealtype": 60, "treslag": 31, "vegetasjonsdekke": 52}, 60),
]


def cell(i: int) -> Polygon:
    """The ``i``-th of a row of small squares inside `BOX`, 20 m apart."""
    return square(10 + 20 * i, 10, 25 + 20 * i, 25)


class TestSplitCodes:
    """Section 3, how a feature is coded, through `read_class_polygons` under `ar50`."""

    def test_each_feature_gets_its_class_or_sub_class_code(self, tmp_path: Path) -> None:
        features = [Feat(f"f{i}", cell(i), props) for i, (props, _) in enumerate(SPLIT_CASES)]
        fs = open_one(write_geojson(tmp_path / "ar50.geojson", features), BOX, "ar50")
        got = {f.fid: f.code for f in fs.features}
        assert got == {i: code for i, (_, code) in enumerate(SPLIT_CASES)}

    def test_the_masks_are_the_first_attributes(self, tmp_path: Path) -> None:
        """A sub-class's borders carry its parent's edge names: `land_cover`."""
        features = [Feat(f"f{i}", cell(i), props) for i, (props, _) in enumerate(SPLIT_CASES)]
        fs = open_one(write_geojson(tmp_path / "ar50.geojson", features), BOX, "ar50")
        assert {f.mask for f in fs.features} == {V.mask("land_cover")}
        assert len(fs.features) == len(SPLIT_CASES)

    def test_a_first_value_equal_to_a_sub_class_code_is_refused(
        self, fi: ModuleType, tmp_path: Path
    ) -> None:
        """`ar50` codes an unlisted arealtype as land, so a stray 5051 would
        otherwise mean two things in one mesh."""
        path = write_geojson(
            tmp_path / "stray.geojson",
            [Feat("ok", cell(0), {"arealtype": 50}), Feat("stray-7", cell(1), {"arealtype": 5051})],
        )
        with pytest.raises(fi.FeatureError) as caught:
            open_one(path, BOX, "ar50")
        for word in ("stray.geojson", "stray-7", "5051"):
            assert word in str(caught.value), word


class TestSplitKey:
    """Section 3, the key of a value: one function keys the loader's values
    and a feature's."""

    @pytest.mark.parametrize("value", [52, 52.0, "52", " 52 ", "52.0"], ids=repr)
    def test_a_whole_number_is_its_integer_text(self, fi: ModuleType, value: Any) -> None:
        assert fi.split_key(value) == "52"

    def test_a_text_with_letters_is_itself(self, fi: ModuleType) -> None:
        assert fi.split_key("5a") == "5a"

    @pytest.mark.parametrize(
        "value", [None, True, False, 52.5, float("nan"), "nan", [52]], ids=repr
    )
    def test_no_key(self, fi: ModuleType, value: Any) -> None:
        assert fi.split_key(value) is None


# -------------------------------- the readers' extra attributes (ar50 splits)


#: Fixture A's treslag, in fid order (`gml32_fixtures.AR50_FIGURES`).
FIXTURE_A_TRESLAG = (32, 31, 31, 31, 31, 31)
#: A box in EPSG:25832 holding every fixture-A feature.
EVERYWHERE = (-1.0e7, -1.0e7, 1.0e7, 1.0e8)


@pytest.fixture
def fixture_a_gpkg(tmp_path: Path) -> Path:
    """Fixture A converted with `gml_to_geopackage`."""
    import importlib

    from gml32_fixtures import AR50_SAMPLE
    from tin_engine.io.repository import create_sqlite

    gml_import = importlib.import_module("tin_engine.gml_import")
    target = tmp_path / "ar50.gpkg"
    with AR50_SAMPLE.open("rb") as stream, closing(create_sqlite(target)) as conn:
        gml_import.gml_to_geopackage(stream, conn)
    return target


class TestExtraFromAGeoPackage:
    """Section 8, test 5: `read_source(..., extra=...)` over fixture A."""

    def read(self, fi: ModuleType, path: Path, **kw: Any) -> Any:
        return fi.read_source(path, None, "arealtype", lambda _: EVERYWHERE, **kw)

    def test_one_tuple_per_row_in_extras_order(self, fi: ModuleType, fixture_a_gpkg: Path) -> None:
        read = self.read(fi, fixture_a_gpkg, extra=("vegetasjonsdekke", "treslag"))
        assert read.extra == tuple((98, t) for t in FIXTURE_A_TRESLAG)

    def test_the_rows_are_a_plain_reads(self, fi: ModuleType, fixture_a_gpkg: Path) -> None:
        plain = self.read(fi, fixture_a_gpkg)
        read = self.read(fi, fixture_a_gpkg, extra=("vegetasjonsdekke", "treslag"))
        assert [(f, g.wkb, v) for f, g, v in read.rows] == [(f, g.wkb, v) for f, g, v in plain.rows]
        assert all(len(row) == 3 for row in read.rows)

    def test_a_column_name_matches_in_any_case(self, fi: ModuleType, fixture_a_gpkg: Path) -> None:
        read = self.read(fi, fixture_a_gpkg, extra=("VEGETASJONSDEKKE",))
        assert read.extra == ((98,),) * 6

    def test_a_missing_column_is_none_not_a_refusal(
        self, fi: ModuleType, fixture_a_gpkg: Path
    ) -> None:
        """Section 6: a missing second attribute is no value."""
        read = self.read(fi, fixture_a_gpkg, extra=("vegetasjonsdekke", "no_such_column"))
        assert read.extra == ((98, None),) * 6

    def test_without_extra_there_is_none(self, fi: ModuleType, fixture_a_gpkg: Path) -> None:
        assert self.read(fi, fixture_a_gpkg).extra == ()

    def test_ar50_codes_the_forest_feature_by_its_tree_type(
        self, fi: ModuleType, fixture_a_gpkg: Path
    ) -> None:
        """Five marsh features (60) keep their code whatever their treslag;
        the forest feature (fid 6, treslag 31) is coniferous, 3031."""
        rows = self.read(fi, fixture_a_gpkg).rows
        x0, y0, x1, y1 = shapely.union_all([g for _, g, _ in rows]).bounds
        around = shapely.box(x0 - 500, y0 - 500, x1 + 500, y1 + 500)
        domain = domain_of(around, "EPSG:25832")
        fs = open_one(fixture_a_gpkg, domain, "ar50", dem_crs="EPSG:25832")
        assert {f.fid: f.code for f in fs.features} == {0: 60, 1: 60, 2: 60, 3: 60, 4: 60, 5: 3031}


class TestExtraFromGeoJson:
    """Section 8, test 6, the GeoJSON branch: an absent property is None."""

    def test_an_absent_property_is_none(self, fi: ModuleType, tmp_path: Path) -> None:
        path = write_geojson(
            tmp_path / "b.geojson",
            [Feat("one", cell(0), {"a": 1, "b": 52}), Feat("two", cell(1), {"a": 1})],
        )
        read = fi.read_source(path, None, "a", lambda _: EVERYWHERE, extra=("b",))
        assert read.extra == ((52,), (None,))
        assert [v for _, _, v in read.rows] == [1, 1]


class TestASplitAttributeAsWritten:
    """Section 3: a split's attribute is read under the name the legend gives
    it. The GeoJSON and GML readers match names exactly, so a capitalised
    split attribute must reach them as written; the shipped lowercase `ar50`
    splits are covered by `TestSplitCodes` and `TestExtraFromAGeoPackage`."""

    def source(self, fi: ModuleType, tmp_path: Path, data: Path, attribute: str) -> Any:
        from legend_fixtures import coded_map
        from legend_fixtures import write as write_text
        from test_legends import split, split_classes

        splits = split_classes(split(50, "Veg", (51, 5051)), split(30, "tre", (31, 3031)))
        toml = write_text(tmp_path / "caps.toml", coded_map(attribute=attribute, classes=splits))
        return data, str(toml)

    def codes(self, fi: ModuleType, source: Any) -> dict[Any, int]:
        fs = open_one(source[0], BOX, source[1])
        return {f.fid: f.code for f in fs.features}

    def test_a_geojson_property_named_with_a_capital(self, fi: ModuleType, tmp_path: Path) -> None:
        data = write_geojson(
            tmp_path / "caps.geojson",
            [
                Feat("open", cell(0), {"cover": 50, "Veg": 51}),
                Feat("forest", cell(1), {"cover": 30, "tre": 31}),
            ],
        )
        source = self.source(fi, tmp_path, data, "cover")
        assert self.codes(fi, source) == {0: 5051, 1: 3031}

    def test_a_gml_element_named_with_a_capital(self, fi: ModuleType, tmp_path: Path) -> None:
        def with_second(fid: str, i: int, code: str, element: str) -> str:
            text = member(fid, moved(cell(i), UTM33, "EPSG:4326"), code)
            return text.replace("<ogr:sl_sdeid>99</ogr:sl_sdeid>", element)

        data = tmp_path / "caps.gml"
        data.write_bytes(
            document(
                with_second("s.1", 0, "50", "<ogr:Veg>51</ogr:Veg>"),
                with_second("s.2", 1, "30", "<ogr:tre>31</ogr:tre>"),
            ).getvalue()
        )
        source = self.source(fi, tmp_path, data, "clc18_kode")
        assert self.codes(fi, source) == {0: 5051, 1: 3031}


# ------------------------- increment 111: what every polygon in the map is

#: P5's two polygons: a shell the read box misses, with a hole in the domain.
P5 = {
    "shell-s-envelope-misses-the-region-s": Polygon(
        square(2_000, 2_000, 2_100, 2_100).exterior.coords, [INNER.exterior.coords]
    ),
    "shell-s-envelope-overlaps-the-region-s": Polygon(ELL, [INNER.exterior.coords]),
}
P5_POLYGONS = [pytest.param(p, id=name) for name, p in P5.items()]
#: A bowtie 50 km off `BOX`: invalid, and its cut to the read box is empty.
FAR_BOWTIE = Polygon(
    [at(50_000, 50_000), at(50_100, 50_100), at(50_100, 50_000), at(50_000, 50_100)]
)
#: A small square of class 312 in `BOX`, clear of `INNER`.
SMALL = square(20, 20, 60, 60)


def every_area_valid(built: Any) -> bool:
    return all(bool(shapely.is_valid(g)) for _, g in built.areas)


class TestWhatEveryPolygonInTheMap:
    """`docs/increments/111-one-land-cover-map.md`, section 3, *What every
    polygon in the map is*, and section 5, test 12. The read region is `BOX`
    grown by 100 m (`source_region`), the read box that region's bounds grown
    by `land_cover_map.MARGIN` (100 m) more. After `build`, ``areas`` holds
    the polygons that reach the read box, each valid; ``outside`` the
    polygons whose cut is empty, as read: whole, never checked, never
    repaired. The edges step reads ``outside`` by rings only, after the
    areas; no polygon of it labels a triangle."""

    def build(self, pairs: list[tuple[int, BaseGeometry]], **kw: Any) -> Any:
        return ff.land_cover_map().build([pairs], [], ff.classes("corine"), ff.region_of(BOX), **kw)

    def test_premises(self) -> None:
        box = ff.land_cover_map().MARGIN + 100.0
        read_box = (at(-box, -box), at(300 + box, 300 + box))
        assert not FAR_BOWTIE.is_valid
        assert shapely.clip_by_rect(FAR_BOWTIE, *read_box[0], *read_box[1]).is_empty
        for polygon in P5.values():
            assert not polygon.is_valid
            assert shapely.clip_by_rect(polygon, *read_box[0], *read_box[1]).is_empty

    def test_a_far_bowtie_is_outside_whole_and_the_areas_are_valid(self) -> None:
        """Today `areas` holds both, the bowtie invalid (0.0 m² by GEOS)."""
        inner = Polygon(INNER.exterior.coords)
        built = self.build([(311, inner), (311, FAR_BOWTIE)])
        assert [code for code, _ in built.areas] == [311]
        assert built.areas[0][1].area == pytest.approx(10_000.0, abs=1e-6)
        assert shapely.equals(built.areas[0][1], INNER)
        assert every_area_valid(built)
        assert built.repaired == 0
        ((code, far),) = built.outside
        assert code == 311 and far is FAR_BOWTIE

    def test_a_map_built_without_fix_passes_outside_through(self) -> None:
        """The `--prefer` build (``fix=False``) is given the base map's
        ``outside`` as a keyword and passes it through; it defaults to empty."""
        base = self.build([(311, INNER), (311, FAR_BOWTIE)])
        lm, corine, region = ff.land_cover_map(), ff.classes("corine"), ff.region_of(BOX)
        lake = [(512, [square(120, 120, 180, 180)])]
        again = lm.build([base.areas], lake, corine, region, fix=False, outside=base.outside)
        assert len(again.outside) == 1 and again.outside[0][1] is FAR_BOWTIE
        assert every_area_valid(again)
        assert tuple(lm.build([base.areas], [], corine, region, fix=False).outside) == ()

    def test_the_p1_bowtie_is_valid_after_repair(self) -> None:
        bowtie = Polygon([at(100, 100), at(200, 200), at(200, 100), at(100, 200)])
        built = self.build([(311, bowtie)])
        assert len(built.areas) == 1 and tuple(built.outside) == ()
        assert every_area_valid(built)

    @pytest.mark.parametrize("polygon", P5_POLYGONS)
    def test_p5_with_cleanup_is_counted_outside(self, tmp_path: Path, polygon: Polygon) -> None:
        """Today the second raises GEOSException ("unable to assign free
        hole to a shell") in the clean-up's clip."""
        fs = one(tmp_path, polygon, switches=CLEANUP)
        assert fs.features == ()
        assert fs.outside == 1

    @pytest.mark.parametrize("polygon", P5_POLYGONS)
    def test_p5_without_cleanup_its_hole_labels_no_triangle(
        self, tmp_path: Path, polygon: Polygon
    ) -> None:
        """Its ring is kept as given (`test_p5_...`), but by its cut it has no
        area in the box: the feature carries no polygon. Today it carries the
        invalid polygon."""
        fs = one(tmp_path, polygon)
        assert len(fs.features) == 1 and fs.outside == 0
        assert [list(line.coords) for line in lines_of(fs)] == [list(INNER.exterior.coords)]
        assert fs.features[0].polygon is None

    @pytest.mark.parametrize("polygon", P5_POLYGONS)
    def test_p5_outside_polygons_come_after_the_areas(
        self, tmp_path: Path, polygon: Polygon
    ) -> None:
        """Given first in the file, p5's polygon still comes after the
        areas in the edges step's feature order. Today it comes first."""
        rows = [Feat("p5", polygon, {"Code_18": "311"}), Feat("s", SMALL, {"Code_18": "312"})]
        fs = open_one(write_geojson(tmp_path / "f.geojson", rows), BOX)
        assert [f.code for f in fs.features] == [312, 311]
        assert [f.polygon is None for f in fs.features] == [False, True]

    @pytest.mark.parametrize("polygon", P5_POLYGONS)
    def test_p5_in_a_map_cut_by_prefer_is_counted_outside(self, polygon: Polygon) -> None:
        """With the map cut by `--prefer` the clean-up runs, and each chosen
        class's outside polygon is counted outside, with no overlay: one more
        than the same map without it, and the same features."""
        fi, lm, corine = ff.feature_input(), ff.land_cover_map(), ff.classes("corine")
        region = ff.region_of(BOX)
        lake = [(512, [INNER])]

        def edges(pairs: list[tuple[int, BaseGeometry]]) -> Any:
            built = lm.build([pairs], lake, corine, region)
            assert every_area_valid(built)
            return fi.edges_of(built, [], ff.every_class(), ff.cleanup(), BOX, V, region=region)

        near = [(311, SMALL), (512, INNER)]
        with_p5, without = edges([*near, (311, polygon)]), edges(near)
        assert with_p5.outside == without.outside + 1
        assert [(f.code, f.lines) for f in with_p5.features] == [
            (f.code, f.lines) for f in without.features
        ]
