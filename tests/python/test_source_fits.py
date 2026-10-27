"""Which options a source fits: the pure rules, and their agreement with the readers.

`docs/increments/57-vector-sources.md`, section 4.4, tests T4 (`fits` over
`Facts`, one row per rule and per exclusion, and the MCP spelling of each fit
against `mcp.schema.derive`) and T5 (for each fixture file, a fit the listing
names is a file the option's reader accepts, and a fit it does not name is
a file the reader refuses).

T5 reads the files through the readers the options call, never the CLI:
`open_features` for `mesh --features` (with the fit's legend),
`read_lake_polygons` for `catchment --lakes`, `read_segments` for either
`--rivers`, `read_nve_lakes`, `read_stations` and `read_references` for
`station-catchments`' others. "Refuses" is the reader's `ValueError`; for
`catchment --lakes`, whose reader skips what is not a polygon, and for
`--features`, whose map may drop a feature, a read that keeps nothing is
counted as a refusal, since the run has nothing to use. A fit is a
header-level necessary condition: these fixtures are well formed, so the
converse is claimed for them only.
"""

from __future__ import annotations

import importlib
import json
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pytest
from pyproj import Transformer
from shapely.geometry import Point, box
from shapely.ops import transform

import batch_fixtures as bf
from catchment_fixtures import EPSG, lake_box
from dataset_fixtures import write_sources
from feature_fixtures import Feat, domain_of, rules, write_geojson
from tin_engine import source_report
from tin_engine.cli import app
from tin_engine.feature_input import (
    FeatureRequest,
    FeatureSource,
    open_features,
    read_lake_polygons,
)
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_nve_lakes, read_references, read_stations
from tin_engine.mcp.schema import derive
from vector_fixtures import LANDCOVER, fixture_v, landcover_gpkg, vector_entry

Fit = tuple[str, str, str | None]

DEMS = {("mesh", "--dem", None), ("catchment", "--dem", None)}
ALL_DEMS = DEMS | {("station-catchments", "--dem", None)}
LAKES = ("catchment", "--lakes", None)
SC_LAKES = ("station-catchments", "--lakes", None)
RIVERS = {("catchment", "--rivers", None), ("station-catchments", "--rivers", None)}
STATIONS = ("station-catchments", "--stations", None)
REFERENCE = ("station-catchments", "--reference", None)
#: A polygon source is a network of units; a point source, their outlet points.
NETWORK = ("network-catchment", "--network", None)
OUTLETS = ("network-catchment", "--outlets", None)
#: Any polygons fit `mesh --reservoirs`; projected ones `mesh --catchments` too,
#: whose border band is in metres.
RESERVOIRS = ("mesh", "--reservoirs", None)
CATCHMENTS = ("mesh", "--catchments", None)
#: What any polygons fit.
POLYGONS_FIT = (LAKES, NETWORK, RESERVOIRS)


def features(*maps: str) -> set[Fit]:
    return {("mesh", "--features", m) for m in maps}


def fits_for(fields: dict[str, Any]) -> Any:
    return source_report.fits(source_report.Facts(**fields), maps())


def maps() -> dict[str, Any]:
    """What `listing` passes: every shipped legend plus `property`."""
    shipped = importlib.import_module("tin_engine.meshing_rules").shipped_names()
    return {name: rules(name) for name in (*shipped, "property")}


def facts(
    values: str = "features",
    *,
    suffix: str = ".geojson",
    geometry: tuple[str, ...] = (),
    attributes: tuple[str, ...] = (),
    crs: bool = True,
    projected: bool | None = None,
    rules: str = "",
) -> dict[str, Any]:
    """`Facts`' fields; built in the test, so the module collects before `Facts` exists."""
    return dict(
        values=values,
        projected=projected,
        suffix=suffix,
        geometry=frozenset(geometry),
        attributes=attributes,
        crs_stated=crs,
        rules=rules,
    )


def gpkg(attributes: tuple[str, ...], geometry: str = "Geometry", **kw: Any) -> dict[str, Any]:
    return facts(suffix=".gpkg", geometry=(geometry,), attributes=attributes, **kw)


POLY = ("Polygon",)
POLYS = ("Polygon", "MultiPolygon")

#: (Facts, the fits expected), one row per rule of section 4.4 and per exclusion.
ROWS: dict[str, tuple[dict[str, Any], set[Fit]]] = {
    "dem projected": (facts("heights", suffix="", crs=False, projected=True), ALL_DEMS),
    "dem crs unknown": (facts("heights", suffix="", crs=False, projected=None), ALL_DEMS),
    "dem in degrees, not station-catchments": (
        facts("heights", suffix="", crs=False, projected=False), DEMS),
    "classes": (facts("classes", suffix="", crs=False), {("mesh", "--land-cover", None)}),
    "land cover gpkg, every map of its attribute": (
        gpkg(("Code_18",)), features("corine", "corine-water") | {*POLYGONS_FIT}),
    "gpkg attribute case-insensitive": (
        gpkg(("code_18", "Remark")), features("corine", "corine-water") | {*POLYGONS_FIT}),
    "geojson attribute exact": (
        facts(geometry=POLY, attributes=("code_18",)), {*POLYGONS_FIT, SC_LAKES}),
    "geojson attribute as the map spells it": (
        facts(geometry=POLY, attributes=("Code_18",)),
        features("corine", "corine-water") | {*POLYGONS_FIT, SC_LAKES}),
    "the entry's rules alone": (
        gpkg(("Code_18",), rules="corine"), features("corine") | {*POLYGONS_FIT}),
    "the entry's rules, its attribute absent": (
        gpkg(("Code_18",), rules="clc18_kode"), {*POLYGONS_FIT}),
    "another map's attribute": (gpkg(("clc18_kode",)), features("clc18_kode") | {*POLYGONS_FIT}),
    "mapbiomas classes": (facts(geometry=POLY, attributes=("class",)),
                          features("mapbiomas-c11") | {*POLYGONS_FIT, SC_LAKES}),
    "an uncoded property map": (facts(geometry=POLY, attributes=("property",), crs=False),
                                features("property") | {*POLYGONS_FIT}),
    "rivers": (facts(geometry=("LineString",), attributes=("objectid", "elvid")), RIVERS),
    "rivers in a gpkg: no rivers": (gpkg(("objectid",), geometry="LineString"), set()),
    "rivers with no crs member: no rivers": (
        facts(geometry=("LineString",), attributes=("objectid",), crs=False), set()),
    "MultiLineString: no rivers": (
        facts(geometry=("MultiLineString",), attributes=("objectid",)), set()),
    "rivers with no objectid: no rivers": (
        facts(geometry=("LineString",), attributes=("elvid",)), set()),
    "lines and polygons: lakes, no rivers": (
        facts(geometry=("LineString", "Polygon"), attributes=("objectid",)), {*POLYGONS_FIT}),
    "stations": (facts(geometry=("Point",), attributes=("station", "name")), {STATIONS, OUTLETS}),
    "stations in .json": (
        facts(suffix=".json", geometry=("Point",), attributes=("station",)), {STATIONS, OUTLETS}),
    "stations with no crs member: outlets only": (
        facts(geometry=("Point",), attributes=("station",), crs=False), {OUTLETS}),
    "stations with no station attribute: outlets only": (
        facts(geometry=("Point",), attributes=("name",)), {OUTLETS}),
    "points and polygons: lakes, network and outlets": (
        facts(geometry=("Point", "Polygon"), attributes=("station",)),
        {*POLYGONS_FIT, OUTLETS}),
    "references": (facts(geometry=POLYS, attributes=("station",)),
                   {*POLYGONS_FIT, SC_LAKES, REFERENCE}),
    "lakes": (facts(geometry=POLY, attributes=("objectid", "navn")), {*POLYGONS_FIT, SC_LAKES}),
    "lakes with no crs member: catchment lakes and network": (
        facts(geometry=POLYS, attributes=("name",), crs=False), {*POLYGONS_FIT}),
    "points in a gpkg: outlets": (gpkg(("unit",), geometry="Point"), {OUTLETS}),
    "a lakes gpkg: catchment lakes and network": (
        gpkg(("navn",), geometry="MultiPolygon"), {*POLYGONS_FIT}),
    "projected polygons: --catchments too": (
        facts(geometry=POLY, attributes=("name",), projected=True),
        {*POLYGONS_FIT, CATCHMENTS, SC_LAKES}),
    "polygons in degrees: --reservoirs, no --catchments": (
        facts(geometry=POLYS, attributes=("name",), projected=False),
        {*POLYGONS_FIT, SC_LAKES}),
    "a projected polygon gpkg: --catchments too": (
        gpkg(("navn",), geometry="MultiPolygon", projected=True),
        {*POLYGONS_FIT, CATCHMENTS}),
    "projected points: neither": (
        facts(geometry=("Point",), attributes=("name",), projected=True), {OUTLETS}),
}  # fmt: skip


class TestT4Fits:
    @pytest.mark.parametrize(("given", "expected"), ROWS.values(), ids=ROWS.keys())
    def test_the_rules(self, given: dict[str, Any], expected: set[Fit]) -> None:
        got = fits_for(given)
        listed = [(f.command, f.option, f.rules) for f in got]
        assert len(listed) == len(set(listed)), listed  # each fit once
        assert set(listed) == expected

    def test_fits_read_no_file(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Pure: opening anything fails the test."""

        def no_open(*args: Any, **kwargs: Any) -> Any:
            raise AssertionError("fits opened a file")

        given = maps()
        built = [source_report.Facts(**row) for row, _expected in ROWS.values()]
        monkeypatch.setattr("builtins.open", no_open)
        monkeypatch.setattr(Path, "open", no_open)
        for row in built:
            source_report.fits(row, given)

    def test_each_fits_mcp_spelling_is_the_servers(self) -> None:
        tools = {d.tool_name: d.command for d in derive(app)}
        seen = set()
        for row, _expected in ROWS.values():
            for fit in fits_for(row):
                command = tools[fit.tool]
                assert command.name == fit.command, fit
                (param,) = [p for p in command.params if p.name == fit.argument]
                assert fit.option in param.opts, (fit, param.opts)
                seen.add((fit.command, fit.option, fit.tool, fit.argument))
        assert ("station-catchments", "--dem", "station_catchments", "dem") in seen
        assert ("mesh", "--features", "mesh", "features") in seen

    def test_only_a_features_fit_names_rules(self) -> None:
        for row, _expected in ROWS.values():
            for fit in fits_for(row):
                assert (fit.rules is not None) == (fit.option == "--features"), fit


# ------------------------------------------------------------ T5

#: A point inside the lake and the land cover, in EPSG:25833: the lake box's centre.
INSIDE_POINT = lake_box().centroid.coords[0]


def accepts_features(path: Path, layer: str | None, name: str) -> bool:
    source = FeatureSource(path=path, rules=rules(name), layer=layer)
    domain = domain_of(box(495_000.0, 6_555_000.0, 525_000.0, 6_605_000.0))
    try:
        found = open_features(FeatureRequest(sources=(source,)), domain, EPSG)
    except ValueError:
        return False
    return len(found.features) > 0


def accepts(read: Callable[[Path], object]) -> Callable[[Path, str | None], bool]:
    def judged(path: Path, layer: str | None) -> bool:
        if layer is not None:
            return False  # these readers take no layer: GeoJSON only
        try:
            read(path)
        except ValueError:
            return False
        return True

    return judged


def accepts_lakes(path: Path, layer: str | None) -> bool:
    try:
        found, _crs = read_lake_polygons(path, layer, INSIDE_POINT, EPSG)
    except ValueError:
        return False
    return len(found) > 0


READERS: Mapping[tuple[str, str], Callable[[Path, str | None], bool]] = {
    ("catchment", "--lakes"): accepts_lakes,
    ("catchment", "--rivers"): accepts(read_segments),
    ("station-catchments", "--rivers"): accepts(read_segments),
    ("station-catchments", "--lakes"): accepts(read_nve_lakes),
    ("station-catchments", "--stations"): accepts(read_stations),
    ("station-catchments", "--reference"): accepts(read_references),
}


def without_crs(path: Path) -> Path:
    doc = json.loads(path.read_text())
    del doc["crs"]
    target = path.with_name("no-crs-" + path.name)
    target.write_text(json.dumps(doc))
    return target


@pytest.fixture
def listed(tmp_path: Path) -> dict[str, tuple[Path, str | None, set[Fit]]]:
    """Every fixture file registered and listed: id -> (path, layer, fits)."""
    lakes, cover = fixture_v(tmp_path / "v")
    files: dict[str, tuple[Path, str | None]] = {
        "lakes": (lakes, None),
        "cover": (cover, LANDCOVER),
        "stations": (bf.write_stations(tmp_path / "stations.geojson"), None),
        "rivers": (bf.write_rivers(tmp_path / "rivers.geojson"), None),
        "reference": (bf.write_references(tmp_path / "reference.geojson"), None),
        # Not "nve-lakes": that id is the package's map-service entry.
        "hrd-lakes": (bf.write_lakes(tmp_path / "nve-lakes.geojson"), None),
        "lakes-no-crs": (without_crs(lakes), None),
    }
    user = write_sources(
        tmp_path / "conf" / "sources.toml",
        {key: vector_entry(path, f"fixture {key}", layer=layer)
         for key, (path, layer) in files.items()},
    )  # fmt: skip
    document = source_report.listing({"RASPUTIN_DATA": str(user.parent)}, None, None)
    found = {s.id: s for s in document.sources}
    out = {}
    for key, (path, layer) in files.items():
        entry = found[key]
        assert entry.problem is None, (key, entry.problem)
        out[key] = (path, layer, {(f.command, f.option, f.rules) for f in entry.fits})
    return out


class TestT5FitsAgreeWithTheReaders:
    def test_every_vector_option_on_every_file(
        self, listed: dict[str, tuple[Path, str | None, set[Fit]]]
    ) -> None:
        disagreements = []
        for key, (path, layer, fits) in listed.items():
            for (command, option), reader in READERS.items():
                named = (command, option, None) in fits
                if named != reader(path, layer):
                    disagreements.append((key, option, command, "fit" if named else "no fit"))
            for name in maps():
                named = ("mesh", "--features", name) in fits
                if named != accepts_features(path, layer, name):
                    disagreements.append((key, "--features", name, "fit" if named else "no fit"))
        assert disagreements == []

    def test_the_fixtures_exercise_both_sides(
        self, listed: dict[str, tuple[Path, str | None, set[Fit]]]
    ) -> None:
        """Each vector option is named for some file and not for another, so
        neither side of the agreement is vacuous."""
        options = {(c, o) for _p, _l, fits in listed.values() for c, o, _m in fits}
        assert set(READERS) | {("mesh", "--features")} <= options, options
        for command, option in READERS:
            assert any((command, option, None) not in f for _p, _l, f in listed.values())

    def test_a_geopackage_and_its_layer(self, tmp_path: Path) -> None:
        """The land cover fits `catchment --lakes` (polygons; the description
        says what they are), and no GeoJSON-only reader."""
        path = landcover_gpkg(tmp_path / "landcover.gpkg")
        assert accepts_lakes(path, LANDCOVER)
        assert accepts_features(path, LANDCOVER, "corine")
        assert not accepts_features(path, LANDCOVER, "clc18_kode")


# ------------------------------------------------------------ 106e, test 8


def test_106e_the_listing_names_mesh_for_polygons(tmp_path: Path) -> None:
    """`docs/increments/106-wednesday-cells.md`, section 8, test 8: a projected
    polygon file fits `mesh --reservoirs` and `mesh --catchments`, one in
    EPSG:4326 `--reservoirs` only, a points file neither; the CRS is the one
    the listing reads from the file."""
    unit = box(500_000.0, 6_600_000.0, 501_000.0, 6_601_000.0)
    to_wgs84 = Transformer.from_crs("EPSG:25833", "EPSG:4326", always_xy=True).transform
    files = {
        "projected": write_geojson(
            tmp_path / "units.geojson", [Feat(1, unit, {"name": "u"})], "EPSG:25833"
        ),
        "degrees": write_geojson(
            tmp_path / "units-deg.geojson",
            [Feat(1, transform(to_wgs84, unit), {"name": "u"})],
            "EPSG:4326",
        ),
        "points": write_geojson(
            tmp_path / "points.geojson",
            [Feat(k, Point(c), {"name": f"p{k}"}) for k, c in enumerate(unit.exterior.coords[:2])],
            "EPSG:25833",
        ),
    }
    user = write_sources(
        tmp_path / "conf" / "sources.toml",
        {key: vector_entry(path, f"fixture {key}") for key, path in files.items()},
    )
    document = source_report.listing({"RASPUTIN_DATA": str(user.parent)}, None, None)
    found = {s.id: s for s in document.sources}
    mesh_fits = {
        key: {(f.command, f.option) for f in found[key].fits if f.command == "mesh"}
        for key in files
    }
    for key in files:
        assert found[key].problem is None, (key, found[key].problem)
    assert mesh_fits == {
        "projected": {("mesh", "--reservoirs"), ("mesh", "--catchments")},
        "degrees": {("mesh", "--reservoirs")},
        "points": set(),
    }
