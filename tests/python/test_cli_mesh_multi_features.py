"""`rasputin mesh` with several `--features` sources in one mesh: increment 16e.

`docs/increments/16e-multi-features.md`: R1 (the flags become repeatable and
pair by position), R2/R3 (`_feature_sources` and its degeneracy policy), R4
(the noder merges cross-source coincident edges, unioning the masks), R5/D1
(labelling spans coded sources; two different non-empty code systems are
refused, a coded plus a code-less source is allowed), R6 (one file record per
source, one distinct notice), and "Tests, and the regression floor" (R7 cases
1-5, 7). Invariants I1-I3, I5.

Reuses `test_cli_mesh_features.py`'s fixtures and helpers (`bumpy`,
`plain_square`, `mesh`, `meshed`, `edges`, `edges_at`, `text_field`,
`FEATURES_FIELD`, the gallery geometries) so the multi-source cases read the
same way the single-source ones do.

This suite is RED before any production code: `--features` and its paired
options are single-valued (`Path | None` / `str | None`, `cli.py:637-660`), so
a second `--features` overrides the first rather than adding a source, and
`_feature_sources` does not exist. The suite fails because the new multi-source
behaviour is absent, not on an unknown option: Typer accepts a repeated
`--features` (it just keeps the last), so these runs execute and give a
single-source result the assertions reject. `_feature_sources`'s import in
`test_feature_input.py`'s companion case fails on the missing symbol.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon

import feature_fixtures as ff
from feature_fixtures import Feat, write_geojson
from geotiff_fixtures import micro_tiff
from gpkg_fixtures import Layer, Row, needs_rtree, write_gpkg
from landcover_fixtures import vtk_labels
from test_cli_mesh_dem import USAGE, invoke, write_tiff
from test_cli_mesh_domain import COLS, ROWS, SQUARE, geojson
from test_cli_mesh_features import (
    FEATURES_FIELD,
    SNAP,
    edges,
    edges_at,
    mesh,
    meshed,
    rel,
    text_field,
)
from tin_engine.features import DEFAULT_VOCABULARY
from vtkread import VtkFile, read_vtk

V = DEFAULT_VOCABULARY

# Two disjoint coded polygons inside `SQUARE`, well off every DEM node: a
# forest (311) in the west, a lake (512, water) in the east. Authored in the
# DEM's own CRS (write_geojson's default EPSG:25833), so no reprojection.
FOREST = Polygon([rel(30.3, -60.1), rel(60.7, -60.3), rel(60.9, -20.2), rel(30.1, -20.4)])
LAKE = Polygon([rel(120.3, -60.2), rel(170.1, -60.4), rel(170.3, -20.1), rel(120.1, -20.3)])
# An uncoded `property` boundary: a wall crossing the middle of the domain, its
# edge a constraint that blocks the flood fill but names no class.
WALL = LineString([rel(15.2, -35.1), rel(185.3, -35.1)])

# R4 / I2: an edge shared bit-for-bit between two sources, both authored in the
# DEM's CRS, so the noder meets identical coordinates and merges them. The
# forest's east side and the lake's west side are one segment.
SHARED_X = rel(90.5, 0.0)[0]
LEFT = Polygon([rel(30.4, -60.3), (SHARED_X, rel(0, -60.3)[1]), (SHARED_X, rel(0, -20.1)[1]),
                rel(30.2, -20.2)])  # fmt: skip
RIGHT = Polygon([(SHARED_X, rel(0, -60.3)[1]), rel(160.6, -60.4), rel(160.8, -20.3),
                 (SHARED_X, rel(0, -20.1)[1])])  # fmt: skip
SEAM = LineString([(SHARED_X, rel(0, -55.1)[1]), (SHARED_X, rel(0, -25.3)[1])])


def coded(fid: object, geometry: Polygon | LineString, code: object) -> Feat:
    return Feat(fid, geometry, {"Code_18": code})


def prop(fid: object, geometry: Polygon | LineString, klass: str) -> Feat:
    return Feat(fid, geometry, {"property": klass})


@pytest.fixture
def bumpy(tmp_path: Path) -> Path:
    array = np.random.default_rng(16).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))


@pytest.fixture
def plain_square(tmp_path: Path) -> Path:
    return geojson(tmp_path / "square.geojson", SQUARE)


@pytest.fixture
def corine_src(tmp_path: Path) -> Path:
    """A CORINE-coded source: forest 311 and lake 512."""
    return write_geojson(
        tmp_path / "corine.geojson", [coded("forest", FOREST, "311"), coded("lake", LAKE, "512")]
    )


@pytest.fixture
def parcel_src(tmp_path: Path) -> Path:
    """A code-less `property` source: one wall boundary."""
    return write_geojson(tmp_path / "parcel.geojson", [prop("wall", WALL, "wall")])


def two_sources(a: Path, b: Path, *, a_map: str = "property", b_map: str = "property") -> list[str]:
    """`--features a --features-map <a_map> --features b --features-map <b_map>`,
    each map written next to its source (R1: no broadcast)."""
    return [
        "--features", str(a), "--features-map", a_map,
        "--features", str(b), "--features-map", b_map,
    ]  # fmt: skip


# --------------------------------------------------------- R7 case 1: two sources


class TestTwoSources:
    """R7 case 1 / I3: a CORINE-coded source and an uncoded boundary source in
    one mesh. The mesh carries both sets of constraint edges, the CORINE
    source's `land_cover_code`, and the boundary source's edges block the fill."""

    @pytest.fixture
    def vtk(self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path,
            parcel_src: Path) -> VtkFile:  # fmt: skip
        return meshed(
            tmp_path, bumpy, plain_square,
            *two_sources(corine_src, parcel_src, a_map="corine", b_map="property"),
        )[0]  # fmt: skip

    def test_both_sources_constraints_are_present(self, vtk: VtkFile) -> None:
        _, masks = edges(vtk)
        found = set(masks.tolist())
        # CORINE: land_cover on the forest, land_cover+water on the lake.
        assert {V.mask("land_cover"), V.mask("land_cover", "water")} <= found
        # The parcel wall carries its own `wall` bit, from the second source.
        assert V.mask("wall") in found

    def test_the_wall_edge_is_a_constraint_alongside_the_coded_edges(self, vtk: VtkFile) -> None:
        """The parcel's wall (second source) and the CORINE edges (first
        source) are both present: neither source overrode the other."""
        xy, masks = edges(vtk)
        wall = xy[(masks & V.mask("wall")) != 0]
        assert len(wall), "the parcel source contributed no wall constraint"
        land = xy[(masks & V.mask("land_cover")) != 0]
        assert len(land), "the CORINE source contributed no land_cover constraint"

    def test_the_coded_source_labels_the_cells(self, vtk: VtkFile) -> None:
        """I3: the CORINE polygons label the cells; the uncoded parcel adds no
        code but its boundary still partitions the fill."""
        _, _, codes = vtk_labels(vtk, [(FOREST, 311), (LAKE, 512)])
        found = set(codes.tolist())
        assert 311 in found and 512 in found

    def test_the_file_records_both_sources(self, vtk: VtkFile) -> None:
        """R6 / I5: the `features` field has one entry per source, in order."""
        text = text_field(vtk, "features")
        entries = [e.strip() for e in text.split(";")]
        assert len(entries) == 2, text
        assert FEATURES_FIELD.fullmatch(entries[0])["map"] == "corine", entries  # type: ignore[index]
        assert FEATURES_FIELD.fullmatch(entries[1])["map"] == "property", entries  # type: ignore[index]


# --------------------------------------------------- R7 case 2: positional pairing


class TestPositionalPairing:
    """R1 / I1: the Nth `--features` takes the Nth of each paired option; a
    paired list shorter than `--features` defaults the unpaired sources."""

    def test_one_map_two_sources_defaults_the_second(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        """`--features a --features b --features-map corine`: `a` gets `corine`,
        `b` takes the default `property` (the shorter list defaults)."""
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square,
            "--features", str(corine_src), "--features", str(parcel_src),
            "--features-map", "corine",
        )  # fmt: skip
        entries = [e.strip() for e in text_field(vtk, "features").split(";")]
        assert len(entries) == 2, entries
        assert FEATURES_FIELD.fullmatch(entries[0])["map"] == "corine", entries  # type: ignore[index]
        assert FEATURES_FIELD.fullmatch(entries[1])["map"] == "property", entries  # type: ignore[index]
        # `corine` applied to the parcel would refuse its `wall` value; a clean
        # run proves `b` took `property`.
        _, masks = edges(vtk)
        assert V.mask("wall") in set(masks.tolist())

    def test_a_crs_given_for_the_second_source_only(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """`--features-crs` for the second source only pairs with the second:
        the first reads its own file CRS, the second is told EPSG:4326."""
        first = write_geojson(tmp_path / "a.geojson", [prop("wall", WALL, "wall")])
        # The second source's geometry authored in EPSG:4326, its file carrying
        # no crs member, so it must be told its CRS by the paired option.
        lonlat = ff.moved(FOREST, "EPSG:25833", "EPSG:4326")
        second = write_geojson(tmp_path / "b.geojson", [prop("f", lonlat, "land_cover")], crs=None)
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square,
            "--features", str(first),
            "--features", str(second), "--features-crs", "EPSG:4326",
        )  # fmt: skip
        crss = [c.strip() for c in text_field(vtk, "features_crs").split(";")]
        assert crss == ["EPSG:25833", "EPSG:4326"], crss

    def test_one_features_with_no_paired_options_is_unchanged(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path
    ) -> None:
        """The regression guard: a single `--features` with a single
        `--features-map` still means one source with that map."""
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square, "--features", str(corine_src), "--features-map", "corine"
        )
        text = text_field(vtk, "features")
        assert ";" not in text, text
        assert FEATURES_FIELD.fullmatch(text)["map"] == "corine", text  # type: ignore[index]


# ------------------------------------------------ R7 case 3: the same file twice


class TestSameFileTwice:
    """R3 / D4: the same path given twice is allowed (two layers or two maps of
    one GeoPackage); the noder merges the duplicated edges."""

    @needs_rtree
    def test_two_layers_of_one_geopackage(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        rows_a = [Row(1, FOREST, {"Code_18": "311"})]
        rows_b = [Row(1, LAKE, {"Code_18": "512"})]
        gpkg = write_gpkg(
            tmp_path / "two.gpkg",
            [
                Layer("forest", 25833, rows_a, columns=("Code_18",)),
                Layer("lake", 25833, rows_b, columns=("Code_18",)),
            ],
        )
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square,
            "--features", str(gpkg), "--features-layer", "forest", "--features-map", "corine",
            "--features", str(gpkg), "--features-layer", "lake", "--features-map", "corine",
        )  # fmt: skip
        _, masks = edges(vtk)
        found = set(masks.tolist())
        assert {V.mask("land_cover"), V.mask("land_cover", "water")} <= found

    def test_same_geojson_twice_two_maps(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path
    ) -> None:
        """The same coded file read once as `corine` and once as
        `corine-water`; both are the same code system (allowed under D1)."""
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square,
            "--features", str(corine_src), "--features-map", "corine",
            "--features", str(corine_src), "--features-map", "corine-water",
        )  # fmt: skip
        entries = [e.strip() for e in text_field(vtk, "features").split(";")]
        assert len(entries) == 2, entries
        # The lake edge is present from both reads, merged, not doubled.
        _, masks = edges(vtk)
        assert V.mask("land_cover", "water") in set(masks.tolist())


# --------------------------- R7 case 4 / I2: cross-source edge merge (R4)


class TestCrossSourceMerge:
    """R4 / I2: two sources with an edge coincident in the DEM's CRS merge into
    one constraint edge whose mask is the union of the two maps' bits. The
    cross-source analogue of 16b's within-source shared-edge merge."""

    def seam_mesh(self, tmp_path: Path, bumpy: Path, plain_square: Path) -> VtkFile:
        """LEFT is a lake (water) in one CORINE source; RIGHT is forest
        (land_cover) in a second. Their shared vertical seam is authored
        identically in both, in the DEM's CRS."""
        left = write_geojson(tmp_path / "left.geojson", [coded("l", LEFT, "512")])
        right = write_geojson(tmp_path / "right.geojson", [coded("r", RIGHT, "311")])
        return meshed(
            tmp_path, bumpy, plain_square,
            "--features", str(left), "--features-map", "corine",
            "--features", str(right), "--features-map", "corine",
        )[0]  # fmt: skip

    def test_the_shared_seam_carries_both_masks(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        vtk = self.seam_mesh(tmp_path, bumpy, plain_square)
        xy, masks = edges(vtk)
        # Every edge lying on the seam carries both sources' bits, merged.
        mids = (xy[:, 0] + xy[:, 1]) / 2
        on_seam = np.asarray(shapely.distance(SEAM, shapely.points(mids))) <= SNAP
        assert on_seam.any(), "no edge found on the shared seam"
        want = V.mask("land_cover", "water")
        assert set(masks[on_seam].tolist()) == {want}, sorted(set(masks[on_seam].tolist()))

    def test_the_seam_is_one_run_of_edges_not_doubled(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """The merge collapses the two coincident chains to one edge run: no
        edge on the seam carries only one source's bit."""
        vtk = self.seam_mesh(tmp_path, bumpy, plain_square)
        masks_here = edges_at(vtk, (SHARED_X, rel(0, -40.2)[1]))
        assert masks_here, "no constraint edge at the middle of the seam"
        assert set(masks_here) == {V.mask("land_cover", "water")}, masks_here


# ------------------------------------------------------ R7 case 5: degeneracy (R3)


class TestDegeneracy:
    """R3: every degenerate case a clear `typer.BadParameter` naming the flag,
    and the 1-based index/path the design says so."""

    def test_a_map_longer_than_features_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path
    ) -> None:
        """Two `--features-map`, one `--features`: refused naming the flag and
        the counts (R3, the one length rule)."""
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            "--features", str(corine_src), "--features-map", "corine", "--features-map", "property",
        )  # fmt: skip
        assert code == USAGE, output
        assert "--features-map" in output and "No such option" not in output, output
        assert re.search(r"\b2\b.*\b1\b", output), output
        assert not target.exists()

    def test_a_crs_longer_than_features_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path
    ) -> None:
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            "--features", str(corine_src),
            "--features-crs", "EPSG:25833", "--features-crs", "EPSG:4326",
        )  # fmt: skip
        assert code == USAGE and "--features-crs" in output, output
        assert "No such option" not in output, output
        assert not target.exists()

    def test_a_map_with_no_features_is_refused(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """A paired option with no `--features` at all: the same length rule
        with M = 0 (generalises the old "applies only with --features")."""
        code, output, target = mesh(tmp_path, bumpy, plain_square, "--features-map", "corine")
        assert code == USAGE and "--features-map" in output, output
        assert "No such option" not in output, output
        assert not target.exists()

    def test_features_with_no_domain_is_refused(
        self, tmp_path: Path, bumpy: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        """R3: `--features` given with no `--domain`, unchanged, with two."""
        code, output = invoke(
            "--dem", str(bumpy),
            "--features", str(corine_src), "--features", str(parcel_src),
            "--out", str(bumpy.parent / "x.vtk"),
        )  # fmt: skip
        assert code == USAGE and "--features" in output, output
        assert "No such option" not in output, output

    def test_an_unknown_map_at_index_two_is_refused_naming_the_index(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            "--features", str(corine_src), "--features-map", "corine",
            "--features", str(parcel_src), "--features-map", "nosuch",
        )  # fmt: skip
        assert code == USAGE, output
        assert "nosuch" in output and "2" in output, output
        assert "No such option" not in output, output
        assert not target.exists()

    def test_a_layer_on_a_geojson_at_index_two_is_refused_naming_index_and_path(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        """A `--features-layer` on a non-`.gpkg` source at index 2: refused
        naming the index and the path (R3)."""
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            "--features", str(corine_src),
            "--features", str(parcel_src), "--features-layer", "x",
        )  # fmt: skip
        assert code == USAGE, output
        assert "--features-layer" in output and "No such option" not in output, output
        assert "2" in output and parcel_src.name in output, output
        assert not target.exists()

    def test_a_shorter_paired_list_is_allowed(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        """The common case: a paired option shorter than `--features` defaults
        the missing entries, no error (R3)."""
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            "--features", str(corine_src), "--features-map", "corine",
            "--features", str(parcel_src),
        )  # fmt: skip
        assert code == 0, output
        assert target.exists()


# --------------------------------------------------- D1: code systems across sources


class TestCodeSystems:
    """D1: two sources with different non-empty code systems are refused (one
    mesh, one `land_cover_code` system); a coded source plus a code-less
    (`property`) source is allowed — the motivating CORINE + parcels case."""

    def test_a_coded_source_and_a_codeless_source_is_allowed(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        """The motivating case: CORINE (coded) plus a parcel boundary
        (code-less `property`). Allowed, and only the coded source labels."""
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            *two_sources(corine_src, parcel_src, a_map="corine", b_map="property"),
        )  # fmt: skip
        assert code == 0, output
        vtk = read_vtk(target.read_bytes())
        assert "land_cover_codes" in vtk.field_data
        assert V.mask("wall") in set(edges(vtk)[1].tolist())

    def test_two_same_system_coded_sources_is_allowed(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """Both CORINE maps carry the same `codes` system; mixing them is fine
        (R5: the refusal never fires on today's maps)."""
        a = write_geojson(tmp_path / "a.geojson", [coded("f", FOREST, "311")])
        b = write_geojson(tmp_path / "b.geojson", [coded("l", LAKE, "512")])
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            "--features", str(a), "--features-map", "corine",
            "--features", str(b), "--features-map", "corine-water",
        )  # fmt: skip
        assert code == 0, output
        vtk = read_vtk(target.read_bytes())
        assert "land_cover_codes" in vtk.field_data
        # Both sources' polygons label the mesh: neither overrode the other.
        _, _, codes = vtk_labels(vtk, [(FOREST, 311), (LAKE, 512)])
        assert {311, 512} <= set(codes.tolist()), sorted(set(codes.tolist()))


# ------------------------------------------------------------ R6 / I5: the record


class TestRecordPerSource:
    """R6 / I5: one `features` and one `features_crs` entry per source in source
    order; `features_notice` lists each distinct notice once; per-source `k
    features` counts are correct."""

    def test_two_corine_sources_give_one_notice(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """Both sources use CORINE maps carrying the same Copernicus notice;
        the record names it once, not twice."""
        a = write_geojson(tmp_path / "a.geojson", [coded("f", FOREST, "311")])
        b = write_geojson(tmp_path / "b.geojson", [coded("l", LAKE, "512")])
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square,
            "--features", str(a), "--features-map", "corine",
            "--features", str(b), "--features-map", "corine-water",
        )  # fmt: skip
        # Both sources are present (two file entries), and the shared notice is
        # named once, not per source.
        assert len(text_field(vtk, "features").split(";")) == 2, text_field(vtk, "features")
        notice = text_field(vtk, "features_notice")
        assert notice.count("Copernicus Land Monitoring Service") == 1, notice

    def test_per_source_feature_counts_are_correct(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """The first source has two features, the second one; each entry's
        `<k> features` count is that source's own, not the run total."""
        a = write_geojson(
            tmp_path / "a.geojson", [coded("f", FOREST, "311"), coded("l", LAKE, "512")]
        )
        b = write_geojson(tmp_path / "b.geojson", [prop("w", WALL, "wall")])
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square,
            "--features", str(a), "--features-map", "corine",
            "--features", str(b), "--features-map", "property",
        )  # fmt: skip
        entries = [e.strip() for e in text_field(vtk, "features").split(";")]
        assert len(entries) == 2, entries
        counts = [int(FEATURES_FIELD.fullmatch(e)["n"]) for e in entries]  # type: ignore[index]
        assert counts == [2, 1], (counts, entries)

    def test_features_crs_has_one_entry_per_source(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        vtk, _ = meshed(
            tmp_path, bumpy, plain_square,
            *two_sources(corine_src, parcel_src, a_map="corine", b_map="property"),
        )  # fmt: skip
        crss = [c.strip() for c in text_field(vtk, "features_crs").split(";")]
        assert crss == ["EPSG:25833", "EPSG:25833"], crss
