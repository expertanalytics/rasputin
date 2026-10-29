"""`rasputin mesh` with several sources of edges in one mesh: increment 16e,
as `docs/increments/111-one-land-cover-map.md` reshapes it.

`docs/increments/16e-multi-features.md`: R2/R3 (the degeneracy policy), R4
(the noder merges cross-source coincident edges, unioning the masks), R5
(labelling spans the land-cover files), R6 (one entry per source, one
distinct notice), and "Tests, and the regression floor". Increment 111,
section 3 and 8 (`test_cli_mesh_multi_features`): `--features`' positional
pairing, two maps over one file and a codeless source beside a coded one
go; two coded files are two `--land-cover` files under one land-cover class
file, read as one map; an uncoded boundary is a `--lines` line; the shared
seam is a land-cover border, and a `--lines` line along it; two CORINE
files give one `land_cover_notice`.

Reuses `test_cli_mesh_features.py`'s fixtures and helpers (`bumpy`,
`plain_square`, `mesh`, `meshed`, `edges`, `text_field`, `report`,
`edges_entries`) so the multi-source cases read the same way the
single-source ones do.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon

from cli_helpers import SQUARE, USAGE, invoke, polygon_file, rough_dem
from feature_fixtures import Feat, write_geojson
from gpkg_fixtures import Layer, Row, write_gpkg
from landcover_fixtures import vtk_labels
from recordread import stats_row
from test_cli_mesh_features import (
    SNAP,
    edges,
    edges_entries,
    mesh,
    meshed,
    rel,
    report,
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
# A boundary that names no class: a wall crossing the middle of the domain,
# a `--lines wall` line that blocks the flood fill.
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


bumpy = rough_dem(16)
plain_square = polygon_file(SQUARE)


@pytest.fixture
def corine_src(tmp_path: Path) -> Path:
    """A CORINE-coded source: forest 311 and lake 512."""
    return write_geojson(
        tmp_path / "corine.geojson", [coded("forest", FOREST, "311"), coded("lake", LAKE, "512")]
    )


@pytest.fixture
def parcel_src(tmp_path: Path) -> Path:
    """A line source: one wall boundary."""
    return write_geojson(tmp_path / "parcel.geojson", [Feat("wall", WALL, {})])


def land_cover(*paths: Path, classes: str = "corine") -> list[str]:
    """Each of ``paths`` as a `--land-cover` file, under one class file."""
    given = [a for path in paths for a in ("--land-cover", str(path))]
    return [*given, "--land-cover-classes", classes]


def with_wall(cover: Path, wall: Path) -> list[str]:
    """`cover` as the land cover, `wall` as `--lines wall`, `--edges all`."""
    return [*land_cover(cover), "--lines", f"wall={wall}", "--edges", "all"]


# --------------------------------------------------------- R7 case 1: two sources


class TestTwoSources:
    """R7 case 1 / I3: a CORINE land cover and a `--lines wall` boundary in
    one mesh. The mesh carries both sets of constraint edges, the land
    cover's `land_cover_code`, and the wall's edges block the fill."""

    @pytest.fixture
    def vtk(self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path,
            parcel_src: Path) -> VtkFile:  # fmt: skip
        return meshed(tmp_path, bumpy, plain_square, *with_wall(corine_src, parcel_src))[0]

    def test_both_sources_constraints_are_present(self, vtk: VtkFile) -> None:
        _, masks = edges(vtk)
        found = set(masks.tolist())
        # CORINE: land_cover on the forest, land_cover+water on the lake.
        assert {V.mask("land_cover"), V.mask("land_cover", "water")} <= found
        # The wall carries its own `wall` bit, from the line source.
        assert V.mask("wall") in found

    def test_the_wall_edge_is_a_constraint_alongside_the_coded_edges(self, vtk: VtkFile) -> None:
        """The wall (the line source) and the CORINE edges (the map) are both
        present: neither source overrode the other."""
        xy, masks = edges(vtk)
        wall = xy[(masks & V.mask("wall")) != 0]
        assert len(wall), "the line source contributed no wall constraint"
        land = xy[(masks & V.mask("land_cover")) != 0]
        assert len(land), "the land cover contributed no land_cover constraint"

    def test_the_coded_source_labels_the_cells(self, vtk: VtkFile) -> None:
        """I3: the CORINE polygons label the cells; the wall adds no code but
        its line still partitions the fill."""
        _, _, codes = vtk_labels(vtk, [(FOREST, 311), (LAKE, 512)])
        found = set(codes.tolist())
        assert 311 in found and 512 in found

    def test_the_record_names_both_sources(self, vtk: VtkFile, tmp_path: Path) -> None:
        """R6 / I5: the `edges` row (a `--stats` row) has one entry per
        source, the map first; the mesh file has none."""
        assert "edges" not in vtk.field_data
        cover, wall = edges_entries(report(tmp_path))
        assert (cover["name"], cover["classes"]) == ("corine.geojson", "corine")
        assert (wall["name"], wall["key"]) == ("parcel.geojson", "wall")

    def test_one_land_cover_file_is_one_entry(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path
    ) -> None:
        """The regression guard: one `--land-cover` file is one entry."""
        meshed(tmp_path, bumpy, plain_square, *land_cover(corine_src), "--edges", "all")
        assert ";" not in stats_row(report(tmp_path), "edges")
        (entry,) = edges_entries(report(tmp_path))
        assert entry["classes"] == "corine"


# --------------------------- R7 case 4 / I2: cross-source edge merge (R4)


class TestCrossSourceMerge:
    """R4 / I2: two sources with an edge coincident in the DEM's CRS give one
    constraint edge whose mask is the union of the two sources' bits.

    LEFT is a lake (water) in one `--land-cover` file, RIGHT forest in a
    second: one map, so their shared side is one land-cover border, carrying
    `land_cover` and `water` (the lake's side). With `SEAM` also given as a
    `--lines road` line, the seam's edges carry the road's bit as well."""

    def seam_mesh(self, tmp_path: Path, bumpy: Path, plain_square: Path, *more: str) -> VtkFile:
        left = write_geojson(tmp_path / "left.geojson", [coded("l", LEFT, "512")])
        right = write_geojson(tmp_path / "right.geojson", [coded("r", RIGHT, "311")])
        return meshed(
            tmp_path, bumpy, plain_square, *land_cover(left, right), "--edges", "all", *more
        )[0]

    def seam_masks(self, vtk: VtkFile) -> tuple[np.ndarray, np.ndarray]:
        xy, masks = edges(vtk)
        mids = (xy[:, 0] + xy[:, 1]) / 2
        on_seam = np.asarray(shapely.distance(SEAM, shapely.points(mids))) <= SNAP
        assert on_seam.any(), "no edge found on the shared seam"
        return xy[on_seam], masks[on_seam]

    def test_the_shared_seam_carries_both_masks(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        _, masks = self.seam_masks(self.seam_mesh(tmp_path, bumpy, plain_square))
        want = V.mask("land_cover", "water")
        assert set(masks.tolist()) == {want}, sorted(set(masks.tolist()))

    def test_a_line_along_the_seam_merges_with_the_border(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        road = write_geojson(tmp_path / "seam.geojson", [Feat("seam", SEAM, {})])
        vtk = self.seam_mesh(tmp_path, bumpy, plain_square, "--lines", f"road={road}")
        _, masks = self.seam_masks(vtk)
        want = V.mask("land_cover", "water", "road")
        assert set(masks.tolist()) == {want}, sorted(set(masks.tolist()))

    @pytest.mark.parametrize("with_road", [False, True], ids=["border", "border-and-road"])
    def test_the_seam_is_one_run_of_edges_not_doubled(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, with_road: bool
    ) -> None:
        """The merge collapses the coincident chains to one edge run: no seam
        edge is doubled (no two edges share both endpoints; a doubled merge
        would repeat one). Read from the seam edges themselves, without
        guessing a vertex coordinate."""
        more: tuple[str, ...] = ()
        if with_road:
            road = write_geojson(tmp_path / "seam.geojson", [Feat("seam", SEAM, {})])
            more = ("--lines", f"road={road}")
        seam_xy, _ = self.seam_masks(self.seam_mesh(tmp_path, bumpy, plain_square, *more))
        keys = {tuple(sorted((tuple(np.round(a, 3)), tuple(np.round(b, 3))))) for a, b in seam_xy}
        assert len(keys) == len(seam_xy), "a seam edge is doubled (a parallel edge on the seam)"


# ------------------------------------------------------ R7 case 5: degeneracy (R3)


class TestDegeneracy:
    """R3: every degenerate case a clear usage error naming the option."""

    def test_edges_with_no_domain_is_refused(
        self, tmp_path: Path, bumpy: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        """R3: `--edges` with two sources and no `--domain`."""
        code, output = invoke(
            "mesh",
            "--dem", str(bumpy),
            *with_wall(corine_src, parcel_src),
            "--out", str(bumpy.parent / "x.vtk"),
        )  # fmt: skip
        assert code == USAGE and "--edges" in output, output
        assert "needs --dem, --domain and --tolerance" in output, output
        assert "No such option" not in output, output

    def test_a_layer_with_a_geojson_land_cover_file_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        """`--land-cover-layer` is for every `--land-cover` file: with a
        GeoPackage and a GeoJSON file it is refused for the GeoJSON one,
        naming the option and the file."""
        forest_rows = [Row(1, FOREST, {"Code_18": "311"})]
        gpkg = write_gpkg(
            tmp_path / "corine.gpkg",
            [Layer("corine2018", 25833, forest_rows, columns=("Code_18",))],
        )
        lake = write_geojson(tmp_path / "lake.geojson", [coded("l", LAKE, "512")])
        code, output, target = mesh(
            tmp_path, bumpy, plain_square,
            *land_cover(gpkg, lake), "--land-cover-layer", "corine2018", "--edges", "all",
        )  # fmt: skip
        assert code == USAGE, output
        assert "--land-cover-layer" in output and "No such option" not in output, output
        assert lake.name in output, output
        assert not target.exists()


# --------------------------------------------------- D1: two land-cover files


class TestTwoLandCoverFiles:
    """Two `--land-cover` files under one class file are one map: both
    files' polygons label the mesh."""

    def test_both_files_label_the_mesh(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        a = write_geojson(tmp_path / "a.geojson", [coded("f", FOREST, "311")])
        b = write_geojson(tmp_path / "b.geojson", [coded("l", LAKE, "512")])
        code, output, target = mesh(
            tmp_path, bumpy, plain_square, *land_cover(a, b), "--edges", "all"
        )
        assert code == 0, output
        vtk = read_vtk(target.read_bytes())
        assert "land_cover_codes" in vtk.field_data
        # Both files' polygons label the mesh: neither overrode the other.
        _, _, codes = vtk_labels(vtk, [(FOREST, 311), (LAKE, 512)])
        assert {311, 512} <= set(codes.tolist()), sorted(set(codes.tolist()))


# ------------------------------------------------------------ R6 / I5: the record


class TestRecordPerSource:
    """R6 / I5: the map is one `edges` entry naming its files, each line key
    one more; `land_cover_crs` has one entry per land-cover file;
    `land_cover_notice` lists each distinct notice once; the per-source
    `<n> features` counts are each source's own."""

    def test_two_corine_files_give_one_notice(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        a = write_geojson(tmp_path / "a.geojson", [coded("f", FOREST, "311")])
        b = write_geojson(tmp_path / "b.geojson", [coded("l", LAKE, "512")])
        vtk, _ = meshed(tmp_path, bumpy, plain_square, *land_cover(a, b), "--edges", "all")
        (entry,) = edges_entries(report(tmp_path))
        assert entry["name"] == "a.geojson + b.geojson", entry.string
        assert int(entry["n"]) == 2, entry.string
        notice = text_field(vtk, "land_cover_notice")
        assert notice.count("Copernicus Land Monitoring Service") == 1, notice

    def test_per_source_feature_counts_are_correct(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, corine_src: Path, parcel_src: Path
    ) -> None:
        """The map has two features, the wall one; each entry's `<n>
        features` count is that source's own, not the run total."""
        meshed(tmp_path, bumpy, plain_square, *with_wall(corine_src, parcel_src))
        counts = [int(m["n"]) for m in edges_entries(report(tmp_path))]
        assert counts == [2, 1], counts

    def test_land_cover_crs_has_one_entry_per_file(
        self, tmp_path: Path, bumpy: Path, plain_square: Path
    ) -> None:
        a = write_geojson(tmp_path / "a.geojson", [coded("f", FOREST, "311")])
        b = write_geojson(tmp_path / "b.geojson", [coded("l", LAKE, "512")])
        meshed(tmp_path, bumpy, plain_square, *land_cover(a, b), "--edges", "all")
        crss = [c.strip() for c in stats_row(report(tmp_path), "land_cover_crs").split(";")]
        assert crss == ["EPSG:25833", "EPSG:25833"], crss
