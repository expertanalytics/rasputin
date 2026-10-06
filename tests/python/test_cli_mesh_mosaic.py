"""`rasputin mesh --dem DIR | --dem F... [--bbox ...]`: increment 15a through the CLI.

`docs/increments/15-dem-mosaic.md` R6 and R11, and the parked design's C1-C6
and T-real, carried. The oracle for equivalence is the single-file path, which
the existing `test_cli_mesh_dem.py` pins and this increment must not move
(I7): a grid cut into tiles must mesh exactly like the grid as one file.

HOW THIS FILE GOES RED. It imports nothing new: until 15a lands, `--dem DIR`
is refused as "not a file" and `--bbox` is "No such option", so every test here
fails on its assertions, not at collection.

Micro-TIFFs are written to `tmp_path`. Two real extracts of Ola's DTM10
archive are read from `tests/fixtures/dtm10/` (see `extract.py` there), and
T-real cuts the committed benchmark tile into quadrants (`needs_codecs`).

Increment 25 (`docs/increments/25-plain-output.md`, D2 and D4): the mesh file
names the files used in `dem_source` (`; `-joined); `dem_tiles`, `dem_seams`
and `dem_grid` are `--stats` rows, read here through `inputs`, and stderr's
mosaic line reads `DEM: <n> files, <c> columns x <r> rows`. A pair's
`dem_seams` entry reads `<a> and <b> disagree at <n> node(s), by up to ...`; agreeing overlaps read
`none: the tiles agree where they overlap`.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np
import pytest
import tifffile

from cli_helpers import USAGE, invoke, runner
from geotiff_fixtures import (
    EPSG_UTM33,
    GT_RASTER_TYPE,
    KARTVERKET,
    METRE,
    PIXEL_IS_AREA,
    PROJECTED_CS_TYPE,
    VERTICAL_UNITS,
    micro_tiff,
    needs_codecs,
    with_keys,
)
from mosaic_fixtures import piece, quadrants, whole
from test_cli_mesh_refine import SEAMS_AGREE, file_field, stats_row
from tin_engine.cli import app
from tin_engine.io.models import DemTile
from vtkread import VtkFile, read_vtk

DTM10 = Path(__file__).resolve().parents[1] / "fixtures" / "dtm10"


def terrain(rows: int, cols: int) -> np.ndarray:
    """Not a plane, so `--tolerance` has work to do; exact in float32."""
    r, c = np.indices((rows, cols))
    return np.asarray(((7 * r + 3 * c) % 11) * 1.5 + 0.25 * r * c).astype(np.float32)


def tiff_of(tile: DemTile, compression: str | None = None) -> io.BytesIO:
    """`tile` as a GeoTIFF with the registration its meta says."""
    m = tile.meta
    half = 0.5 if m.pixel_is_area else 0.0
    keys = with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA}) if m.pixel_is_area else None
    return micro_tiff(
        np.asarray(tile.array),
        tiepoint=(0.0, 0.0, 0.0, m.x_min - half * m.delta_x, m.y_max + half * m.delta_y, 0.0),
        scale=(m.delta_x, m.delta_y, 0.0),
        geokeys=keys,
        nodata=None if m.nodata is None else f"{m.nodata:g}",
        compression=compression,
    )


def write_tiles(directory: Path, tiles: dict[str, DemTile]) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    for name, tile in tiles.items():
        (directory / name).write_bytes(tiff_of(tile).getvalue())
    return sorted(directory.iterdir())


def report(tmp_path: Path, *args: str) -> str:
    """`--stats -`'s Markdown, read from stdout apart from stderr and unflattened."""
    result = runner.invoke(app, ["mesh", *args, "--out", str(tmp_path / "m.vtk"), "--stats", "-"])
    assert result.exit_code == 0, result.output
    return result.stdout


def write_one(path: Path, tile: DemTile) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(tiff_of(tile).getvalue())
    return path


def run_vtk(out: Path, *args: str) -> VtkFile:
    """Mesh to ``out`` with ``--stats`` beside it; ``inputs`` reads that report."""
    code, output = invoke("mesh", *args, "--out", str(out), "--stats", str(out.with_suffix(".md")))
    assert code == 0, output
    return read_vtk(out.read_bytes())


def field(vtk: VtkFile, name: str) -> str:
    return file_field(vtk, name)


def inputs(out: Path, name: str) -> str:
    """The ``--stats`` row ``name`` of the ``run_vtk`` that wrote ``out``."""
    return stats_row(out.with_suffix(".md").read_text(encoding="utf-8"), name)


def grid_of(out: Path) -> tuple[int, int]:
    """(columns, rows) from ``dem_grid``: ``<c> columns x <r> rows, ...``."""
    match = re.match(r"(\d+) columns x (\d+) rows\b", inputs(out, "dem_grid"))
    assert match is not None, inputs(out, "dem_grid")
    return int(match.group(1)), int(match.group(2))


def same_mesh(a: VtkFile, b: VtkFile) -> None:
    assert np.array_equal(a.points, b.points)
    assert len(a.polygons) == len(b.polygons)
    assert all(np.array_equal(p, q) for p, q in zip(a.polygons, b.polygons, strict=True))
    assert len(a.lines) == len(b.lines)
    assert all(np.array_equal(p, q) for p, q in zip(a.lines, b.lines, strict=True))
    assert a.scalars.keys() == b.scalars.keys()
    for name in a.scalars:
        assert np.array_equal(
            np.asarray(a.scalars[name].values), np.asarray(b.scalars[name].values)
        ), name


@pytest.fixture
def source() -> DemTile:
    return whole(9, 13, array=terrain(9, 13))


@pytest.fixture
def mosaic_dir(tmp_path: Path, source: DemTile) -> Path:
    write_tiles(tmp_path / "tiles", quadrants(source, row_cut=4, col_cut=6, overlap=1))
    return tmp_path / "tiles"


class TestC1Fields:
    """C1 and R11: what a mosaic's `.vtk` records, and the stderr line."""

    def test_a_directory_writes_the_mosaic_fields(self, tmp_path: Path, mosaic_dir: Path) -> None:
        out = tmp_path / "m.vtk"
        code, output = invoke("mesh", "--dem", str(mosaic_dir), "--out", str(out))
        assert code == 0, output
        vtk = read_vtk(out.read_bytes())
        assert field(vtk, "crs") == "EPSG:25833"
        assert field(vtk, "dem_source") == "ne.tif; nw.tif; se.tif; sw.tif"
        for moved in ("elevation_source", "dem_tiles", "dem_seams"):
            assert moved not in vtk.field_data, moved
        assert "DEM: 4 files, 13 columns x 9 rows" in output
        assert "mosaic of" not in output

    def test_the_stats_rows(self, tmp_path: Path, mosaic_dir: Path) -> None:
        out = tmp_path / "m.vtk"
        run_vtk(out, "--dem", str(mosaic_dir))
        assert inputs(out, "dem_tiles") == "ne.tif; nw.tif; se.tif; sw.tif"
        assert grid_of(out) == (13, 9)

    def test_one_file_records_no_tile_list(self, tmp_path: Path, source: DemTile) -> None:
        out = tmp_path / "one.vtk"
        vtk = run_vtk(out, "--dem", str(write_one(tmp_path / "one.tif", source)))
        assert "dem_tiles" not in vtk.field_data
        assert field(vtk, "dem_source") == "one.tif"
        assert grid_of(out) == (13, 9)


class TestC2Equivalence:
    """C2 (T-equiv) and I5, I7: a tile split into four meshes exactly like the tile."""

    @pytest.mark.parametrize(
        ("area", "overlap", "shape"),
        [(False, 1, (9, 13)), (True, 0, (8, 12)), (True, 3, (8, 12))],
        ids=["point-shared-line", "area-abutting", "area-overlap-3"],
    )
    @pytest.mark.parametrize("extra", [(), ("--tolerance", "0.5")], ids=["stride", "tolerance"])
    def test_split_meshes_like_whole(
        self,
        tmp_path: Path,
        area: bool,
        overlap: int,
        shape: tuple[int, int],
        extra: tuple[str, ...],
    ) -> None:
        grid = whole(*shape, array=terrain(*shape), area=area)
        single = write_one(tmp_path / "whole.tif", grid)
        tiles = write_tiles(
            tmp_path / "split", quadrants(grid, row_cut=4, col_cut=6, overlap=overlap)
        )
        assert len(tiles) == 4
        expected = run_vtk(tmp_path / "whole.vtk", "--dem", str(single), *extra)
        got = run_vtk(tmp_path / "split.vtk", "--dem", str(tmp_path / "split"), *extra)
        same_mesh(got, expected)


class TestC3RepeatedDem:
    """C3: `--dem a --dem b ...` writes the same mesh as `--dem DIR`."""

    def test_files_equal_directory(self, tmp_path: Path, mosaic_dir: Path) -> None:
        files = [str(p) for p in sorted(mosaic_dir.glob("*.tif"), reverse=True)]
        by_dir = run_vtk(tmp_path / "dir.vtk", "--dem", str(mosaic_dir))
        by_files = run_vtk(tmp_path / "files.vtk", *[arg for f in files for arg in ("--dem", f)])
        same_mesh(by_files, by_dir)
        assert field(by_files, "dem_source") == field(by_dir, "dem_source")
        assert inputs(tmp_path / "files.vtk", "dem_tiles") == inputs(
            tmp_path / "dir.vtk", "dem_tiles"
        )


class TestC4Refusals:
    """C4: mosaic refusals are usage errors naming the files; nothing is written."""

    def refused(self, tmp_path: Path, *args: str, says: tuple[str, ...]) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke("mesh", *args, "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output, output
        assert "Traceback" not in output
        for word in says:
            assert word in output, f"{word!r} not in {output!r}"
        assert not out.exists()

    def test_mixed_crs(self, tmp_path: Path, source: DemTile) -> None:
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        directory = tmp_path / "mixed"
        write_tiles(directory, tiles)
        ne = tiles["ne.tif"]
        m = ne.meta
        stream = micro_tiff(
            np.asarray(ne.array),
            tiepoint=(0.0, 0.0, 0.0, m.x_min, m.y_max, 0.0),
            scale=(m.delta_x, m.delta_y, 0.0),
            geokeys=with_keys({PROJECTED_CS_TYPE: 25832}),
        )
        (directory / "ne.tif").write_bytes(stream.getvalue())
        self.refused(tmp_path, "--dem", str(directory), says=("ne.tif", "25832", str(EPSG_UTM33)))

    def test_half_a_cell_off(self, tmp_path: Path, source: DemTile) -> None:
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        tiles["se.tif"] = piece(source, 4, 9, 6, 13, x_min=tiles["se.tif"].meta.x_min + 5.0)
        write_tiles(tmp_path / "shifted", tiles)
        self.refused(tmp_path, "--dem", str(tmp_path / "shifted"), says=("se.tif", "0.5 cell"))

    def test_a_hole_in_the_tile_set(self, tmp_path: Path) -> None:
        from mosaic_fixtures import blocks

        tiles = blocks(whole(6, 8, array=terrain(6, 8), area=True), 3, 4, skip={(1, 1)})
        write_tiles(tmp_path / "holed", tiles)
        self.refused(tmp_path, "--dem", str(tmp_path / "holed"), says=("500040",))

    def test_an_empty_directory_names_it(self, tmp_path: Path) -> None:
        (tmp_path / "empty_dir").mkdir()
        self.refused(tmp_path, "--dem", str(tmp_path / "empty_dir"), says=("empty_dir",))


class TestQ1Seams:
    """Ola's Q1 revised (2026-09-28), through the CLI: a disagreeing overlap
    meshes, and each disagreeing pair is recorded in the `--stats` row
    `dem_seams` (increment 25 moved it from the mesh file) and in a
    "DEM seams" section.

    `dem_seams` is recorded whenever `dem_tiles` is: `none: the tiles agree
    where they overlap` when every overlapping node agrees, else one entry per
    disagreeing pair, sorted, `; `-joined, escaped like `dem_tiles`:
    `<first> and <second> disagree at <n> node(s), by up to <largest:g> m
    (median <median:g> m)` (increment 25, "Settled after the red step", 8).
    """

    @staticmethod
    def disagreeing(tmp_path: Path, source: DemTile, name: str = "nw.tif", by: float = 4.0) -> Path:
        """The quadrants, `ne.tif` `by` higher at one node of the shared column
        (global (2, 6), in `ne.tif` and the tile named `name` only). Was
        `TestC4Refusals.test_an_overlap_disagreement`'s directory."""
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        tiles[name] = tiles.pop("nw.tif")
        changed = np.array(tiles["ne.tif"].array)
        changed[2, 0] += np.float32(by)
        tiles["ne.tif"] = piece(source, 0, 5, 6, 13, array=changed)
        write_tiles(tmp_path / "disagree", tiles)
        return tmp_path / "disagree"

    def test_an_overlap_disagreement_meshes_and_is_recorded(
        self, tmp_path: Path, source: DemTile
    ) -> None:
        """Was `TestC4Refusals.test_an_overlap_disagreement`, which expected a
        refusal naming both tiles and the 4. The same names and number are now
        in the record, and the mesh is written."""
        out = tmp_path / "m.vtk"
        vtk = run_vtk(out, "--dem", str(self.disagreeing(tmp_path, source)))
        assert (
            inputs(out, "dem_seams")
            == "ne.tif and nw.tif disagree at 1 node, by up to 4 m (median 4 m)"
        )
        assert inputs(out, "dem_tiles") == "ne.tif; nw.tif; se.tif; sw.tif"
        assert "dem_seams" not in vtk.field_data

    def test_agreeing_tiles_record_none(self, tmp_path: Path, mosaic_dir: Path) -> None:
        run_vtk(tmp_path / "m.vtk", "--dem", str(mosaic_dir))
        assert inputs(tmp_path / "m.vtk", "dem_seams") == SEAMS_AGREE

    def test_one_file_records_no_seams(self, tmp_path: Path, source: DemTile) -> None:
        out = tmp_path / "one.vtk"
        vtk = run_vtk(out, "--dem", str(write_one(tmp_path / "one.tif", source)))
        assert "dem_seams" not in vtk.field_data
        report = out.with_suffix(".md").read_text(encoding="utf-8")
        assert "| dem_seams |" not in report and "`dem_seams`" not in report

    def test_several_pairs_are_sorted_and_joined(self, tmp_path: Path) -> None:
        """Each quadrant 1000 times its rank higher than the grid: every
        overlap disagrees, by a constant. Six pairs, sorted by name; counts are
        the overlaps' node counts (overlap 1 on a 9 x 13 grid: 5, 5, 7, 7, and
        the corner node in the two diagonal pairs)."""
        grid = whole(9, 13, array=terrain(9, 13))
        tiles = quadrants(grid, row_cut=4, col_cut=6, overlap=1)
        for rank, name in enumerate(sorted(tiles), start=1):
            t = tiles[name]
            tiles[name] = DemTile(meta=t.meta, array=np.asarray(t.array) + 1000 * rank)
        write_tiles(tmp_path / "all", tiles)
        run_vtk(tmp_path / "m.vtk", "--dem", str(tmp_path / "all"))
        assert inputs(tmp_path / "m.vtk", "dem_seams") == "; ".join(
            [
                "ne.tif and nw.tif disagree at 5 nodes, by up to 1000 m (median 1000 m)",
                "ne.tif and se.tif disagree at 7 nodes, by up to 2000 m (median 2000 m)",
                "ne.tif and sw.tif disagree at 1 node, by up to 3000 m (median 3000 m)",
                "nw.tif and se.tif disagree at 1 node, by up to 1000 m (median 1000 m)",
                "nw.tif and sw.tif disagree at 7 nodes, by up to 2000 m (median 2000 m)",
                "se.tif and sw.tif disagree at 5 nodes, by up to 1000 m (median 1000 m)",
            ]
        )

    def test_a_non_ascii_name_is_escaped(self, tmp_path: Path, source: DemTile) -> None:
        directory = self.disagreeing(tmp_path, source, name="Ålesund.tif")
        listed = next(p.name for p in directory.iterdir() if p.name.endswith("lesund.tif"))
        escaped = listed.encode("ascii", "backslashreplace").decode("ascii")
        run_vtk(tmp_path / "m.vtk", "--dem", str(directory))
        recorded = inputs(tmp_path / "m.vtk", "dem_seams")
        assert recorded.isascii()
        first, second = (
            n.encode("ascii", "backslashreplace").decode("ascii")
            for n in sorted([listed, "ne.tif"])
        )  # the file system may store the name decomposed, which sorts first
        assert recorded == f"{first} and {second} disagree at 1 node, by up to 4 m (median 4 m)"
        assert escaped in recorded

    def test_stats_has_a_seams_section_listing_only_disagreeing_pairs(
        self, tmp_path: Path, source: DemTile
    ) -> None:
        """After "Sizes", before "Quality": one table row per disagreeing pair."""
        lines = report(tmp_path, "--dem", str(self.disagreeing(tmp_path, source))).splitlines()
        at = lines.index("## DEM seams")
        assert lines.index("## Sizes") < at < lines.index("## Quality (plan view, x/y)")
        rows = []
        for line in lines[at + 1 :]:
            if line.startswith("## "):
                break
            if line.startswith("|"):
                rows.append([cell.strip() for cell in line.strip("|").split("|")])
        assert rows == [
            ["tile", "tile", "nodes", "max", "median"],
            ["---", "---", "---", "---", "---"],
            ["ne.tif", "nw.tif", "1", "4", "4"],
        ]

    def test_stats_has_no_seams_section_when_every_overlap_agrees(
        self, tmp_path: Path, mosaic_dir: Path
    ) -> None:
        lines = report(tmp_path, "--dem", str(mosaic_dir)).splitlines()
        assert "## Sizes" in lines
        assert "## DEM seams" not in lines

    def test_a_sub_millimetre_disagreement_records_none_and_no_section(
        self, tmp_path: Path, source: DemTile
    ) -> None:
        """Ola, 2026-09-28: "Ignore below 1mm". The one differing node is
        0.5 mm off, so no pair qualifies: `none`, and `--stats` has no
        "DEM seams" section."""
        directory = self.disagreeing(tmp_path, source, by=0.0005)
        run_vtk(tmp_path / "m.vtk", "--dem", str(directory))
        assert inputs(tmp_path / "m.vtk", "dem_seams") == SEAMS_AGREE
        lines = report(tmp_path, "--dem", str(directory)).splitlines()
        assert "## Sizes" in lines
        assert "## DEM seams" not in lines

    def test_a_pipe_in_a_tile_name_is_escaped_in_the_stats_table(
        self, tmp_path: Path, source: DemTile
    ) -> None:
        """Review suggestion on 15a: a `|` in a tile name would end its
        Markdown table cell early. In the "DEM seams" table it is written
        `\\|`, so the row still has five cells (split on unescaped pipes).
        The `dem_seams` row is escaped the same way, and reads, unescaped, the
        name as listed (`ne.tif` sorts first: `e` < `|`)."""
        directory = self.disagreeing(tmp_path, source, name="n|w.tif")
        lines = report(tmp_path, "--dem", str(directory)).splitlines()
        at = lines.index("## DEM seams")
        row = next(line for line in lines[at + 1 :] if line.startswith("| ne.tif"))
        cells = [cell.strip() for cell in re.split(r"(?<!\\)\|", row.strip().strip("|"))]
        assert cells == ["ne.tif", "n\\|w.tif", "1", "4", "4"]
        run_vtk(tmp_path / "m.vtk", "--dem", str(directory))
        seams = inputs(tmp_path / "m.vtk", "dem_seams")
        assert seams == "ne.tif and n|w.tif disagree at 1 node, by up to 4 m (median 4 m)"


class TestC5Usage:
    """C5 and R11: `--dem` forms, and `--bbox`."""

    refused = TestC4Refusals.refused

    def test_a_directory_with_files(self, tmp_path: Path, mosaic_dir: Path) -> None:
        self.refused(
            tmp_path, "--dem", str(mosaic_dir), "--dem", str(mosaic_dir / "nw.tif"), says=("--dem",)
        )

    def test_two_directories(self, tmp_path: Path, mosaic_dir: Path) -> None:
        other = tmp_path / "other"
        write_tiles(other, {"a.tif": whole(3, 4, array=terrain(3, 4))})
        self.refused(tmp_path, "--dem", str(mosaic_dir), "--dem", str(other), says=("--dem",))

    def test_bbox_without_dem(self, tmp_path: Path) -> None:
        self.refused(tmp_path, "catchment", "--bbox", "0", "0", "1", "1", says=("--bbox",))

    @pytest.mark.parametrize(
        "box", [("1", "0", "0", "1"), ("0", "0", "nan", "1")], ids=["inverted", "nan"]
    )  # fmt: skip
    def test_an_invalid_box(self, tmp_path: Path, mosaic_dir: Path, box: tuple[str, ...]) -> None:
        self.refused(tmp_path, "--dem", str(mosaic_dir), "--bbox", *box, says=("--bbox",))

    def test_bbox_restricts_the_mesh_to_the_snapped_window(
        self, tmp_path: Path, mosaic_dir: Path
    ) -> None:
        # cols floor(1.2)=1 .. ceil(8.3)=9 (x 500010..500090); rows floor(0.8)=0 .. ceil(6.6)=7
        vtk = run_vtk(
            tmp_path / "box.vtk",
            "--dem",
            str(mosaic_dir),
            "--bbox",
            "500012",
            "6599967",
            "500083",
            "6599996",
        )
        xs, ys = vtk.points[:, 0], vtk.points[:, 1]
        assert (xs.min(), xs.max()) == (500010.0, 500090.0)
        assert (ys.min(), ys.max()) == (6599965.0, 6600000.0)
        assert grid_of(tmp_path / "box.vtk") == (9, 8)
        assert field(vtk, "dem_source") == "ne.tif; nw.tif; se.tif; sw.tif"

    def test_a_box_reaching_into_another_lattices_strip_meshes_on_the_covering_one(
        self, tmp_path: Path, mosaic_dir: Path, source: DemTile
    ) -> None:
        """Ola's Q5 reading, end to end. `odd.tif` is half a cell east-west off
        the quadrants and overlaps `se.tif`; the box selects both, but only the
        main lattice covers it, so it meshes there and records `se.tif` alone."""
        odd = whole(5, 7, array=terrain(5, 7), x_min=500105.0, y_max=6599970.0)
        write_tiles(mosaic_dir, {"odd.tif": odd})
        # cols floor(7.2)=7 .. ceil(11.2)=12; rows floor(5.2)=5 .. ceil(7.6)=8
        vtk = run_vtk(
            tmp_path / "box.vtk",
            "--dem", str(mosaic_dir),
            "--bbox", "500072", "6599962", "500112", "6599974",
        )  # fmt: skip
        assert field(vtk, "dem_source") == "se.tif"
        assert inputs(tmp_path / "box.vtk", "dem_tiles") == "se.tif"
        xs, ys = vtk.points[:, 0], vtk.points[:, 1]
        assert (xs.min(), xs.max()) == (500070.0, 500120.0)
        assert (ys.min(), ys.max()) == (6599960.0, 6599975.0)

    def test_bbox_on_one_file(self, tmp_path: Path, source: DemTile) -> None:
        """Q2: `--bbox` applies to a single file too; it is then a window of it."""
        single = write_one(tmp_path / "one.tif", source)
        vtk = run_vtk(
            tmp_path / "box.vtk",
            "--dem",
            str(single),
            "--bbox",
            "500012",
            "6599967",
            "500083",
            "6599996",
        )
        assert (vtk.points[:, 0].min(), vtk.points[:, 0].max()) == (500010.0, 500090.0)


class TestC6NonAscii:
    """C6: a non-ASCII tile name is recorded escaped, and the file is written."""

    def test_the_name_is_backslash_escaped(self, tmp_path: Path, source: DemTile) -> None:
        tiles = quadrants(source, row_cut=4, col_cut=6, overlap=1)
        tiles["Ålesund.tif"] = tiles.pop("nw.tif")
        write_tiles(tmp_path / "names", tiles)
        listed = next(
            p.name for p in (tmp_path / "names").iterdir() if p.name.endswith("lesund.tif")
        )
        vtk = run_vtk(tmp_path / "n.vtk", "--dem", str(tmp_path / "names"))
        escaped = listed.encode("ascii", "backslashreplace").decode("ascii")
        for recorded in (field(vtk, "dem_source"), inputs(tmp_path / "n.vtk", "dem_tiles")):
            assert recorded.isascii()
            assert escaped in recorded


class TestRealDtm10:
    """The 15a acceptance on real extracts of one release (N5): Q1 and Q5."""

    def test_the_real_seam_meshes(self, tmp_path: Path) -> None:
        out = tmp_path / "seam.vtk"
        vtk = run_vtk(out, "--dem", str(DTM10 / "seam"), "--tolerance", "1")
        assert field(vtk, "dem_source") == "6400_1_10m_z33.tif; 6400_4_10m_z33.tif"
        assert inputs(out, "dem_tiles") == "6400_1_10m_z33.tif; 6400_4_10m_z33.tif"
        assert inputs(out, "dem_grid") == "563 columns x 256 rows, 10 m apart"  # the design's
        assert inputs(out, "dem_seams") == SEAMS_AGREE  # Ola's Q1 revised: it agrees
        assert np.isfinite(vtk.points).all()

    def test_the_agreeing_real_seam_meshes_like_its_stitched_tile(self, tmp_path: Path) -> None:
        """Ola's Q1 revised changes nothing where the overlap agrees: the seam
        meshes exactly like one file holding west's 307 columns and east's
        last 256, with west's georeferencing (the mesh before the change)."""
        west_path = DTM10 / "seam" / "6400_4_10m_z33.tif"
        east_path = DTM10 / "seam" / "6400_1_10m_z33.tif"
        with tifffile.TiffFile(west_path) as tif:
            page = tif.pages.first
            tie, scale, west = page.tags[33922].value, page.tags[33550].value, page.asarray()
        east = tifffile.imread(east_path)
        assert np.array_equal(west[:, -51:], east[:, :51])
        stitched = micro_tiff(
            np.ascontiguousarray(np.hstack([west, east[:, 51:]])),
            tiepoint=(0.0, 0.0, 0.0, tie[3], tie[4], 0.0),
            scale=tuple(scale),
            geokeys=with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA}),
            nodata="-32767",
            compression="deflate",
        )
        single = tmp_path / "stitched.tif"
        single.write_bytes(stitched.getvalue())
        expected = run_vtk(tmp_path / "one.vtk", "--dem", str(single), "--tolerance", "1")
        got = run_vtk(tmp_path / "seam.vtk", "--dem", str(DTM10 / "seam"), "--tolerance", "1")
        same_mesh(got, expected)

    def test_the_two_lattices_are_refused_naming_both(self, tmp_path: Path) -> None:
        out = tmp_path / "x.vtk"
        code, output = invoke("mesh", "--dem", str(DTM10 / "lattices"), "--out", str(out))
        assert code == USAGE, output
        for token in ("7707_1_10m_z33.tif", "7707_2_10m_z33.tif", "0.5 cell"):
            assert token in output, f"{token!r} not in {output!r}"
        assert not out.exists()

    def test_a_box_inside_one_lattice_meshes(self, tmp_path: Path) -> None:
        vtk = run_vtk(
            tmp_path / "box.vtk",
            "--dem", str(DTM10 / "lattices"),
            "--bbox", "769800", "7749350", "770000", "7749700",
            "--tolerance", "1",
        )  # fmt: skip
        # One tile selected out of a directory: the file used is on record
        # (the label is the directory's).
        assert field(vtk, "dem_source") == "7707_2_10m_z33.tif"
        assert inputs(tmp_path / "box.vtk", "dem_tiles") == "7707_2_10m_z33.tif"


@needs_codecs
def test_t_real_the_benchmark_tile_in_quadrants_meshes_like_the_tile(tmp_path: Path) -> None:
    """T-real: the parked design's, stated honestly. It proves the stitch on
    real data, real NoData padding and real size. It does not prove that two
    real neighbours agree; `TestRealDtm10` and `test_dem_input.py` do that."""
    with tifffile.TiffFile(KARTVERKET) as tif:
        page = tif.pages.first
        tie, scale, array = page.tags[33922].value, page.tags[33550].value, page.asarray()
    rows, cols = array.shape
    row_cut, col_cut = rows // 2, cols // 2
    directory = tmp_path / "quadrants"
    directory.mkdir()
    for name, (r0, r1, c0, c1) in {
        "nw.tif": (0, row_cut + 1, 0, col_cut + 1),
        "ne.tif": (0, row_cut + 1, col_cut, cols),
        "sw.tif": (row_cut, rows, 0, col_cut + 1),
        "se.tif": (row_cut, rows, col_cut, cols),
    }.items():
        stream = micro_tiff(
            np.ascontiguousarray(array[r0:r1, c0:c1]),
            tiepoint=(0.0, 0.0, 0.0, tie[3] + c0 * scale[0], tie[4] - r0 * scale[1], 0.0),
            scale=tuple(scale),
            geokeys=with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA, VERTICAL_UNITS: METRE}),
            nodata="-32767",
            compression="deflate",
        )
        (directory / name).write_bytes(stream.getvalue())
    expected = run_vtk(tmp_path / "tile.vtk", "--dem", str(KARTVERKET), "--tolerance", "1")
    got = run_vtk(tmp_path / "split.vtk", "--dem", str(directory), "--tolerance", "1")
    same_mesh(got, expected)
