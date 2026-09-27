"""`tin_engine.dem_input`: `DemRequest` in, `open_dem` out (increment 15a, R1, R11).

The declarative front of `docs/increments/15-dem-mosaic.md`: the CLI parses
flags into a `DemRequest` and calls `open_dem`, and a GUI or API worker builds
the same request without Typer. This file tests it through files, synthetic
and real:

- micro-TIFF mosaics written to `tmp_path`;
- two real extracts of Ola's DTM10 archive (2022-09-24 release), cut by
  `tests/fixtures/dtm10/extract.py`: `seam/` (6400_4 | 6400_1, a 51-column
  overlap that agrees bit for bit, N3) and `lattices/` (7707_1, half a cell
  east-west off the main lattice, over its aligned neighbour 7707_2, N2).

HOW THIS FILE GOES RED: `tin_engine.dem_input` and `tin_engine.mosaic` are
imported in fixtures; each test fails on its own with `ModuleNotFoundError`.
"""

from __future__ import annotations

import importlib
import io
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import tifffile

from geotiff_fixtures import (
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    micro_tiff,
    truncated_strip,
    with_keys,
)
from mosaic_fixtures import DX, X0, Y0, quadrants, same_array, whole
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import DemTile

DTM10 = Path(__file__).resolve().parents[1] / "fixtures" / "dtm10"
SEAM = DTM10 / "seam"
LATTICES = DTM10 / "lattices"
WEST, EAST = SEAM / "6400_4_10m_z33.tif", SEAM / "6400_1_10m_z33.tif"
ODD, MAIN = LATTICES / "7707_1_10m_z33.tif", LATTICES / "7707_2_10m_z33.tif"


@pytest.fixture(scope="module")
def di() -> ModuleType:
    return importlib.import_module("tin_engine.dem_input")


@pytest.fixture(scope="module")
def mz() -> ModuleType:
    return importlib.import_module("tin_engine.mosaic")


def tiff_of(tile: DemTile, **kwargs: Any) -> io.BytesIO:
    """A point-registered micro-TIFF of `tile`: tie point at its first node."""
    m = tile.meta
    return micro_tiff(
        np.asarray(tile.array),
        tiepoint=(0.0, 0.0, 0.0, m.x_min, m.y_max, 0.0),
        scale=(m.delta_x, m.delta_y, 0.0),
        **kwargs,
    )


def write_tiles(directory: Path, tiles: dict[str, DemTile]) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, tile in tiles.items():
        path = directory / name
        path.write_bytes(tiff_of(tile).getvalue())
        paths.append(path)
    return paths


def decoded(path: Path) -> DemTile:
    with path.open("rb") as stream:
        return decode_dem(stream)


def request(
    di: ModuleType,
    *sources: Path,
    box: tuple[float, ...] | None = None,
    nodata: float | None = None,
) -> Any:
    bounds = None
    if box is not None:
        mz = importlib.import_module("tin_engine.mosaic")
        bounds = mz.Bounds(x_min=box[0], y_min=box[1], x_max=box[2], y_max=box[3])
    return di.DemRequest(sources=tuple(sources), bounds=bounds, nodata=nodata)


@pytest.fixture
def source() -> DemTile:
    return whole()  # 9 x 13, point-registered


@pytest.fixture
def mosaic_dir(tmp_path: Path, source: DemTile) -> Path:
    write_tiles(tmp_path / "tiles", quadrants(source, row_cut=4, col_cut=6, overlap=1))
    return tmp_path / "tiles"


class TestOpenDem:
    def test_one_file_is_the_file(self, di: ModuleType, tmp_path: Path) -> None:
        path = tmp_path / "one.tif"
        path.write_bytes(micro_tiff().getvalue())
        opened = di.open_dem(request(di, path))
        expected = decoded(path)
        assert opened.tile.meta == expected.meta
        assert same_array(opened.tile.array, expected.array)
        assert [p.name for p in opened.plan.tiles] == ["one.tif"]
        assert opened.label == "one"

    def test_a_directory_is_stitched(
        self, di: ModuleType, mosaic_dir: Path, source: DemTile
    ) -> None:
        opened = di.open_dem(request(di, mosaic_dir))
        assert opened.tile.meta == source.meta
        assert same_array(opened.tile.array, source.array)
        assert [p.name for p in opened.plan.tiles] == ["ne.tif", "nw.tif", "se.tif", "sw.tif"]
        assert opened.label == "tiles"

    def test_explicit_files_equal_the_directory(self, di: ModuleType, mosaic_dir: Path) -> None:
        files = sorted(mosaic_dir.glob("*.tif"), reverse=True)
        by_files = di.open_dem(request(di, *files))
        by_dir = di.open_dem(request(di, mosaic_dir))
        assert by_files.tile.meta == by_dir.tile.meta
        assert same_array(by_files.tile.array, by_dir.tile.array)
        assert by_files.label == files[0].stem

    def test_bounds_cut_the_window(self, di: ModuleType, mosaic_dir: Path, source: DemTile) -> None:
        opened = di.open_dem(
            request(di, mosaic_dir, box=(500012.0, 6599987.0, 500033.0, 6599996.0))
        )
        assert (opened.tile.meta.x_min, opened.tile.meta.y_max) == (X0 + DX, Y0)
        assert same_array(opened.tile.array, source.array[0:4, 1:5])

    def test_a_caller_nodata_applies_to_every_tile(self, di: ModuleType, mosaic_dir: Path) -> None:
        opened = di.open_dem(request(di, mosaic_dir, nodata=-9999.0))
        assert (opened.tile.meta.nodata, opened.tile.meta.nodata_source) == (-9999.0, "caller")

    def test_the_cap_refuses_before_any_pixel_is_read(
        self, di: ModuleType, mz: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """M5 through a file: its strip is cut, so a load would raise
        `GeoTiffError`; the refusal must be the cap's `MosaicError` instead."""
        path = tmp_path / "cut.tif"
        path.write_bytes(truncated_strip().getvalue())
        monkeypatch.setattr(mz, "physical_memory", lambda: 8)  # cap 4 bytes: 3 x 4 nodes exceed it
        with pytest.raises(mz.MosaicError, match="--bbox"):
            di.open_dem(request(di, path))


class TestDemRequest:
    """R11: exactly one directory, or one or more files."""

    def test_a_directory_with_files_is_refused(self, di: ModuleType, mosaic_dir: Path) -> None:
        with pytest.raises(ValueError, match="director"):
            di.open_dem(request(di, mosaic_dir, mosaic_dir / "nw.tif"))

    def test_two_directories_are_refused(
        self, di: ModuleType, mosaic_dir: Path, tmp_path: Path
    ) -> None:
        other = tmp_path / "other"
        other.mkdir()
        with pytest.raises(ValueError, match="director"):
            di.open_dem(request(di, mosaic_dir, other))

    def test_no_source_is_refused(self, di: ModuleType) -> None:
        with pytest.raises(ValueError):
            di.open_dem(request(di))

    def test_it_is_frozen(self, di: ModuleType, mosaic_dir: Path) -> None:
        req = request(di, mosaic_dir)
        with pytest.raises(ValueError):
            req.nodata = 1.0


def _shifted_copy(path: Path, target: Path, metres: float) -> Path:
    """`path` rewritten with its tie point `metres` further east; same values."""
    with tifffile.TiffFile(path) as tif:
        page = tif.pages.first
        tie, scale, array = page.tags[33922].value, page.tags[33550].value, page.asarray()
    target.parent.mkdir(parents=True, exist_ok=True)
    stream = micro_tiff(
        array,
        tiepoint=(0.0, 0.0, 0.0, tie[3] + metres, tie[4], 0.0),
        scale=tuple(scale),
        geokeys=with_keys({GT_RASTER_TYPE: PIXEL_IS_AREA}),
        nodata="-32767",
        compression="deflate",
    )
    target.write_bytes(stream.getvalue())
    return target


class TestRealDtm10:
    """Q1 and Q5 on real data, one release (N5)."""

    def test_the_extracts_are_where_the_design_says(self) -> None:
        west, east = decoded(WEST), decoded(EAST)
        assert (west.meta.epsg, west.meta.pixel_is_area, west.meta.nodata) == (
            25833,
            True,
            -32767.0,
        )
        assert east.meta.x_min == west.meta.x_min + (307 - 51) * 10.0  # a 51-column overlap
        assert np.array_equal(west.array[:, -51:], east.array[:, :51])
        odd, main = decoded(ODD), decoded(MAIN)
        assert (odd.meta.x_min - main.meta.x_min) / 10.0 == -0.5  # N2: 5 m east-west

    def test_q1_the_agreeing_seam_is_accepted(self, di: ModuleType) -> None:
        west, east = decoded(WEST), decoded(EAST)
        opened = di.open_dem(request(di, SEAM))
        assert (opened.tile.meta.rows, opened.tile.meta.cols) == (256, 307 + 307 - 51)
        assert opened.tile.meta.x_min == west.meta.x_min
        assert same_array(opened.tile.array[:, :307], west.array)
        assert same_array(opened.tile.array[:, 256:], east.array)

    def test_q1_the_seam_shifted_by_one_cell_is_refused(
        self, di: ModuleType, mz: ModuleType, tmp_path: Path
    ) -> None:
        """The control (N3): one column off, the real overlap disagrees."""
        (tmp_path / "shifted").mkdir()
        (tmp_path / "shifted" / WEST.name).write_bytes(WEST.read_bytes())
        _shifted_copy(EAST, tmp_path / "shifted" / EAST.name, 10.0)
        with pytest.raises(mz.MosaicError) as info:
            di.open_dem(request(di, tmp_path / "shifted"))
        message = str(info.value)
        assert WEST.name in message and EAST.name in message

    def test_q5_the_two_lattices_together_are_refused_naming_both(
        self, di: ModuleType, mz: ModuleType
    ) -> None:
        with pytest.raises(mz.MosaicError) as info:
            di.open_dem(request(di, LATTICES))
        message = str(info.value)
        for token in (ODD.name, MAIN.name, "0.5 cell", "east-west"):
            assert token in message, f"{token!r} not in {message!r}"

    @pytest.mark.parametrize("path", [ODD, MAIN], ids=["7707_1", "7707_2"])
    def test_q5_each_lattice_alone_opens(self, di: ModuleType, path: Path) -> None:
        opened = di.open_dem(request(di, path))
        assert opened.tile.meta == decoded(path).meta

    def test_q5_reading_a_box_on_7707_2_reaching_into_7707_1_opens_on_7707_2(
        self, di: ModuleType
    ) -> None:
        """Ola's Q5 reading on real extracts. 7707_1 (y 7750690..7749740) lies
        over 7707_2 (y 7750250..7749300); this box reaches 60 m into their
        overlap strip, so it selects both, and only 7707_2 covers every node."""
        box = (769800.0, 7749350.0, 770000.0, 7749800.0)
        opened = di.open_dem(request(di, LATTICES, box=box))
        assert [p.name for p in opened.plan.tiles] == [MAIN.name]
        main = decoded(MAIN)
        assert opened.tile.meta.x_min == 769800.0
        assert opened.tile.meta.y_max == 7749800.0
        # rows (7750250 - 7749800) / 10 = 45 .. 90, columns 5 .. 25 of 7707_2
        assert same_array(opened.tile.array, main.array[45:91, 5:26])

    def test_q5_reading_a_box_inside_the_strip_takes_the_tie_by_name(self, di: ModuleType) -> None:
        """Both lattices cover a box wholly inside the strip and hold one tile
        each here, so the tie goes to the first tile's name: 7707_1."""
        box = (769800.0, 7749800.0, 770000.0, 7750200.0)
        opened = di.open_dem(request(di, LATTICES, box=box))
        assert [p.name for p in opened.plan.tiles] == [ODD.name]
        assert opened.tile.meta.x_min % 10.0 == 5.0  # on 7707_1's lattice

    def test_q5_a_request_inside_one_lattice_opens_with_both_listed(self, di: ModuleType) -> None:
        """South of 7707_1's last row (y 7749740), only 7707_2 is selected."""
        opened = di.open_dem(request(di, LATTICES, box=(769800.0, 7749350.0, 770000.0, 7749700.0)))
        assert [p.name for p in opened.plan.tiles] == [MAIN.name]
