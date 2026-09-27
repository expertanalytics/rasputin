"""`tin_engine.io.repository`: the directory-backed DEM repository (increment 15a).

`docs/increments/15-dem-mosaic.md` R1-R3 and Ola's Q3 ruling: the repository
lives in `io/repository.py`, **the one module in `io/` that opens files**,
read-only. Tests F1-F4 are the parked design's, carried; the Q3 guard and
the laziness tests are new. Files are micro-TIFFs under `tmp_path`.

HOW THIS FILE GOES RED: the module is imported in a fixture, so each test
fails on its own with `ModuleNotFoundError` and collection is unaffected. The
Q3 guard's own scanner is tested on planted source, so it is shown able to
fail before it is trusted on the package.
"""

from __future__ import annotations

import ast
import importlib
import io
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

from geotiff_fixtures import (
    DELTA_X,
    TIE_X,
    TIE_Y,
    elevations,
    micro_tiff,
    truncated_strip,
)
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import GeoTiffError

IO = Path(__file__).resolve().parents[2] / "src_python" / "tin_engine" / "io"


@pytest.fixture(scope="module")
def repo() -> ModuleType:
    return importlib.import_module("tin_engine.io.repository")


def write(path: Path, stream: io.BytesIO | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((stream or micro_tiff()).getvalue())
    return path


def shifted(columns: int) -> io.BytesIO:
    """The baseline micro-TIFF moved `columns` nodes east, so names map to places."""
    return micro_tiff(tiepoint=(0.0, 0.0, 0.0, TIE_X + columns * DELTA_X, TIE_Y, 0.0))


def decoded(path: Path) -> object:
    with path.open("rb") as stream:
        return decode_dem(stream)


class TestF1Listing:
    """F1: `from_directory` lists `*.tif` and `*.tiff`, any case, not recursively; sorted."""

    def test_lists_tif_and_tiff_only_sorted(self, repo: ModuleType, tmp_path: Path) -> None:
        write(tmp_path / "b.TIFF", shifted(3))
        write(tmp_path / "a.tif", shifted(0))
        write(tmp_path / "c.Tif", shifted(6))
        (tmp_path / "a.tfw").write_text("10\n0\n0\n-10\n0\n0\n")
        (tmp_path / "a.tif.aux.xml").write_text("<PAMDataset/>")
        (tmp_path / "notes.txt").write_text("not a tile")
        write(tmp_path / "sub" / "d.tif", shifted(9))
        names = [f.name for f in repo.TiffDemRepository.from_directory(tmp_path).footprints()]
        assert names == ["a.tif", "b.TIFF", "c.Tif"]

    def test_a_footprint_is_the_files_name_and_header(
        self, repo: ModuleType, tmp_path: Path
    ) -> None:
        path = write(tmp_path / "7908_3_10m_z33.tif", shifted(0))
        (footprint,) = repo.TiffDemRepository.from_directory(tmp_path).footprints()
        assert footprint.name == "7908_3_10m_z33.tif"
        assert footprint.meta == decoded(path).meta  # type: ignore[attr-defined]

    def test_an_empty_directory_is_refused_naming_it(
        self, repo: ModuleType, tmp_path: Path
    ) -> None:
        empty = tmp_path / "no_tiles_here"
        empty.mkdir()
        (empty / "readme.txt").write_text("x")
        with pytest.raises(ValueError, match="no_tiles_here"):
            repo.TiffDemRepository.from_directory(empty).footprints()

    def test_a_subdirectory_named_like_a_tile_is_not_a_tile(
        self, repo: ModuleType, tmp_path: Path
    ) -> None:
        write(tmp_path / "a.tif")
        (tmp_path / "trap.tif").mkdir()
        names = [f.name for f in repo.TiffDemRepository.from_directory(tmp_path).footprints()]
        assert names == ["a.tif"]


class TestF2HeaderOnly:
    """F2 and I6: footprints read headers only, and read them once, lazily."""

    def test_a_file_with_truncated_pixels_lists_and_refuses_on_load(
        self, repo: ModuleType, tmp_path: Path
    ) -> None:
        """With the pixel read in the footprint path, listing would raise here."""
        write(tmp_path / "cut.tif", truncated_strip())
        repository = repo.TiffDemRepository.from_directory(tmp_path)
        (footprint,) = repository.footprints()
        assert (footprint.meta.rows, footprint.meta.cols) == (3, 4)
        with pytest.raises(GeoTiffError, match="pixel data"):
            repository.load("cut.tif")

    def test_construction_reads_nothing(self, repo: ModuleType, tmp_path: Path) -> None:
        write(tmp_path / "junk.tif", io.BytesIO(b"not a tiff"))
        repository = repo.TiffDemRepository.from_directory(tmp_path)  # no raise yet
        with pytest.raises(GeoTiffError):
            repository.footprints()

    def test_footprints_are_cached(self, repo: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "a.tif")
        repository = repo.TiffDemRepository.from_directory(tmp_path)
        first = repository.footprints()
        path.unlink()  # a second read of the header would now fail
        assert repository.footprints() == first

    def test_load_decodes_the_whole_tile(self, repo: ModuleType, tmp_path: Path) -> None:
        path = write(tmp_path / "a.tif", micro_tiff(elevations(rows=5, cols=7)))
        tile = repo.TiffDemRepository.from_directory(tmp_path).load("a.tif")
        expected = decoded(path)
        assert tile.meta == expected.meta  # type: ignore[attr-defined]
        assert np.array_equal(tile.array, expected.array)  # type: ignore[attr-defined]

    def test_load_of_an_unknown_name(self, repo: ModuleType, tmp_path: Path) -> None:
        write(tmp_path / "a.tif")
        with pytest.raises(KeyError):
            repo.TiffDemRepository.from_directory(tmp_path).load("b.tif")

    def test_a_caller_nodata_reaches_both_phases(self, repo: ModuleType, tmp_path: Path) -> None:
        """`DemRequest.nodata` (R1) is the caller asserting a sentinel for every tile."""
        write(tmp_path / "a.tif")
        repository = repo.TiffDemRepository.from_directory(tmp_path, nodata=-9999.0)
        (footprint,) = repository.footprints()
        assert (footprint.meta.nodata, footprint.meta.nodata_source) == (-9999.0, "caller")
        assert repository.load("a.tif").meta == footprint.meta


class TestF3Refusals:
    """F3: a listed file that fails `read_meta` refuses the listing, naming the file."""

    def test_a_tif_that_is_not_a_geotiff(self, repo: ModuleType, tmp_path: Path) -> None:
        write(tmp_path / "a.tif")
        write(tmp_path / "stray.tif", micro_tiff(geokeys={}))
        with pytest.raises(GeoTiffError, match=r"stray\.tif"):
            repo.TiffDemRepository.from_directory(tmp_path).footprints()


class TestF4ExplicitPaths:
    """F4: explicit paths, resolved and deduplicated, list like the directory."""

    def test_duplicates_collapse_and_the_listing_matches(
        self, repo: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        a = write(tmp_path / "a.tif", shifted(0))
        b = write(tmp_path / "b.tif", shifted(3))
        monkeypatch.chdir(tmp_path)
        explicit = repo.TiffDemRepository([b, a, Path("a.tif"), tmp_path / "." / "b.tif"])
        listed = repo.TiffDemRepository.from_directory(tmp_path)
        assert explicit.footprints() == listed.footprints()
        assert [f.name for f in explicit.footprints()] == ["a.tif", "b.tif"]

    def test_two_files_with_one_name_are_refused(self, repo: ModuleType, tmp_path: Path) -> None:
        """A footprint's name is the file name, so it must be unique."""
        one = write(tmp_path / "x" / "t.tif")
        two = write(tmp_path / "y" / "t.tif", shifted(3))
        with pytest.raises(ValueError, match=r"\bt\.tif"):
            repo.TiffDemRepository([one, two]).footprints()

    def test_it_satisfies_the_protocol(self, repo: ModuleType, tmp_path: Path) -> None:
        write(tmp_path / "a.tif")
        repository: object = repo.TiffDemRepository.from_directory(tmp_path)
        assert callable(getattr(repository, "footprints", None))
        assert callable(getattr(repository, "load", None))
        assert set(repo.DemRepository.__dict__) >= {"footprints", "load"}


# ---------------------------------------------------------------------------
# Q3: io/repository.py is the one module in io/ that opens files
# ---------------------------------------------------------------------------

#: Calls that reach the filesystem. A call by attribute is matched by name
#: alone, so `path.open()` and `tifffile.imread()` are both found.
OPENERS = frozenset(
    {
        "open",
        "read_bytes",
        "read_text",
        "write_bytes",
        "write_text",
        "imread",
        "imwrite",
        "memmap",
        "fromfile",
    }
)


def file_openers(source: str) -> list[str]:
    """Every call in `source` that opens a file, as `name:line`."""
    found = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call):
            func = node.func
            name = (
                func.id
                if isinstance(func, ast.Name)
                else func.attr
                if isinstance(func, ast.Attribute)
                else ""
            )
            if name in OPENERS:
                found.append(f"{name}:{node.lineno}")
    return found


def open_modes(source: str) -> list[str]:
    """The mode of every `open` call in `source` ("r" when none is given)."""
    modes = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call):
            func = node.func
            name = (
                func.id
                if isinstance(func, ast.Name)
                else func.attr
                if isinstance(func, ast.Attribute)
                else ""
            )
            if name == "open":
                positional = node.args[1:] if isinstance(func, ast.Name) else node.args
                mode = (
                    positional[0]
                    if positional
                    else next((k.value for k in node.keywords if k.arg == "mode"), None)
                )
                modes.append(
                    mode.value if isinstance(mode, ast.Constant) else "r" if mode is None else "?"
                )
    return modes


class TestQ3OneModuleOpensFiles:
    @pytest.mark.parametrize(
        "planted",
        [
            "open(p)",
            "with path.open('rb') as f: pass",
            "Path(p).read_bytes()",
            "tifffile.imread(p)",
            "np.memmap(p)",
        ],
    )
    def test_the_scanner_finds_a_planted_opener(self, planted: str) -> None:
        assert file_openers(planted), planted

    def test_the_scanner_ignores_streams_and_prose(self) -> None:
        clean = '"""Nothing here opens a file."""\ntif = tifffile.TiffFile(stream)\nstream.read()\n'
        assert file_openers(clean) == []

    def test_no_other_io_module_opens_a_file(self) -> None:
        offenders = {
            path.name: found
            for path in sorted(IO.glob("*.py"))
            if path.name != "repository.py"
            and (found := file_openers(path.read_text(encoding="utf-8")))
        }
        assert offenders == {}

    def test_the_repository_opens_files_read_only(self) -> None:
        source = (IO / "repository.py").read_text(encoding="utf-8")
        assert file_openers(source), "the repository must be the module that opens the tiles"
        modes = open_modes(source)
        assert modes and all(mode == "rb" for mode in modes), modes

    def test_the_package_docstring_says_so(self) -> None:
        """Q3 (a): `io/`'s rule narrowed, in the package's own docstring."""
        doc = ast.get_docstring(ast.parse((IO / "__init__.py").read_text(encoding="utf-8"))) or ""
        assert "repository.py" in doc
        assert "Nothing here opens a file." not in doc
