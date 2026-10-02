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
from typing import Any

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


class TestS2DecodedDtype:
    """S2: a footprint carries the dtype its tile decodes to, from the header
    through `PROMOTION`, so the cap can count the canvas the plan will build."""

    @pytest.mark.parametrize(
        ("stored", "decoded_as"),
        [
            (np.int16, np.float32),
            (np.uint8, np.float32),
            (np.float32, np.float32),
            (np.int32, np.float64),
            (np.uint32, np.float64),
            (np.float64, np.float64),
        ],
        ids=["int16", "uint8", "float32", "int32", "uint32", "float64"],
    )
    def test_the_footprint_dtype_is_the_decoded_one(
        self, repo: ModuleType, tmp_path: Path, stored: Any, decoded_as: Any
    ) -> None:
        path = write(tmp_path / "a.tif", micro_tiff(elevations().astype(stored)))
        repository = repo.TiffDemRepository.from_directory(tmp_path)
        (footprint,) = repository.footprints()
        assert np.dtype(footprint.dtype) == decoded_as
        assert np.dtype(footprint.dtype) == decoded(path).array.dtype  # type: ignore[attr-defined]


ARCHIVE = Path(__file__).resolve().parents[3] / "rasputin_data" / "DTM10_UTM33_20220924"
needs_archive = pytest.mark.skipif(
    not ARCHIVE.is_dir(), reason=f"Ola's DTM10 archive is not at {ARCHIVE}"
)


@needs_archive
class TestB1RealArchive:
    """Ola's Q5 reading on the design's own 15a acceptance box, headers only.

    The 2 x 2 block 7908_3, 7908_2, 7808_4, 7808_1 is on the main lattice, and
    the half-cell tile 7807_2's 51-node overlap reaches into its south-west
    corner. Until the reading, the box was refused naming 7807_2. CI has no
    archive; the synthetic equivalent is `test_mosaic.py::TestB1LatticeByCoverage`
    and the real-extract one `test_dem_input.py::TestRealDtm10`.
    """

    BOX = (799750.0, 7849750.0, 900250.0, 7950250.0)
    BLOCK = ("7908_3", "7908_2", "7808_4", "7808_1")

    def test_the_acceptance_box_plans_on_the_main_lattice(self, repo: ModuleType) -> None:
        mz = importlib.import_module("tin_engine.mosaic")
        prints = repo.TiffDemRepository.from_directory(ARCHIVE).footprints()
        x_min, y_min, x_max, y_max = self.BOX
        box = mz.Bounds(x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)
        result = mz.plan_mosaic(prints, box, None)
        assert result.reference == (-100250.0, 7950250.0)  # the main lattice's (N2)
        assert (result.meta.rows, result.meta.cols) == (10051, 10051)
        assert (result.meta.x_min, result.meta.y_max) == (799750.0, 7950250.0)
        names = {p.name for p in result.tiles}
        assert {f"{b}_10m_z33.tif" for b in self.BLOCK} <= names
        assert "7807_2_10m_z33.tif" not in names
        # Dropping the other lattice's tiles is the same as their absence.
        odd = {  # N2: the half-cell tiles are 5 m east-west off
            f.name for f in prints if (f.meta.x_min - result.reference[0]) / f.meta.delta_x % 1 != 0
        }
        assert "7807_2_10m_z33.tif" in odd
        main = [f for f in prints if f.name not in odd]
        assert result == mz.plan_mosaic(main, box, None)


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


#: The one class in `repository.py` that writes (23a-2: the cache's write side).
WRITER = "CacheWriter"

#: Calls that change the filesystem, matched by name as `OPENERS` are.
WRITERS = frozenset(
    {"write_bytes", "write_text", "mkdir", "unlink", "rmdir", "rmtree", "touch", "rename",
     "replace", "symlink_to", "chmod"}
)  # fmt: skip


def split_out(source: str, name: str) -> tuple[str, str]:
    """`source` without the top-level class `name`, and that class alone."""
    tree = ast.parse(source)
    inside = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == name]
    rest = [n for n in tree.body if n not in inside]
    unparse = lambda nodes: ast.unparse(ast.Module(body=nodes, type_ignores=[]))  # noqa: E731
    return unparse(rest), unparse(inside)


def file_writers(source: str) -> list[str]:
    """Every call in `source` that changes the filesystem, as `name:line`.
    `str.replace` is told apart by its receiver: only `os.replace` counts."""
    found = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            name, receiver = node.func.attr, node.func.value
            if name == "replace" and not (isinstance(receiver, ast.Name) and receiver.id == "os"):
                continue
            if name in WRITERS:
                found.append(f"{name}:{node.lineno}")
    return found


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
        """Outside `CacheWriter` (23a-2, the cache's write side), every open is
        "rb" and nothing writes; `CacheWriter` itself is where the writes are."""
        source = (IO / "repository.py").read_text(encoding="utf-8")
        assert file_openers(source), "the repository must be the module that opens the tiles"
        reading, writer = split_out(source, WRITER)
        modes = open_modes(reading)
        assert modes and all(mode == "rb" for mode in modes), modes
        assert file_writers(reading) == []
        assert file_writers(writer), f"{WRITER} must be where the cache is written"

    @pytest.mark.parametrize(
        "planted",
        [
            "def f(p):\n    return open(p, 'wb')\n",
            "def f(p):\n    return os.open(p, os.O_CREAT)\n",
            "def f(p):\n    p.write_bytes(b'')\n",
            "class Other:\n    def f(self, p):\n        os.replace(p, p)\n",
            "def f(p):\n    shutil.rmtree(p)\n",
        ],
    )
    def test_the_read_only_check_finds_a_planted_write(self, planted: str) -> None:
        source = f"class {WRITER}:\n    def put(self, p):\n        p.write_bytes(b'')\n{planted}"
        reading, writer = split_out(source, WRITER)
        assert file_writers(writer)
        modes = open_modes(reading)
        assert file_writers(reading) or any(mode != "rb" for mode in modes), planted

    def test_the_package_docstring_says_so(self) -> None:
        """Q3 (a): `io/`'s rule narrowed, in the package's own docstring."""
        doc = ast.get_docstring(ast.parse((IO / "__init__.py").read_text(encoding="utf-8"))) or ""
        assert "repository.py" in doc
        assert "Nothing here opens a file." not in doc
