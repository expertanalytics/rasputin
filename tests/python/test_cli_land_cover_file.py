"""`rasputin mesh --land-cover FILE.tif [--land-cover-classes NAME|FILE.toml]`.

`docs/increments/44-catalogue-as-data.md`, section 4.5 and section 7's T11
(the CLI by path) and T8's CLI halves (a class raster in US survey feet; one
whose GeoKeys say its heights are in US survey feet). The scene is
`test_cli_mesh_land_cover_fractions.py`'s: its DEM in SIRGAS 2000 / UTM 23S
over `class_fixtures`' raster (EPSG:4326, cells of 1/128 degree), the class
source installed in `SOURCES` under `class_fixtures.KEY` with two years in a
cache, the network closed.

THE ORACLES: the cached route itself, array for array (the same 2025 bytes,
read by path, give the same table), and that suite's
`fraction_fixtures.oracle_fractions`, which recomputes every triangle's
fractions from the file's own vertices with pyproj, GEOS and `Geod`
(bound 1e-6 on a fraction, as there: float32 rounding is 6e-8). "Before any
phase is timed" is read off a spy on `PhaseClock.add`, which every timed
phase goes through.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- Each refusal is a usage error (exit 2), and none times a phase.
- `--land-cover-year` with a file says "no years"; a missing path names the
  path and lists the catalogue keys; `--land-cover-classes` alone says "needs
  --land-cover"; a float32 raster names the file and `uint8`.
- A file run needs neither `--cache` nor `RASPUTIN_DATA`.
- A feet raster over another place than the mesh gets past the option check
  and is refused by the "outside" check `assign` makes, not by the
  height-unit check.
"""

from __future__ import annotations

import socket
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import feature_fixtures as ff
from class_fixtures import (
    AREA_KEYS,
    KEY,
    LAT0,
    LON0,
    PATH_TEMPLATE,
    STEP,
    TILE,
    class_grid,
    class_tiff,
    per_year_entry,
)
from cli_helpers import USAGE, invoke, squashed
from fraction_fixtures import (
    oracle_fractions,
    rows_of,
    same_rows,
    write_class_cache,
    write_dem,
    write_domain,
)
from geotiff_fixtures import micro_tiff
from recordread import file_field
from test_class_raster_crs import US_SURVEY_FOOT, feet_classes
from test_cli_fetch import install
from test_cli_mesh_cache import no_network
from tin_engine import sources
from tin_engine.stats import PhaseClock
from vtkread import VtkFile, read_vtk

TOLERANCE = "3"
ARRAYS = ("land_cover_offsets", "land_cover_classes", "land_cover_fractions")
CREDIT_FIELDS = ("land_cover_credit", "land_cover_licence_note", "land_cover_cite")
FRACTION = 1e-6
FILE_NAME = "lc2025.tif"


@pytest.fixture
def classes(monkeypatch: pytest.MonkeyPatch) -> Any:
    template = "https://example.invalid/" + PATH_TEMPLATE
    return install(sources, monkeypatch, per_year_entry(sources, template))


@pytest.fixture
def offline(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)


@pytest.fixture
def scene(tmp_path: Path) -> dict[str, Path]:
    lc = tmp_path / FILE_NAME
    lc.write_bytes(class_tiff(2025))
    return {
        "cache": write_class_cache(tmp_path / "cache"),
        "dem": write_dem(tmp_path / "dem.tif"),
        "domain": write_domain(tmp_path / "domain.geojson"),
        "file": lc,
    }


@pytest.fixture
def phases(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """The names of every phase timed in the run."""
    seen: list[str] = []
    original = PhaseClock.add

    def spy(self: PhaseClock, name: str, seconds: float) -> None:
        seen.append(name)
        original(self, name, seconds)

    monkeypatch.setattr(PhaseClock, "add", spy)
    return seen


def mesh(scene: dict[str, Path], target: Path, *extra: str) -> tuple[int, str]:
    return invoke(
        "mesh",
        *("--dem", str(scene["dem"]), "--domain", str(scene["domain"]), "--tolerance", TOLERANCE),
        *("--out", str(target)),
        *extra,
    )


def table(vtk: VtkFile) -> tuple[np.ndarray, ...]:
    for name in ARRAYS:
        assert name in vtk.field_data, f"no {name}; the file has {sorted(vtk.field_data)}"
    return tuple(np.asarray(vtk.field_data[n].values) for n in ARRAYS)


# ------------------------------------------------------------ the file route


class TestByPath:
    def test_a_file_gives_the_cached_routes_table_and_names_itself(
        self, classes: Any, offline: None, tmp_path: Path, scene: dict[str, Path]
    ) -> None:
        cached, by_path = tmp_path / "cached.vtk", tmp_path / "path.vtk"
        code, output = mesh(scene, cached, "--cache", str(scene["cache"]), "--land-cover", KEY)
        assert code == 0, output
        code, output = mesh(scene, by_path, "--land-cover", str(scene["file"]))
        assert code == 0, output
        want, got = read_vtk(cached.read_bytes()), read_vtk(by_path.read_bytes())
        for ours, theirs in zip(table(got), table(want), strict=True):
            assert_array_equal(ours, theirs)
        assert_array_equal(
            got.cell_array("land_cover_code").values, want.cell_array("land_cover_code").values
        )
        offsets, held, fractions = table(got)
        triangles = np.array([np.asarray(p) for p in got.polygons], dtype=np.int64)
        rows = rows_of(offsets, held, fractions.astype(np.float32))
        assert (
            same_rows(rows, oracle_fractions(got.points, triangles, class_grid(2025)), FRACTION)
            == []
        )
        assert file_field(got, "land_cover_source") == FILE_NAME
        assert not set(CREDIT_FIELDS) & set(got.field_data), sorted(got.field_data)

    def test_the_record_strings_name_the_file(
        self, offline: None, tmp_path: Path, scene: dict[str, Path]
    ) -> None:
        target = tmp_path / "x.vtk"
        code, output = mesh(scene, target, "--land-cover", str(scene["file"]))
        assert code == 0, output
        vtk = read_vtk(target.read_bytes())
        assert file_field(vtk, "land_cover_source") == FILE_NAME
        assert file_field(vtk, "land_cover_codes") == (
            f"{FILE_NAME} class number; 0 = no data, and every constraint line"
        )

    def test_a_land_cover_map_names_and_colours_the_gpkg(
        self, offline: None, tmp_path: Path, scene: dict[str, Path]
    ) -> None:
        target = tmp_path / "named.gpkg"
        code, output = mesh(
            scene,
            target,
            "--land-cover",
            str(scene["file"]),
            "--land-cover-classes",
            "mapbiomas-c11",
        )
        assert code == 0, output
        palette = ff.classes("mapbiomas-c11").legend.colours
        with closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT code, land_cover, colour FROM triangles").fetchall()
        assert len({r[0] for r in rows}) > 2
        for code_, name, colour in rows:
            assert (name, colour) == palette[code_], code_


# ------------------------------------------------------------ refusals


def float32_tiff() -> bytes:
    return micro_tiff(
        np.zeros((40, 48), dtype=np.float32), tiepoint=(0, 0, 0, LON0, LAT0, 0),
        scale=(STEP, STEP, 0), geokeys=AREA_KEYS, tile=TILE,
    ).getvalue()  # fmt: skip


def says(output: str, *words: str) -> None:
    """Each of `words` in `output`, whitespace ignored on both sides: Rich
    breaks a long path anywhere inside its panel, at a column that depends on
    the terminal's width and on the length of the temporary directory."""
    for word in words:
        assert squashed(word) in squashed(output), (word, output)


class TestRefusals:
    """Each a usage error naming the option, before any phase is timed."""

    def refused(
        self, scene: dict[str, Path], tmp_path: Path, phases: list[str], *extra: str
    ) -> str:
        target = tmp_path / "x.vtk"
        code, output = mesh(scene, target, *extra)
        assert code == USAGE, output
        assert "Nosuchoption" not in squashed(output), output
        assert "Traceback" not in squashed(output), output
        assert not target.exists()
        assert phases == [], phases
        return output

    def test_a_year_with_a_file(
        self, offline: None, tmp_path: Path, scene: dict[str, Path], phases: list[str]
    ) -> None:
        output = self.refused(
            scene, tmp_path, phases,
            "--land-cover", str(scene["file"]), "--land-cover-year", "2025",
        )  # fmt: skip
        says(output, "--land-cover-year", "no years")

    def test_a_missing_path_lists_the_catalogue(
        self, offline: None, tmp_path: Path, scene: dict[str, Path], phases: list[str]
    ) -> None:
        missing = tmp_path / "nope.tif"
        output = self.refused(scene, tmp_path, phases, "--land-cover", str(missing))
        says(output, "--land-cover", "nope.tif", "anadem-v1", "glo30", "mapbiomas-c11")

    def test_a_float32_raster(
        self, offline: None, tmp_path: Path, scene: dict[str, Path], phases: list[str]
    ) -> None:
        heights = tmp_path / "heights.tif"
        heights.write_bytes(float32_tiff())
        output = self.refused(scene, tmp_path, phases, "--land-cover", str(heights))
        says(output, "--land-cover", "heights.tif", "uint8")

    def test_a_map_without_codes(
        self, offline: None, tmp_path: Path, scene: dict[str, Path], phases: list[str]
    ) -> None:
        output = self.refused(
            scene, tmp_path, phases,
            "--land-cover", str(scene["file"]), "--land-cover-classes", "property",
        )  # fmt: skip
        says(output, "--land-cover-classes", "property")

    def test_a_map_without_land_cover(
        self, offline: None, tmp_path: Path, scene: dict[str, Path], phases: list[str]
    ) -> None:
        output = self.refused(scene, tmp_path, phases, "--land-cover-classes", "mapbiomas-c11")
        says(output, "--land-cover-classes", "needs --land-cover")

    def test_a_file_stating_feet_heights_is_refused_naming_it(
        self, offline: None, tmp_path: Path, scene: dict[str, Path], phases: list[str]
    ) -> None:
        """4.7: a class raster whose GeoKeys say its heights are in US survey
        feet, read with `z_unit="metre"`, is refused."""
        feet = tmp_path / "feet_heights.tif"
        feet.write_bytes(feet_classes(vertical=US_SURVEY_FOOT))
        output = self.refused(scene, tmp_path, phases, "--land-cover", str(feet))
        says(output, "feet_heights.tif", "US survey feet")


def test_a_feet_raster_passes_the_option_check(
    offline: None, tmp_path: Path, scene: dict[str, Path]
) -> None:
    """T8: no height-unit refusal for a class raster in a feet CRS. The raster
    lies in California and the mesh in Brazil, so the run reaches the
    fractions and is refused there, as any mesh past its raster is."""
    feet = tmp_path / "feet_classes.tif"
    feet.write_bytes(feet_classes())
    target = tmp_path / "x.vtk"
    code, output = mesh(scene, target, "--land-cover", str(feet))
    assert code == USAGE, output
    says(output, "outside")
    assert "--z-unit" not in squashed(output) and "neither" not in squashed(output), output
