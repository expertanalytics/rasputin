"""`--z-unit` and what the user sees: increment 43, section 6, tests 15 and 16.

`docs/increments/43-units-and-projections.md`, sections 4.3 and 4.4. A 4 x 4
DEM in EPSG:2229 (US survey feet) whose heights are unstated is refused until
`--z-unit` says what they are; given, the mesh is in metres, heights and
coordinates. A metre DEM with silent heights is read as today, with one line
on stderr (Ola's ruling, section 8). A cached source whose catalogue entry
states metres says "metres" with no assumption.

HOW THIS FILE GOES RED. At the design's base there is no `--z-unit`, a feet
CRS is refused for its axes, nothing is printed for an assumed unit, and
`RemoteSource` has no `z_unit`.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from typer.testing import Result

from cli_driver import plain, refused, runner, write_tiff
from cog_fixtures import VARIANTS, build, write_cache
from geotiff_fixtures import (
    GT_MODEL_TYPE,
    GT_RASTER_TYPE,
    PIXEL_IS_AREA,
    PROJECTED_CS_TYPE,
    VERTICAL_UNITS,
    micro_tiff,
    with_keys,
)
from recordread import stats_row
from test_cli_mesh_mosaic import terrain
from tin_engine import sources as sources_module
from tin_engine.cli import app
from vtkread import read_vtk

FT_US = 1200 / 3937
METRE = 9001
#: Section 4.3: the line `mesh` prints when the DEM does not state its height unit.
ASSUMED_LINE = (
    "DEM heights read as metres; the file does not say. If they are feet, give "
    "--z-unit ft or us-ft."
)

#: A 4 x 4 grid of distinct heights, in US survey feet in the feet file.
HEIGHTS = (10 * np.arange(16).reshape(4, 4) + 100).astype(np.float32)
#: EPSG:2229, area-registered, 3 ftUS cells, no 4099: test 8's file at 4 x 4.
FEET_KEYS = {GT_MODEL_TYPE: 1, GT_RASTER_TYPE: PIXEL_IS_AREA, PROJECTED_CS_TYPE: 2229}
FEET_TIE = (0.0, 0.0, 0.0, 6_500_000.0, 2_100_000.0, 0.0)


def mesh(*args: str) -> Result:
    """`rasputin mesh`, stdout and stderr kept apart."""
    return runner.invoke(app, ["mesh", *args])


@pytest.fixture
def feet_dem(tmp_path: Path) -> Path:
    stream = micro_tiff(HEIGHTS, tiepoint=FEET_TIE, scale=(3.0, 3.0, 0.0), geokeys=FEET_KEYS)
    return write_tiff(tmp_path / "feet.tif", stream)


@pytest.fixture
def metre_dem(tmp_path: Path) -> Path:
    return write_tiff(tmp_path / "metres.tif", micro_tiff(HEIGHTS))


class TestTheFeetDem:
    def test_silent_heights_are_refused_naming_the_option(
        self, tmp_path: Path, feet_dem: Path
    ) -> None:
        refused(tmp_path, "mesh", "--dem", str(feet_dem), says=("--z-unit",))

    def test_given_us_ft_the_mesh_is_in_metres(self, tmp_path: Path, feet_dem: Path) -> None:
        out, md = tmp_path / "feet.vtk", tmp_path / "feet.md"
        result = mesh(
            "--dem", str(feet_dem), "--z-unit", "us-ft", "--out", str(out), "--stats", str(md)
        )
        assert result.exit_code == 0, plain(result.output)
        points = read_vtk(out.read_bytes()).points
        # Float64 products cast once to float32, as `to_metres` stores them.
        expected = (HEIGHTS.astype(np.float64) * FT_US).astype(np.float32).astype(np.float64)
        assert np.array_equal(np.sort(points[:, 2]), np.sort(expected.ravel()))
        # Node (0, 0) sits half a 3 ftUS cell inside the tie point, in metres.
        assert points[:, 0].min() == (6_500_000.0 + 1.5) * FT_US
        assert points[:, 1].max() == (2_100_000.0 - 1.5) * FT_US
        report = md.read_text(encoding="utf-8")
        assert stats_row(report, "dem_vertical_unit") == (
            "metres, converted from US survey feet (given with --z-unit)"
        )
        assert stats_row(report, "dem_file_crs") == (
            "EPSG:2229, in US survey feet; read as EPSG:26945, the same projection in metres"
        )
        assert ASSUMED_LINE not in result.stderr

    def test_a_unit_contradicting_the_file_is_refused(self, tmp_path: Path) -> None:
        stated = micro_tiff(HEIGHTS, geokeys=with_keys({VERTICAL_UNITS: METRE}))
        dem = write_tiff(tmp_path / "stated.tif", stated)
        refused(tmp_path, "mesh", "--dem", str(dem), "--z-unit", "ft", says=("--z-unit",))


class TestTheAssumptionLine:
    def test_a_silent_metre_dem_prints_it_once(self, tmp_path: Path, metre_dem: Path) -> None:
        result = mesh("--dem", str(metre_dem), "--out", str(tmp_path / "m.vtk"))
        assert result.exit_code == 0, plain(result.output)
        assert plain(result.stderr).count(ASSUMED_LINE) == 1, result.stderr

    def test_a_stated_metre_dem_does_not(self, tmp_path: Path) -> None:
        stated = micro_tiff(HEIGHTS, geokeys=with_keys({VERTICAL_UNITS: METRE}))
        dem = write_tiff(tmp_path / "stated.tif", stated)
        result = mesh("--dem", str(dem), "--out", str(tmp_path / "m.vtk"))
        assert result.exit_code == 0, plain(result.output)
        assert ASSUMED_LINE not in plain(result.stderr)


@pytest.mark.parametrize("command", ["mesh", "catchment", "station-catchments"])
def test_every_command_with_dem_lists_the_option(command: str) -> None:
    result = runner.invoke(app, [command, "--help"])
    assert result.exit_code == 0
    assert "--z-unit" in plain(result.output)


# ---------------------------------------------------------------------------
# Test 16: the catalogue states metres
# ---------------------------------------------------------------------------

KEY = "test-metres"


class TestTheCatalogue:
    def test_the_two_ruled_sources_state_metres(self) -> None:
        assert sources_module.SOURCES["glo30"].z_unit == "metre"
        assert sources_module.SOURCES["anadem-v1"].z_unit == "metre"

    def test_a_class_source_states_nothing(self) -> None:
        classes = [s for s in sources_module.SOURCES.values() if s.values == "classes"]
        assert classes, "no class source in the catalogue"
        assert all(s.z_unit is None for s in classes)

    @pytest.fixture
    def catalogue(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A projected test entry stating metres, put wherever a module bound `SOURCES`."""
        entry = sources_module.RemoteSource(
            id=KEY,
            kind="one-cog",
            url="https://example.invalid/test-metres.tif",
            crs="EPSG:25833",
            nodata=-9999.0,
            credit="Test elevation",
            licence_note="test fixture only",
            z_unit="metre",
        )
        original = sources_module.SOURCES
        patched = {**original, KEY: entry}
        for module in list(sys.modules.values()):
            name = getattr(module, "__name__", "") or ""
            if name.startswith("tin_engine") and getattr(module, "SOURCES", None) is original:
                monkeypatch.setattr(module, "SOURCES", patched)

    def test_a_cached_run_says_metres_and_prints_no_assumption(
        self, tmp_path: Path, catalogue: Any
    ) -> None:
        cache = tmp_path / "cache"
        write_cache(cache, KEY, {"dem": build(VARIANTS[0], array=terrain(50, 70))})
        md = tmp_path / "key.md"
        result = mesh(
            "--dem",
            KEY,
            "--cache",
            str(cache),
            "--out",
            str(tmp_path / "k.vtk"),
            "--stats",
            str(md),
        )
        assert result.exit_code == 0, plain(result.output)
        assert stats_row(md.read_text(encoding="utf-8"), "dem_vertical_unit") == "metres"
        assert ASSUMED_LINE not in plain(result.stderr)
