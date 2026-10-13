"""`rasputin mesh|catchment|station-catchments --threads N --memory-budget B`:
increment 38, T6 to T8 and T10's command-line cases.

`docs/increments/38-threads-and-memory.md`, sections 4.2 to 4.4 and 8. The
count a run resolves (`--threads`, else `RASPUTIN_THREADS`, else every core)
must reach every `_core` call that takes `threads` and every Python pool, and
be printed as "Threads used: N." Each site has its own spy, because a dropped
keyword falls back to the C++ default (every core) or the library default
(one worker) and the mesh comes out the same: only the spy sees it.

Spies wrap the real function in the module that calls it and call through,
so every run still writes its mesh, which each test asserts. Pools are
observed as `test_io_cog.py` observes one: the pool class monkeypatched in
the module under test, recording `max_workers`.

Nothing here times anything: the machine's load is not a test's business.

T10's figures are the fixture's (`test_cli_catchment.dem_dir`, the bowl in
four tiles, through the lake's box plus the margin): the first window is 49 x
49 nodes, 9 604 bytes at float32 for the mosaic cap, 14 406 to flood. They
were measured on this tree (`5d2e706d`) with the cap faked to twice each
budget, through the CLI as well as `delineate`, and match the design's
table (review round 2).
"""

from __future__ import annotations

import importlib
import json
import re
from pathlib import Path
from typing import Any

import pytest

import tin_engine.cli as cli
import tin_engine.edge_strip as edge_strip
import tin_engine.feature_input as feature_input
import tin_engine.final_check as final_check
import tin_engine.io.cog as cog
import tin_engine.land_cover_fractions as land_cover_fractions
import tin_engine.target_grid as target_grid
from cli_driver import USAGE, UTM33, invoke, plain, polygon_file, ran, refused, rough_dem, runner
from cog_fixtures import write_cache
from feature_fixtures import Feat, moved, write_geojson
from fraction_fixtures import write_class_cache, write_dem, write_domain
from recordread import stats_row
from test_cli_catchment import dem_dir, lakes, utm_seed
from test_cli_mesh_cache import KEY as DEM_KEY
from test_cli_mesh_cache import catalogue, data, sources, tif
from test_cli_mesh_features import GALLERY, SQUARE
from test_cli_mesh_geographic import TARGET, VELHAS
from test_cli_mesh_land_cover_fractions import classes, covered, offline
from test_cli_mesh_stats import section
from test_cli_tolerance_slope import STEEP, box, v2
from tin_engine.cli import app
from tin_engine.crs import same_crs
from vtkread import read_vtk

__all__ = [  # fixtures, used by name
    "box", "catalogue", "classes", "data", "dem_dir", "lakes", "offline", "sources", "tif", "v2",
]  # fmt: skip

bumpy = rough_dem(16)
plain_square = polygon_file(SQUARE)

#: The three commands that take the two options.
TAKERS = ("mesh", "catchment", "station-catchments")
VELHAS_RUN = (
    *("--dem", str(VELHAS / "anadem_velhas.tif")),
    *("--domain", str(VELHAS / "catchment.geojson")),
    *("--out-crs", TARGET, "--tolerance", "5"),
)


@pytest.fixture(autouse=True)
def no_threads_variable(monkeypatch: pytest.MonkeyPatch) -> None:
    """The developer's shell is not the test's default."""
    monkeypatch.delenv("RASPUTIN_THREADS", raising=False)


@pytest.fixture
def dem_cache(tmp_path: Path, data: bytes) -> Path:
    """`test_cli_mesh_cache`'s cached DEM, under its own name here."""
    write_cache(tmp_path / "dem-cache", DEM_KEY, {"dem": data})
    return tmp_path / "dem-cache"


@pytest.fixture
def class_scene(tmp_path: Path) -> tuple[Path, Path, Path]:
    """The land-cover fractions scene: (dem, domain, cache)."""
    return (
        write_dem(tmp_path / "dem.tif"),
        write_domain(tmp_path / "domain.geojson"),
        write_class_cache(tmp_path / "class-cache"),
    )


class Calls:
    """Every call of one function as a module uses it: the positional and
    keyword arguments, and what it returned; calls through."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch, module: Any, name: str) -> None:
        self.args: list[tuple[Any, ...]] = []
        self.kwargs: list[dict[str, Any]] = []
        self.out: list[Any] = []
        real = getattr(module, name)

        def spy(*args: Any, **kwargs: Any) -> Any:
            self.args.append(args)
            self.kwargs.append(kwargs)
            self.out.append(real(*args, **kwargs))
            return self.out[-1]

        monkeypatch.setattr(module, name, spy)

    def threads(self) -> list[Any]:
        """The `threads` keyword of each call, None where it was not given."""
        return [kw.get("threads") for kw in self.kwargs]


class Pools:
    """The worker count of every `ThreadPoolExecutor` one module builds."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch, module: Any) -> None:
        self.workers: list[int | None] = []
        real = module.ThreadPoolExecutor

        def spy(*args: Any, **kwargs: Any) -> Any:
            self.workers.append(kwargs.get("max_workers", args[0] if args else None))
            return real(*args, **kwargs)

        monkeypatch.setattr(module, "ThreadPoolExecutor", spy)


def mesh(*args: str, env: dict[str, str] | None = None) -> Any:
    """`rasputin mesh <args>`, required to succeed; the runner's result."""
    result = runner.invoke(app, ["mesh", *args], env=env)
    assert result.exit_code == 0, result.output
    return result


def written(out: Path) -> None:
    """The run wrote a mesh with triangles: the spies called through."""
    assert out.is_file(), out
    assert len(read_vtk(out.read_bytes()).polygons) > 0


def machine_cores() -> int:
    return int(importlib.import_module("tin_engine.resources").machine_cores())


# ---------------------------------------------------------------- T6


class TestT6TheOptionAndTheVariable:
    def test_the_option_sets_the_stats_row(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "x.vtk"
        args = ("--dem", str(bumpy), "--tolerance", "1", "--threads", "3", "--out", str(out))
        result = mesh(*args, "--stats", "-")
        assert "Threads used: 3." in result.stdout
        assert "hardware concurrency" not in result.stdout

    def test_the_variable_is_the_default(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "x.vtk"
        args = ("--dem", str(bumpy), "--tolerance", "1", "--out", str(out), "--stats", "-")
        result = mesh(*args, env={"RASPUTIN_THREADS": "5"})
        assert "Threads used: 5." in result.stdout

    def test_the_option_beats_the_variable(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "x.vtk"
        args = ("--dem", str(bumpy), "--tolerance", "1", "--threads", "3", "--out", str(out))
        result = mesh(*args, "--stats", "-", env={"RASPUTIN_THREADS": "5"})
        assert "Threads used: 3." in result.stdout

    @pytest.mark.parametrize("how", ["--threads 0", "RASPUTIN_THREADS=0", "nothing"])
    def test_zero_or_nothing_is_every_core(self, tmp_path: Path, bumpy: Path, how: str) -> None:
        out = tmp_path / "x.vtk"
        args = ["--dem", str(bumpy), "--tolerance", "1", "--out", str(out), "--stats", "-"]
        env = {"RASPUTIN_THREADS": "0"} if how == "RASPUTIN_THREADS=0" else None
        if how == "--threads 0":
            args += ["--threads", "0"]
        result = mesh(*args, env=env)
        assert f"Threads used: {machine_cores()}." in result.stdout

    def test_a_stride_run_prints_the_row_too(self, tmp_path: Path, bumpy: Path) -> None:
        out = tmp_path / "x.vtk"
        result = mesh("--dem", str(bumpy), "--threads", "3", "--out", str(out), "--stats", "-")
        assert "Threads used: 3." in result.stdout

    def test_the_record_and_the_inputs_table_carry_the_count(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        out, md, rec = tmp_path / "x.vtk", tmp_path / "x.md", tmp_path / "x.json"
        mesh(*("--dem", str(bumpy), "--tolerance", "1", "--threads", "3", "--out", str(out)),
             "--stats", str(md), "--record", str(rec))  # fmt: skip
        record = json.loads(rec.read_text(encoding="ascii"))
        assert record["threads"] == 3 and type(record["threads"]) is int
        report = md.read_text(encoding="utf-8")
        assert stats_row(report, "threads") == "3"
        inputs = section(report, "Inputs")
        assert re.search(r"^\| Threads used \| 3 \| `?threads`? \|$", inputs, re.M), inputs
        vtk = read_vtk(out.read_bytes())
        assert "threads" not in vtk.field_data, "how a run was scheduled is not a mesh field"

    def test_a_negative_count_is_refused_naming_the_option(
        self, tmp_path: Path, bumpy: Path
    ) -> None:
        args = ("--dem", str(bumpy), "--tolerance", "1", "--threads", "-1")
        refused(tmp_path, "mesh", *args, says=("--threads",))

    @pytest.mark.parametrize("text", ["abc", "-1"])
    def test_a_variable_that_is_not_a_count_is_refused_naming_it(
        self, tmp_path: Path, bumpy: Path, text: str
    ) -> None:
        out = tmp_path / "x.vtk"
        args = ["mesh", "--dem", str(bumpy), "--tolerance", "1", "--out", str(out)]
        result = runner.invoke(app, args, env={"RASPUTIN_THREADS": text})
        output = plain(result.output)
        assert result.exit_code == USAGE, output
        assert "Traceback" not in output
        assert "RASPUTIN_THREADS" in output and text in output, output
        assert "0 means every core" in output
        assert not out.exists()

    @pytest.mark.parametrize("command", TAKERS)
    def test_the_help_lists_both_options(self, command: str) -> None:
        text = ran(command, "--help")
        assert "--threads" in text and "--memory-budget" in text
        assert "RASPUTIN_THREADS" in text

    def test_catchment_runs_with_a_count(self, tmp_path: Path, dem_dir: Path, lakes: Path) -> None:
        out = tmp_path / "c.geojson"
        args = ("--dem", str(dem_dir), *utm_seed(), "--lakes", str(lakes), "--threads", "3")
        ran("catchment", *args, "--out", str(out))
        assert out.is_file()


# ---------------------------------------------------------------- T7


class TestT7EveryCoreCallGetsTheCount:
    """Five `_core` entry points take `threads`; each call site passes the count."""

    def test_refine_and_refine_strip(
        self, tmp_path: Path, bumpy: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Calls(monkeypatch, cli, "refine")
        strip = Calls(monkeypatch, edge_strip, "refine_strip")
        out = tmp_path / "x.vtk"
        mesh("--dem", str(bumpy), "--tolerance", "1", "--threads", "3", "--out", str(out))
        assert phase.threads() == [3]
        assert strip.threads() == [3]
        written(out)

    def test_refine_points_on_the_resampled_path(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        points = Calls(monkeypatch, final_check, "refine_points")
        out = tmp_path / "x.vtk"
        mesh(*VELHAS_RUN, "--threads", "3", "--out", str(out))
        assert points.threads() == [3]
        written(out)

    def test_slope_tolerance(
        self, tmp_path: Path, v2: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        slope = Calls(monkeypatch, cli, "SlopeTolerance")
        out = tmp_path / "x.vtk"
        args = ("--dem", str(v2), "--domain", str(box), "--tolerance", "10", *STEEP)
        mesh(*args, "--threads", "3", "--out", str(out))
        assert slope.threads() == [3]
        written(out)

    def test_class_areas(
        self,
        classes: Any,
        offline: None,
        tmp_path: Path,
        class_scene: tuple[Path, Path, Path],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        areas = Calls(monkeypatch, land_cover_fractions._core, "class_areas")
        run = covered(tmp_path, *class_scene, "--threads", "3")
        assert areas.threads() and set(areas.threads()) == {3}
        assert len(run.vtk.polygons) > 0


# ---------------------------------------------------------------- T8


class TestT8EveryPoolGetsTheCount:
    """Three pools, and the three reprojector sites whose arrays reach one."""

    def test_the_block_decoder_on_a_cached_source(
        self, catalogue: Any, dem_cache: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        pools = Pools(monkeypatch, cog)
        out = tmp_path / "x.vtk"
        args = ("--dem", DEM_KEY, "--cache", str(dem_cache), "--tolerance", "1", "--threads", "3")
        mesh(*args, "--out", str(out))
        assert pools.workers and set(pools.workers) == {3}, pools.workers
        written(out)

    def test_the_block_decoder_on_a_tile(
        self, tif: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        pools = Pools(monkeypatch, cog)
        out = tmp_path / "x.vtk"
        mesh("--dem", str(tif), "--tolerance", "1", "--threads", "3", "--out", str(out))
        assert pools.workers and set(pools.workers) == {3}, pools.workers
        written(out)

    def test_the_target_grids_pool_and_the_serial_reprojector_inside_it(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """`resample` and `check_point_blocks` each open a pool of the count;
        the reprojector their workers call takes `threads=1` (no pool in a
        pool). The spy also sees the rim's call (`source_region`), which
        passes nothing and so is serial by the library default."""
        pools = Pools(monkeypatch, target_grid)
        moves = Calls(monkeypatch, target_grid, "reprojector")
        out = tmp_path / "x.vtk"
        mesh(*VELHAS_RUN, "--threads", "3", "--out", str(out))
        assert len(pools.workers) >= 2 and set(pools.workers) == {3}, pools.workers
        explicit = [kw["threads"] for kw in moves.kwargs if "threads" in kw]
        assert len(explicit) >= 2 and set(explicit) == {1}, moves.kwargs
        assert all(kw.get("threads", 1) == 1 for kw in moves.kwargs)
        written(out)

    def test_the_feature_reprojector(
        self, tmp_path: Path, bumpy: Path, plain_square: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """`feature_input.py`'s `_Tally._take` builds a reprojector from the
        file's CRS to the DEM's when they differ (the gallery in EPSG:4326
        over the EPSG:25833 DEM); that call carries the count. The hull's
        call runs the other way (DEM to file) and is not the one asserted."""
        moves = Calls(monkeypatch, feature_input, "reprojector")
        lonlat = [Feat(f.fid, moved(f.geometry, UTM33, "EPSG:4326"), f.properties) for f in GALLERY]  # type: ignore[arg-type]
        path = write_geojson(tmp_path / "wgs.geojson", lonlat, crs=None)
        out = tmp_path / "x.vtk"
        args = ("--dem", str(bumpy), "--domain", str(plain_square), "--tolerance", "1")
        mesh(*args, "--features", str(path), "--threads", "3", "--out", str(out))
        from_the_file = [
            kw
            for args_, kw in zip(moves.args, moves.kwargs, strict=True)
            if same_crs(args_[0], "EPSG:4326") and same_crs(args_[1], UTM33)
        ]
        assert from_the_file, moves.args
        assert all(kw.get("threads") == 3 for kw in from_the_file), from_the_file
        written(out)

    def test_the_class_fraction_reprojector(
        self,
        classes: Any,
        offline: None,
        tmp_path: Path,
        class_scene: tuple[Path, Path, Path],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """`land_cover_fractions.index_space` moves every mesh vertex to the
        class raster's longitude and latitude: the one site in the module."""
        moves = Calls(monkeypatch, land_cover_fractions, "reprojector")
        run = covered(tmp_path, *class_scene, "--threads", "3")
        assert moves.threads() and set(moves.threads()) == {3}, moves.kwargs
        assert len(run.vtk.polygons) > 0


# ---------------------------------------------------------------- T10, the command line


class TestT10TheMemoryBudgetOnTheCommandLine:
    def test_8k_is_refused_by_the_mosaic_cap_naming_the_figures(
        self, tmp_path: Path, dem_dir: Path, lakes: Path
    ) -> None:
        out = tmp_path / "c.geojson"
        args = ("--dem", str(dem_dir), *utm_seed(), "--lakes", str(lakes), "--memory-budget", "8K")
        code, output = invoke("catchment", *args, "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output and "Traceback" not in output
        for word in ("49 x 49", "9604", "8192", "--memory-budget"):
            assert word in output, f"{word!r} not in {output!r}"
        assert not out.exists()

    def test_a_wrong_suffix_is_refused_naming_the_text(
        self, tmp_path: Path, dem_dir: Path, lakes: Path
    ) -> None:
        out = tmp_path / "c.geojson"
        args = ("--dem", str(dem_dir), *utm_seed(), "--lakes", str(lakes), "--memory-budget", "8KB")
        code, output = invoke("catchment", *args, "--out", str(out))
        assert code == USAGE, output
        assert "No such option" not in output and "Traceback" not in output
        assert "8KB" in output and "--memory-budget" in output, output
        assert not out.exists()

    def test_mesh_takes_the_budget_too(self, tmp_path: Path, bumpy: Path) -> None:
        """Ola's ruling 1: `mesh`'s mosaic cap is the budget. `rough_dem` is
        17 x 21 nodes, 1 428 bytes at float32; one byte refuses, 8G passes."""
        args = ("--dem", str(bumpy), "--tolerance", "1", "--memory-budget", "1")
        refused(tmp_path, "mesh", *args, says=("17 x 21", "1428", "--memory-budget"))
        out = tmp_path / "x.vtk"
        mesh("--dem", str(bumpy), "--tolerance", "1", "--memory-budget", "8G", "--out", str(out))
        written(out)
