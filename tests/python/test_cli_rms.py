"""`rasputin mesh --rms` and the budget run's error summary: tests B3 to B6.

`docs/increments/104-triangle-budget.md`, *4.3 The error summary*, *4.4 The
record and the words* and *10. The red suite*. One pass over the written mesh,
on every budget run (the chosen trial) and on a tolerance run that asks with
`--rms`, gives `rms_error_m` and `error_nodes_measured` in the record; neither
goes into the mesh file.

The oracle (`source_summary`) recomputes both from the written `.vtk` and the
fixture's own nodes: the source DEM's valid nodes moved to the target CRS by
pyproj on the reprojected path, the DEM's nodes on the native one; each node
located by barycentric weights with a 1e-12 slack, its error against the
plane of the mesh's own vertices; a node on a vertex is not a point (as in
the scan). It never reads the record's figures or the scan's.

Invariant-critical (each node once, agreement with the record).
"""

from __future__ import annotations

import json
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
import tifffile
from numpy.testing import assert_allclose

from cli_helpers import COLS, ROWS, plain, rough_dem, runner, squashed
from geographic_fixtures import nodes, project
from geotiff_fixtures import DELTA_X, DELTA_Y, TIE_X, TIE_Y
from recordread import stats_names
from test_cli_mesh_stats import seconds
from tin_engine import _core
from tin_engine.cli import app
from vtkread import VtkFile, read_vtk

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
VELHAS = FIXTURES / "velhas"
TARGET = "EPSG:31983"
VELHAS_ARGS = (
    *("--dem", str(VELHAS / "anadem_velhas.tif")),
    *("--domain", str(VELHAS / "catchment.geojson")),
    *("--out-crs", TARGET),
)
RMS_KEYS = ("rms_error_m", "error_nodes_measured")
SUMMARY_ROW = "error summary"

#: Velhas at each tolerance, measured by `source_summary` on the meshes this
#: tree writes (2026-10-25): the source nodes inside the mesh that are not
#: vertices, and their root-mean-square error in metres.
VELHAS_FIGURES = {"1": (12_666, 0.42147), "2": (12_667, 0.77775), "0.5": (12_653, 0.21506)}

bumpy = rough_dem(14)


@dataclass(frozen=True)
class Run:
    """One successful `rasputin mesh` run: its mesh, record, `--stats` and stderr."""

    vtk: VtkFile
    record: dict[str, Any]
    stats: str
    stderr: str


def mesh_run(directory: Path, *args: str) -> Run:
    out, rec, md = (directory / f"x{suffix}" for suffix in (".vtk", ".json", ".md"))
    result = runner.invoke(
        app, ["mesh", *args, "--out", str(out), "--record", str(rec), "--stats", str(md)]
    )
    assert result.exit_code == 0, result.output
    record = json.loads(rec.read_text(encoding="ascii"))
    return Run(read_vtk(out.read_bytes()), record, md.read_text(encoding="utf-8"), result.stderr)


@dataclass(frozen=True)
class Measured:
    count: int
    rms: float
    largest: float


def source_summary(
    vtk: VtkFile, xy: npt.NDArray[np.float64], zp: npt.NDArray[np.float64]
) -> Measured:
    """Every node `xy` (in the mesh's CRS) with height `zp` that lies in a
    closed triangle of `vtk` and on none of its vertices: how many, their
    root-mean-square and largest |plane - z|, the plane from the file's vertices."""
    points = vtk.points
    tri = np.asarray(vtk.polygons, np.int64)
    origin = points[:, :2].min(axis=0)
    v = points[:, :2] - origin
    a, b, c = v[tri[:, 0]], v[tri[:, 1]], v[tri[:, 2]]
    za, zb, zc = points[tri[:, 0], 2], points[tri[:, 1], 2], points[tri[:, 2], 2]
    two_a = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    assert (np.abs(two_a) > 0).all(), "a zero-area output triangle"
    q = xy - origin
    error = np.full(len(q), np.nan)
    on_vertex = np.zeros(len(q), bool)
    for lo in range(0, len(q), 256):
        x, y = q[lo : lo + 256, 0, None], q[lo : lo + 256, 1, None]

        def weight(p: Any, r: Any, x: Any = x, y: Any = y) -> Any:
            return (
                (r[:, 0] - p[:, 0]) * (y - p[:, 1]) - (r[:, 1] - p[:, 1]) * (x - p[:, 0])
            ) / two_a

        wa, wb, wc = weight(b, c), weight(c, a), weight(a, b)
        inside = (wa >= -1e-12) & (wb >= -1e-12) & (wc >= -1e-12)
        plane = wa * za + wb * zb + wc * zc
        e = np.where(inside, np.abs(plane - zp[lo : lo + 256, None]), -1.0)
        error[lo : lo + 256] = np.where(inside.any(axis=1), e.max(axis=1), np.nan)
        top = np.maximum(np.maximum(wa, wb), wc)
        on_vertex[lo : lo + 256] = (inside & (top > 1.0 - 1e-9)).any(axis=1)
    measured = ~np.isnan(error) & ~on_vertex
    found = error[measured]
    return Measured(int(found.size), float(np.sqrt(np.mean(found**2))), float(found.max()))


def velhas_nodes() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """The fixture's valid nodes in TARGET: area-registered, as ANADEM is."""
    with tifffile.TiffFile(VELHAS / "anadem_velhas.tif") as tif:
        page = tif.pages.first
        (_, _, _, x0, y0, _), (dx, dy, _) = page.tags[33922].value, page.tags[33550].value
        nodata = float(page.tags[42113].value)
        array = page.asarray().astype(np.float64)
    lon, lat = nodes(x0 + 0.5 * dx, y0 - 0.5 * dy, dx, dy, array.shape)
    z = array.ravel()
    valid = np.isfinite(z) & (z != nodata)
    return project("EPSG:4674", TARGET, lon[valid], lat[valid]), z[valid]


def bumpy_nodes() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """`rough_dem(14)`'s nodes, point-registered, in its own CRS."""
    array = np.random.default_rng(14).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    x, y = nodes(TIE_X, TIE_Y, DELTA_X, DELTA_Y, array.shape)
    return np.column_stack([x, y]), array.ravel().astype(np.float64)


def agrees_with_the_oracle(run: Run, found: Measured) -> None:
    """The record's two figures are the oracle's, and its largest error the
    oracle's largest (the store keeps a node's position to 2 um)."""
    assert run.record["error_nodes_measured"] == found.count
    assert_allclose(run.record["rms_error_m"], found.rms, rtol=1e-6)
    assert_allclose(run.record["max_error_m"], found.largest, rtol=0, atol=1e-6)


class Summaries:
    """Every call of `_core.error_summary`, wherever a module bound it: the
    keyword arguments of each; calls through."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.kwargs: list[dict[str, Any]] = []
        real: Callable[..., Any] = _core.error_summary  # type: ignore[attr-defined]

        def spy(*args: Any, **kwargs: Any) -> Any:
            self.kwargs.append(kwargs)
            return real(*args, **kwargs)

        for name, module in list(sys.modules.items()):
            if name.startswith("tin_engine") and module is not None:
                for attr, value in list(vars(module).items()):
                    if value is real:
                        monkeypatch.setattr(module, attr, spy)

    @property
    def calls(self) -> int:
        return len(self.kwargs)

    def threads(self) -> list[Any]:
        return [kw.get("threads") for kw in self.kwargs]


@pytest.fixture
def summaries(monkeypatch: pytest.MonkeyPatch) -> Summaries:
    return Summaries(monkeypatch)


@pytest.fixture(scope="module")
def velhas_rms(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Run]:
    """Velhas `--tolerance T --rms` at 1, 2 and 0.5 m (four threads), and at
    1 m on one thread."""
    runs = {
        t: mesh_run(tmp_path_factory.mktemp(f"t{t}"), *VELHAS_ARGS, "--tolerance", t, "--rms",
                    "--threads", "4")
        for t in VELHAS_FIGURES
    }  # fmt: skip
    runs["1/1"] = mesh_run(
        tmp_path_factory.mktemp("t1-1"), *VELHAS_ARGS, "--tolerance", "1", "--rms", "--threads", "1"
    )
    return runs


@pytest.fixture(scope="module")
def velhas_oracle() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    return velhas_nodes()


# ---------------------------------------------------------------- B3: agreement with the record


class TestB3AgreesWithTheRecord:
    @pytest.mark.parametrize("tolerance", list(VELHAS_FIGURES))
    def test_the_figures_at_each_tolerance(
        self, velhas_rms: dict[str, Run], tolerance: str
    ) -> None:
        count, rms = VELHAS_FIGURES[tolerance]
        record = velhas_rms[tolerance].record
        assert record["error_nodes_measured"] == count
        assert record["rms_error_m"] == pytest.approx(rms, abs=1e-4)

    @pytest.mark.parametrize("tolerance", list(VELHAS_FIGURES))
    def test_the_oracle_over_the_written_mesh(
        self, velhas_rms: dict[str, Run], velhas_oracle: Any, tolerance: str
    ) -> None:
        run = velhas_rms[tolerance]
        agrees_with_the_oracle(run, source_summary(run.vtk, *velhas_oracle))

    def test_the_native_path(self, tmp_path: Path, bumpy: Path) -> None:
        """The strided start, the edge strip, the DEM's own nodes."""
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--stride", "8", "--tolerance", "2", "--rms")
        found = source_summary(run.vtk, *bumpy_nodes())
        assert found.count == ROWS * COLS - len(run.vtk.points)
        assert run.record["error_nodes_measured"] == found.count
        assert_allclose(run.record["rms_error_m"], found.rms, rtol=1e-9)

    def test_every_node_a_vertex_measures_none(self, tmp_path: Path, bumpy: Path) -> None:
        """Without a stride every node is a vertex: no node to measure, and the
        root-mean-square of none is 0."""
        run = mesh_run(tmp_path, "--dem", str(bumpy), "--tolerance", "2", "--rms")
        assert run.record["error_nodes_measured"] == 0
        assert run.record["rms_error_m"] == 0.0


# ---------------------------------------------------------------- B4: threads


class TestB4Threads:
    def test_one_thread_and_four_give_the_same_bits(self, velhas_rms: dict[str, Run]) -> None:
        one, four = velhas_rms["1/1"].record, velhas_rms["1"].record
        assert one["rms_error_m"] == four["rms_error_m"]
        assert one["error_nodes_measured"] == four["error_nodes_measured"]


# ---------------------------------------------------------------- B5: the words


class TestB5Words:
    def test_the_summary_line_ends_with_the_rms(self, velhas_rms: dict[str, Run]) -> None:
        stderr = squashed(plain(velhas_rms["1"].stderr))
        assert squashed("root-mean-square 0.42 m).") in stderr, stderr
        assert squashed("(largest difference 0.99789 m, root-mean-square") in stderr, stderr

    def test_stats_has_a_row_for_each_key_and_the_timing(self, velhas_rms: dict[str, Run]) -> None:
        stats = velhas_rms["1"].stats
        assert set(RMS_KEYS) <= set(stats_names(stats))
        assert SUMMARY_ROW in seconds(stats), sorted(seconds(stats))
        assert seconds(stats)[SUMMARY_ROW] >= 0.0

    def test_neither_key_is_in_the_mesh_file(self, velhas_rms: dict[str, Run]) -> None:
        fields = velhas_rms["1"].vtk.field_data
        for key in RMS_KEYS:
            assert key not in fields, key

    def test_without_rms_nothing_is_measured(self, tmp_path: Path, summaries: Summaries) -> None:
        run = mesh_run(tmp_path, *VELHAS_ARGS, "--tolerance", "1")
        assert summaries.calls == 0
        assert not set(RMS_KEYS) & set(run.record), sorted(run.record)
        assert "root-mean-square" not in plain(run.stderr)
        assert SUMMARY_ROW not in seconds(run.stats)
        assert not set(RMS_KEYS) & set(stats_names(run.stats))

    def test_with_rms_one_pass(self, tmp_path: Path, summaries: Summaries) -> None:
        mesh_run(tmp_path, *VELHAS_ARGS, "--tolerance", "1", "--rms", "--threads", "3")
        assert summaries.threads() == [3]

    def test_the_mcp_mesh_tool_lists_rms(self) -> None:
        pytest.importorskip("mcp")
        from tin_engine.mcp import schema as mcp_schema

        (tool,) = (d for d in mcp_schema.derive(app) if d.tool_name == "mesh")
        assert "rms" in tool.tool().input_schema["properties"]


# ---------------------------------------------------------------- B6: the budget run


class TestB6BudgetRun:
    def test_the_chosen_trial_is_measured_once_without_asking(
        self, tmp_path: Path, summaries: Summaries, velhas_oracle: Any
    ) -> None:
        """Five trials, one pass, on the mesh written."""
        run = mesh_run(tmp_path, *VELHAS_ARGS, "--triangles", "7386")
        assert summaries.calls == 1
        agrees_with_the_oracle(run, source_summary(run.vtk, *velhas_oracle))
        assert "root-mean-square" in plain(run.stderr)
        assert set(RMS_KEYS) <= set(stats_names(run.stats))
