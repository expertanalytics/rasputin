"""`edge_strip.py` alone, with a fake `_core`: increment 15f-3, PY2.

`docs/increments/15f-edge-strip.md`, D5 (`edge_strip.py`'s and
`final_check.py`'s rows), D6 and D7 (the `--stats` rows). The module is
`generate(view, start, clock) -> ConstraintCheckPoints` and `run(view, strip,
start, tolerance, clock) -> PointRefineOutcome`: the two calls and their clock
rows, no geometry. `final_check.run` gains `strip=None`, passed on to
`refine_points`.

The fake replaces the two `_core` functions wherever `edge_strip` reads them,
as module attributes or through a `_core` attribute; which import style the
module uses is not ruled, and the test does not pin it.

CHOSEN HERE, where D7 is silent: `generate` records exactly one row, `edge
strip: generate`; `run` records the outcome's own `scan_seconds` and
`split_seconds` as `edge strip: scan (parallel)` and `edge strip: split +
flip (serial)`, unchanged, as `final_check.run` does for its two rows.

Went red at `4157dab` because `tin_engine.edge_strip` did not exist (the
fixture's import failed) and `final_check.run` took no `strip`.

Not invariant-critical; no mutation round.
"""

from __future__ import annotations

import ast
import importlib
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest

import tin_engine.final_check as final_check
from tin_engine.stats import PhaseClock
from tin_engine.target_grid import TargetGrid


@pytest.fixture
def edge_strip() -> ModuleType:
    return importlib.import_module("tin_engine.edge_strip")


def ticking() -> PhaseClock:
    """A clock whose `now` advances 1 ms a reading, so a phase is never 0."""
    counter: Iterator[int] = iter(range(0, 10**12, 1_000_000))
    return PhaseClock(now=lambda: next(counter))


class FakeCore:
    """Records every call, in order, with its arguments."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = []
        self.strip = SimpleNamespace(size=7, no_data=1, duplicates=0, edge_count=4)
        self.outcome = SimpleNamespace(
            ok=lambda: True, scan_seconds=1.25, split_seconds=0.5, strip_inserted=3
        )

    def constraint_check_points(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append(("constraint_check_points", args, kwargs))
        return self.strip

    def refine_strip(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append(("refine_strip", args, kwargs))
        return self.outcome


@pytest.fixture
def fake(edge_strip: ModuleType, monkeypatch: pytest.MonkeyPatch) -> FakeCore:
    core = FakeCore()
    monkeypatch.setattr(edge_strip, "_core", core, raising=False)
    for name in ("constraint_check_points", "refine_strip"):
        monkeypatch.setattr(edge_strip, name, getattr(core, name), raising=False)
    return core


def start_mesh() -> SimpleNamespace:
    """What `refine` returns, as far as `edge_strip` reads it."""
    return SimpleNamespace(
        vertices=np.array([[0.0, 0.0], [30.0, 0.0], [0.0, -30.0]]),
        triangles=np.array([[0, 2, 1]], np.uint32),
        z=np.array([1.0, 2.0, 3.0]),
        valid=np.ones(3, bool),
        edges=np.array([[0, 1], [1, 2], [0, 2]], np.uint32),
        masks=np.array([0, 0, 4], np.uint32),
    )


def argument(call: tuple[str, tuple[Any, ...], dict[str, Any]], at: int, name: str) -> Any:
    _, args, kwargs = call
    return kwargs[name] if name in kwargs else args[at]


class TestGenerate:
    def test_calls_the_generator_with_the_starts_vertices_and_edges(
        self, edge_strip: ModuleType, fake: FakeCore
    ) -> None:
        view, start = object(), start_mesh()
        strip = edge_strip.generate(view, start, ticking())
        assert strip is fake.strip
        ((name, _, _),) = fake.calls
        assert name == "constraint_check_points"
        call = fake.calls[0]
        assert argument(call, 0, "view") is view
        np.testing.assert_array_equal(np.asarray(argument(call, 1, "vertices")), start.vertices)
        np.testing.assert_array_equal(np.asarray(argument(call, 2, "edges")), start.edges)

    def test_records_one_row(self, edge_strip: ModuleType, fake: FakeCore) -> None:
        clock = ticking()
        edge_strip.generate(object(), start_mesh(), clock)
        rows = clock.phases()
        assert [name for name, _ in rows] == ["edge strip: generate"]
        assert rows[0][1] > 0


class TestRun:
    def test_calls_refine_strip_with_the_strip_and_the_start(
        self, edge_strip: ModuleType, fake: FakeCore
    ) -> None:
        view, strip, start = object(), object(), start_mesh()
        out = edge_strip.run(view, strip, start, 0.75, ticking())
        assert out is fake.outcome
        ((name, _, kwargs),) = fake.calls
        assert name == "refine_strip"
        call = fake.calls[0]
        assert argument(call, 0, "view") is view
        assert argument(call, 1, "strip") is strip
        names = ("vertices", "triangles", "z", "valid", "edges", "masks")
        for at, field in enumerate(names, start=2):
            given = np.asarray(argument(call, at, field))
            np.testing.assert_array_equal(given, getattr(start, field), field)
        assert kwargs["tolerance"] == 0.75

    def test_records_the_outcomes_own_times(self, edge_strip: ModuleType, fake: FakeCore) -> None:
        clock = ticking()
        edge_strip.run(object(), object(), start_mesh(), 1.0, clock)
        rows = dict(clock.phases())
        assert rows["edge strip: scan (parallel)"] == 1.25
        assert rows["edge strip: split + flip (serial)"] == 0.5

    def test_generate_then_run_is_two_calls_in_order(
        self, edge_strip: ModuleType, fake: FakeCore
    ) -> None:
        clock, view, start = ticking(), object(), start_mesh()
        strip = edge_strip.generate(view, start, clock)
        edge_strip.run(view, strip, start, 1.0, clock)
        assert [name for name, _, _ in fake.calls] == ["constraint_check_points", "refine_strip"]
        assert argument(fake.calls[1], 1, "strip") is fake.strip
        assert [name for name, _ in clock.phases()][:1] == ["edge strip: generate"]


def test_the_module_holds_no_geometry_and_no_crs(edge_strip: ModuleType) -> None:
    """D1: Python never builds a point, and CRS never crosses into `_core`.
    The module imports neither shapely nor pyproj."""
    tree = ast.parse(Path(edge_strip.__file__).read_text(encoding="utf-8"))
    imported = {
        alias.name.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    } | {
        (node.module or "").split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
    }
    assert not imported & {"shapely", "pyproj"}, imported


class TestFinalCheckPassesTheStripOn:
    """D5, `final_check.py:28`: `run(..., strip=None)`, passed to `refine_points`."""

    @pytest.fixture
    def seen(self, monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
        calls: list[dict[str, Any]] = []

        def refine_points(*args: Any, **kwargs: Any) -> Any:
            calls.append({"args": args, **kwargs})
            return SimpleNamespace(scan_seconds=0.0, split_seconds=0.0)

        monkeypatch.setattr(final_check, "refine_points", refine_points)
        return calls

    @staticmethod
    def grid() -> TargetGrid:
        return TargetGrid(crs="EPSG:31983", spacing=30, row0=-10, col0=5, rows=4, cols=4)

    @staticmethod
    def strip_of(call: dict[str, Any]) -> Any:
        args = call["args"]
        return call["strip"] if "strip" in call else (args[7] if len(args) > 7 else None)

    def test_with_a_strip(self, seen: list[dict[str, Any]]) -> None:
        strip = object()
        final_check.run(start_mesh(), self.grid(), [], 1.0, ticking(), strip=strip)  # type: ignore[arg-type, call-arg]
        ((call,),) = (seen,)
        assert self.strip_of(call) is strip

    def test_without_one(self, seen: list[dict[str, Any]]) -> None:
        final_check.run(start_mesh(), self.grid(), [], 1.0, ticking())  # type: ignore[arg-type]
        ((call,),) = (seen,)
        assert self.strip_of(call) is None
