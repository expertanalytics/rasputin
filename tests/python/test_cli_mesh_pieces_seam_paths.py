"""A cut run's seams end to end, on paths the synthetic scenes never take.

`docs/increments/23-basin-scale.md`, "The seam protocol" (steps 2, 3 and 5,
K4) and "Output" (the index's piece windows). The scenes' seam edges all
reach the cut in the seam pass's own end order, the seam pass's z equals
refine's at every seam vertex, both pieces always agree, and the run's
window starts within a node of the DEM's corner; so the tests here bend
one of those, through a wrapper around the seam pass or the cut, or a
domain set well inside a larger DEM, and check that the files still say
what the design says.

Heights: the hills are 30-90 m (`pieces_fixtures.hills`), so a shift of
`DELTA` (1 mm) is some 10^11 ulps there, never lost to rounding.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

import pieces_fixtures as pf
from tin_engine import basin_run

DELTA = 1e-3  # m; see the module docstring for the scale


@pytest.fixture(scope="module")
def scene(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    where = tmp_path_factory.mktemp("seam_paths")
    return (
        pf.write_dem(where / "dem.tif", pf.hills()),
        pf.write_domain(where / "square.geojson", pf.SQUARE),
    )


def cut_args(scene: tuple[Path, Path], *extra: str) -> tuple[str, ...]:
    dem, domain = scene
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1", "--pieces", "4", *extra)


def shifted(outcome: Any, by: float) -> SimpleNamespace:
    """The seam pass's answer with every height moved by `by`."""

    def move(z: float | None) -> float | None:
        return None if z is None else z + by

    return SimpleNamespace(
        a=tuple(outcome.a),
        points=np.asarray(outcome.points),
        z=np.asarray(outcome.z) + by,
        z_a=move(outcome.z_a),
        z_b=move(outcome.z_b),
    )


def xyz_set(done: pf.Cut) -> set[frozenset[pf.XYZ]]:
    u = pf.union(done.pieces)
    return {frozenset(pf.xyz(u.points[i]) for i in t) for t in u.triangles.tolist()}


def test_the_direction_a_seam_edge_is_stored_in_does_not_reach_a_file(
    tmp_path: Path, scene: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every constraint edge handed to the cut reversed: the slice turns the
    seam pass's points round and the seam strip's box is the same, so the
    union is the same set of triangles, z included."""
    today = pf.cut(tmp_path / "as-noded", *cut_args(scene))
    original = basin_run.cut_run
    called = []

    def reversed_edges(load: Any, start: Any, edges: np.ndarray, *rest: Any) -> Any:
        called.append(1)
        return original(load, start, np.ascontiguousarray(edges[:, ::-1]), *rest)

    # `cli.py` imports the cut path inside the function that runs it (D8 of
    # 23c-3), so the wrapper goes on `basin_run` itself.
    monkeypatch.setattr(basin_run, "cut_run", reversed_edges)
    turned = pf.cut(tmp_path / "reversed", *cut_args(scene))
    assert called == [1], "the wrapper did not reach the run"
    assert [p.id for p in turned.pieces] == [p.id for p in today.pieces]
    assert xyz_set(turned) == xyz_set(today)
    assert pf.conformity_problems(turned.pieces) == []


def test_seam_vertices_carry_the_seam_pass_heights(
    tmp_path: Path, scene: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A seam vertex's z in a piece file is the seam pass's, written over
    refine's own: with the seam pass's heights moved by DELTA, every file
    vertex at a seam-pass point has the moved height."""
    original = basin_run.refine_seam
    told: dict[tuple[float, float], float] = {}

    def moved(strip: Any, a: Any, b: Any, **kwargs: Any) -> SimpleNamespace:
        out = shifted(original(strip, a, b, **kwargs), DELTA)
        ends = [(out.a, out.z_a), ((b if tuple(a) == out.a else a), out.z_b)]
        inner = zip(map(tuple, out.points.tolist()), out.z.tolist(), strict=True)
        for xy, z in [*ends, *inner]:
            if z is not None:
                told[(float(xy[0]), float(xy[1]))] = z
        return out

    monkeypatch.setattr(basin_run, "refine_seam", moved)
    done = pf.cut(tmp_path / "run", *cut_args(scene))
    found = [
        (xy, float(p[2]), told[xy])
        for piece in done.pieces
        for p in piece.points
        if (xy := (float(p[0]), float(p[1]))) in told
    ]
    assert len({xy for xy, _, _ in found}) == len(told), "a seam-pass point missing from the files"
    assert [f for f in found if f[1] != f[2]] == []
    assert pf.conformity_problems(done.pieces) == []


def test_a_seam_two_pieces_write_differently_fails_the_run(
    tmp_path: Path, scene: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """K4: the second piece to compute a seam gets heights moved by DELTA,
    so the two pieces' seam records differ; the run exits 1 naming the
    seam, writes no index, and keeps the piece files written (23c-3, D6)."""
    original = basin_run.refine_seam
    calls: Counter[tuple[tuple[float, ...], tuple[float, ...]]] = Counter()
    inner = 0

    def second_differs(strip: Any, a: Any, b: Any, **kwargs: Any) -> Any:
        nonlocal inner
        out = original(strip, a, b, **kwargs)
        calls[(tuple(a), tuple(b))] += 1
        inner += len(np.asarray(out.z))
        return out if calls[(tuple(a), tuple(b))] == 1 else shifted(out, DELTA)

    monkeypatch.setattr(basin_run, "refine_seam", second_differs)
    done = pf.run(tmp_path / "run", *cut_args(scene))
    assert set(calls.values()) == {2}, "each seam edge computed once by each of its two pieces"
    assert inner > 0, "no seam-pass point: only the ends could differ"
    assert done.code == 1, done.output
    assert "Error: seam edge " in done.output
    assert "the pieces do not conform" in done.output
    # 23c-3, D6: each piece is written when done, and K4 is judged after the
    # last, so the piece files stay; no index says the run is unfinished.
    assert done.dir is not None and not (done.dir / "index.json").exists()
    assert sorted(p.name for p in done.dir.iterdir()) == [
        "0-0-0.vtk", "0-1-0.vtk", "1-0-0.vtk", "1-1-0.vtk"
    ]  # fmt: skip
    assert re.search(r"\bno index", " ".join(done.output.split())), done.output


def test_each_piece_window_in_the_index_holds_its_cell_and_its_vertices(tmp_path: Path) -> None:
    """The index places each piece's window on the lattice of the partition
    record's `origin`, `spacing`, `row0` and `col0`. The domain sits some
    twenty nodes inside a larger DEM, so a window counted from the run's
    window instead of the lattice's node (0, 0) is off by that much."""
    inset = 200.0  # m: 20 nodes in from the DEM's corner
    dem = pf.write_dem(tmp_path / "dem.tif", pf.hills(n=81))
    square = [(x + inset, y - inset) for x, y in pf.SQUARE]
    domain = pf.write_domain(tmp_path / "square.geojson", square)
    done = pf.cut(tmp_path / "run", *cut_args((dem, domain)))
    part = done.index["partition"]
    (ox, oy), h, dx, dy = part["origin"], part["spacing"], part["dx"], part["dy"]
    assert part["col0"] >= 15 and part["row0"] >= 15, part
    by_id = {p.id: p for p in done.pieces}
    assert len(by_id) == 4
    for entry in done.index["pieces"]:
        j, i, _ = entry["id"]
        w = entry["window"]
        x0, x1 = ox + w["col0"] * h, ox + (w["col0"] + w["cols"] - 1) * h
        y1, y0 = oy - w["row0"] * h, oy - (w["row0"] + w["rows"] - 1) * h
        # The cell, cut to the run's window, lies inside the piece's window.
        cols, rows = part["cols"], part["rows"]
        c0, c1 = part["col0"] + i * dx, part["col0"] + min((i + 1) * dx, cols - 1)
        r0, r1 = part["row0"] + j * dy, part["row0"] + min((j + 1) * dy, rows - 1)
        assert w["col0"] <= c0 and c1 <= w["col0"] + w["cols"] - 1, (entry, c0, c1)
        assert w["row0"] <= r0 and r1 <= w["row0"] + w["rows"] - 1, (entry, r0, r1)
        points = by_id[tuple(entry["id"])].points
        assert (x0 <= points[:, 0]).all() and (points[:, 0] <= x1).all(), (entry, x0, x1)
        assert (y0 <= points[:, 1]).all() and (points[:, 1] <= y1).all(), (entry, y0, y1)
