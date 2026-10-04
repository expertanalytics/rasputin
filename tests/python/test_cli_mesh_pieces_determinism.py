"""DC7 and DC8: a cut run's files depend on the input and the options alone.

`docs/increments/23-basin-scale.md`, K5 (determinism) and K6 (locality), and
"Tests @tester can write red", 23c, DC7 and DC8.

- DC7: refine's thread count does not reach a file. 23c has no `--threads`
  option (the design's `run_plan(plan, jobs, threads)` is driven by 23d's
  `--jobs`), so the count is forced from outside: every `refine` and
  `refine_points` the run calls, wherever it imported them from, is wrapped
  to pass `threads=1` or `threads=4`, and the wrapper counts its calls so
  the probe cannot pass by never running. Every file in the pieces
  directory is then compared byte for byte.
- DC8: "one piece re-run alone from its job spec writes the same bytes".
  23c pins no API for re-running one piece (`PieceJob`'s fields and the
  runner's signature are the developer's), so the same property is tested
  through what a piece's spec may depend on: the DEM changed deep inside
  another cell (five nodes or more from every partition line) leaves the
  first piece's file byte-identical, and changes the piece it lies in (the
  control). Confirmed under "Settled after 23c's red step (65e3990)", item 10.

The design's other half of DC7, "pieces in reversed order", moved to 23d's
PJ4 (item 12 there): 23c runs pieces in one fixed order and offers no hook.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import pieces_fixtures as pf
from tin_engine import _core


def files(directory: Path | None) -> dict[str, bytes]:
    assert directory is not None
    return {p.name: p.read_bytes() for p in sorted(directory.iterdir()) if p.is_file()}


def force_threads(monkeypatch: pytest.MonkeyPatch, threads: int) -> dict[str, int]:
    """Wrap `refine` and `refine_points` everywhere `tin_engine` holds them."""
    calls = {"refine": 0, "refine_points": 0}
    for name in calls:
        original: Callable[..., Any] = getattr(_core, name)

        def wrapped(*args: Any, _n: str = name, _f: Any = original, **kwargs: Any) -> Any:
            calls[_n] += 1
            kwargs["threads"] = threads
            return _f(*args, **kwargs)

        monkeypatch.setattr(_core, name, wrapped)
        for module in [m for k, m in list(sys.modules.items()) if k.startswith("tin_engine")]:
            if getattr(module, name, None) is original:
                monkeypatch.setattr(module, name, wrapped)
    return calls


@pytest.fixture(scope="module")
def scene(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, np.ndarray]:
    where = tmp_path_factory.mktemp("scene")
    array = pf.hills()
    return (
        pf.write_dem(where / "dem.tif", array),
        pf.write_domain(where / "sq.geojson", pf.SQUARE),
        array,
    )


def args(dem: Path, domain: Path, *extra: str) -> tuple[str, ...]:
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1", "--pieces", "4", *extra)


class TestDC7Threads:
    @pytest.mark.parametrize("extra", [(), ("--no-constraint-feet", "--start-min-angle", "0")])
    def test_one_and_four_threads_write_the_same_files(
        self,
        tmp_path: Path,
        scene: Any,
        monkeypatch: pytest.MonkeyPatch,
        extra: tuple[str, ...],
    ) -> None:
        dem, domain, _ = scene
        found = {}
        for threads in (1, 4):
            with monkeypatch.context() as patch:
                calls = force_threads(patch, threads)
                done = pf.cut(tmp_path / f"t{threads}", *args(dem, domain, *extra))
                assert calls["refine"] >= len(done.pieces), calls
            found[threads] = files(done.dir)
        assert found[1].keys() == found[4].keys()
        assert found[1] == found[4]

    def test_two_runs_write_the_same_files(self, tmp_path: Path, scene: Any) -> None:
        dem, domain, _ = scene
        first = pf.cut(tmp_path / "a", *args(dem, domain))
        second = pf.cut(tmp_path / "b", *args(dem, domain))
        assert files(first.dir) == files(second.dir)


class TestDC8Locality:
    def test_a_change_in_another_cell_leaves_a_piece_byte_identical(
        self, tmp_path: Path, scene: Any
    ) -> None:
        dem, domain, array = scene
        before = pf.cut(tmp_path / "before", *args(dem, domain))
        xs, ys = pf.seam_lines(before.pieces)
        (x,), (y,) = xs, ys
        col, row = int((x - pf.TIE_X) / pf.H), int((pf.TIE_Y - y) / pf.H)
        changed = array.copy()
        changed[row + 5 : row + 12, col + 5 : col + 12] += 9.0  # deep in cell (1, 1)
        other = pf.write_dem(tmp_path / "changed" / "dem.tif", changed)
        after = pf.cut(tmp_path / "after", *args(other, domain))

        def by_id(done: pf.Cut) -> dict[tuple[int, ...], bytes]:
            return {p.id: p.path.read_bytes() for p in done.pieces}

        old, new = by_id(before), by_id(after)
        assert old.keys() == new.keys() == {(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)}
        assert new[(0, 0, 0)] == old[(0, 0, 0)]
        assert new[(1, 1, 0)] != old[(1, 1, 0)], "the change reached no piece: a dead probe"
        sha = {tuple(e["id"]): e["sha256"] for e in before.index["pieces"]}
        sha_after = {tuple(e["id"]): e["sha256"] for e in after.index["pieces"]}
        assert sha[(0, 0, 0)] == sha_after[(0, 0, 0)]
