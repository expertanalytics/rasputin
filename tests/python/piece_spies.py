"""Spies for 23c-3's suites: what a cut run decodes, holds and writes, piece by piece.

`docs/increments/23-basin-scale.md`, "23c-3: each piece's memory its own",
"Red tests". Test support only; nothing in `src_python` uses it. Every name
23c-3 adds is looked up inside a function, so a missing one fails the test
that reaches it, not the collection of this helper.

- `patch_everywhere`: a function replaced on its module and in every
  `tin_engine` module that holds the same object (a name imported with
  `from ... import`), as `pieces_fixtures.node_spy` does for `node`.
- `Decodes`: every `load_window` and `load` a `TiffDemRepository` makes,
  with the index of the `load_piece` call it fell in (None outside any).
- `Liveness`: `weakref.finalize` on each piece's arrays (MP3); a group per
  `load_piece` call, judged at the next one after `gc.collect()`.
- `refine_calls`: the spy of `test_the_refusal_comes_before_any_piece_is_refined`
  (the seam pass, `refine`, `refine_points`, `basin_run.cut_run`).
- `window_bounds`: a piece window as inclusive `(r0, r1, c0, c1)`, whether it
  is an `IndexWindow` or 23c-2's `tile_window` tuple.
"""

from __future__ import annotations

import gc
import inspect
import sys
import weakref
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pytest


def patch_everywhere(
    monkeypatch: pytest.MonkeyPatch, owner: Any, name: str, replacement: Any
) -> Any:
    """`owner.name` replaced by `replacement` there and wherever `tin_engine`
    holds the same object; the original. An absent name raises
    `AttributeError`."""
    original = getattr(owner, name)
    monkeypatch.setattr(owner, name, replacement)
    for key, module in list(sys.modules.items()):
        held = key.startswith("tin_engine") and module is not owner
        if held and getattr(module, name, None) is original:
            monkeypatch.setattr(module, name, replacement)
    return original


def window_bounds(window: Any) -> tuple[int, int, int, int]:
    """Inclusive `(r0, r1, c0, c1)` of a window in the run window's indices."""
    if hasattr(window, "row0"):
        r0, c0 = int(window.row0), int(window.col0)
        return r0, r0 + int(window.rows) - 1, c0, c0 + int(window.cols) - 1
    r0, r1, c0, c1 = (int(v) for v in window)
    return r0, r1, c0, c1


# ------------------------------------------------------------------ the load and the write


@dataclass
class PieceCalls:
    """Each `load_piece` call's result in order, and the decodes inside each."""

    loads: list[Any] = field(default_factory=list)
    current: int | None = None


def spy_load_piece(
    monkeypatch: pytest.MonkeyPatch,
    before: Callable[[int], None] | None = None,
    after: Callable[[int, Any], None] | None = None,
) -> PieceCalls:
    """`dem_input.load_piece` wrapped: `before(k)` runs before the k-th call,
    `after(k, piece_dem)` after it returns. The spy keeps no `PieceDem`."""
    from tin_engine import dem_input

    calls = PieceCalls()
    original = None

    def spy(*args: Any, **kwargs: Any) -> Any:
        k = len(calls.loads)
        if before is not None:
            before(k)
        calls.current = k
        try:
            got = original(*args, **kwargs)  # type: ignore[misc]
        finally:
            calls.current = None
        calls.loads.append(None)  # the count only: holding it would keep it alive
        if after is not None:
            after(k, got)
        return got

    original = patch_everywhere(monkeypatch, dem_input, "load_piece", spy)
    return calls


def spy_write(monkeypatch: pytest.MonkeyPatch, after: Callable[[int, Any], None]) -> list[int]:
    """`basin_run.cut_run`'s `write` callable wrapped: `after(k, piece)` once
    the k-th piece's `write` returns. Returns the list of written indices."""
    from tin_engine import basin_run

    written: list[int] = []
    original = basin_run.cut_run
    signature = inspect.signature(original)

    def spy(*args: Any, **kwargs: Any) -> Any:
        bound = signature.bind(*args, **kwargs)
        assert "write" in bound.arguments, "cut_run takes no `write` callable"
        write = bound.arguments["write"]

        def wrapped(piece: Any) -> Any:
            entry = write(piece)
            k = len(written)
            written.append(k)
            after(k, piece)
            return entry

        bound.arguments["write"] = wrapped
        return original(*bound.args, **bound.kwargs)

    patch_everywhere(monkeypatch, basin_run, "cut_run", spy)
    return written


# ------------------------------------------------------------------ MP1, MP2: the decodes


@dataclass(frozen=True)
class Decode:
    name: str
    row0: int  # in the tile's own indices
    col0: int
    rows: int
    cols: int
    load: int | None  # the `load_piece` call it fell in
    whole: bool = False  # `load`, the whole tile


def spy_decodes(monkeypatch: pytest.MonkeyPatch, calls: PieceCalls | None) -> list[Decode]:
    """Every decode a `TiffDemRepository` makes, `load_window` and `load`."""
    from tin_engine.io.repository import TiffDemRepository

    seen: list[Decode] = []
    by_window = TiffDemRepository.load_window
    by_name = TiffDemRepository.load

    def load_window(self: Any, name: str, window: Any) -> Any:
        at = calls.current if calls is not None else None
        seen.append(Decode(name, window.row0, window.col0, window.rows, window.cols, at))
        return by_window(self, name, window)

    def load(self: Any, name: str) -> Any:
        tile = by_name(self, name)
        at = calls.current if calls is not None else None
        m = tile.meta
        seen.append(Decode(name, 0, 0, m.rows, m.cols, at, whole=True))
        return tile

    monkeypatch.setattr(TiffDemRepository, "load_window", load_window)
    monkeypatch.setattr(TiffDemRepository, "load", load)
    return seen


# ------------------------------------------------------------------ MP3: what is alive


def _owners(array: Any) -> list[Any]:
    """The array and every ndarray it is a view of: the memory is freed only
    when the last of them goes."""
    out = []
    while isinstance(array, np.ndarray):
        out.append(array)
        array = array.base
    return out


@dataclass
class Liveness:
    """Finalizers on each piece's arrays, grouped by the `load_piece` call."""

    groups: list[list[tuple[str, weakref.finalize]]] = field(default_factory=list)
    pending: list[tuple[str, weakref.finalize]] = field(default_factory=list)
    blocks: list[tuple[str, weakref.finalize]] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)

    def watch(self, into: list[tuple[str, weakref.finalize]], what: str, array: Any) -> None:
        for owner in _owners(array):
            into.append((what, weakref.finalize(owner, lambda: None)))

    def alive(self, groups: list[list[tuple[str, weakref.finalize]]]) -> list[str]:
        gc.collect()
        return [
            f"piece {g}: {what}" for g, group in enumerate(groups) for what, f in group if f.alive
        ]

    def before_load(self, k: int) -> None:
        if k:
            found = self.alive(self.groups[:k])
            if found:
                self.problems.append(f"at load {k}, still alive: {sorted(set(found))}")
        self.pending = []

    def after_load(self, k: int, got: Any) -> None:
        group = list(self.pending)
        self.watch(group, "tile", got.tile.array)
        self.groups.append(group)
        self.pending = []

    def after_write(self, k: int, piece: Any) -> None:
        group = self.groups[-1] if self.groups else self.pending
        made = piece.trimmed
        for name in ("vertices", "triangles", "edges", "edge_masks"):
            self.watch(group, f"written {name}", getattr(made, name))
        gc.collect()
        found = [what for what, f in self.blocks if f.alive]
        if found:
            self.problems.append(
                f"when piece {k}'s write returned, {len(found)} block arrays alive"
            )


def watch_source_canvases(monkeypatch: pytest.MonkeyPatch, live: Liveness) -> None:
    """Each `HeldWindows`' held tile, watched in the piece that builds it."""
    from tin_engine.target_grid import HeldWindows

    init = HeldWindows.__init__
    signature = inspect.signature(init)

    def spy(self: Any, *args: Any, **kwargs: Any) -> None:
        init(self, *args, **kwargs)
        held = signature.bind(self, *args, **kwargs).arguments["tile"]
        live.watch(live.pending, "source canvas", getattr(held, "array", held))

    monkeypatch.setattr(HeldWindows, "__init__", spy)


def watch_blocks(monkeypatch: pytest.MonkeyPatch, live: Liveness) -> None:
    """Every block `check_point_blocks` yields: its `xy` and `z` watched."""
    from tin_engine import target_grid

    original = None

    def spy(*args: Any, **kwargs: Any) -> Iterator[Any]:
        for block in original(*args, **kwargs):  # type: ignore[misc]
            live.watch(live.blocks, "xy", block[0])
            live.watch(live.blocks, "z", block[1])
            yield block
            del block

    original = patch_everywhere(monkeypatch, target_grid, "check_point_blocks", spy)


# ------------------------------------------------------------------ MP4: nothing refined


def refine_calls(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Spies on the seam pass, `refine`, `refine_points` and `basin_run.cut_run`,
    wherever `tin_engine` holds them; the names called, in order."""
    from tin_engine import _core, basin_run

    calls: list[str] = []
    for owner, name in (
        (_core, "refine"),
        (_core, "refine_points"),
        (_core, "refine_seam"),
        (basin_run, "cut_run"),
    ):
        original = getattr(owner, name)

        def spy(*a: object, _n: str = name, _f: object = original, **k: object) -> object:
            calls.append(_n)
            return _f(*a, **k)  # type: ignore[operator]

        patch_everywhere(monkeypatch, owner, name, spy)
    return calls


def piece_files(where: Path) -> dict[str, bytes]:
    """Every file of a pieces directory, `index.json` included, by name."""
    return {p.name: p.read_bytes() for p in sorted(where.iterdir())}


# ------------------------------------------------------------------ the three scenes


def square_scene(where: Path) -> tuple[str, ...]:
    """`pieces_fixtures.hills()` as one file, the SQUARE domain, `--tolerance 1`."""
    import pieces_fixtures as pf

    dem = pf.write_dem(where / "dem.tif", pf.hills())
    domain = pf.write_domain(where / "square.geojson", pf.SQUARE)
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1")


#: The four-tile scene's cut: tiles split after lattice row and column 23,
#: sharing that line (`mosaic_fixtures.quadrants(row_cut=23, col_cut=23, overlap=1)`).
CUT = 23


def four_tile_scene(where: Path, raised: float = 0.0) -> tuple[str, ...]:
    """The same DEM as four tiles. With `raised`, `ne`'s copy of the column it
    shares with `nw` alone (rows 0 to 22 of lattice column 23) is moved by it."""
    import pieces_fixtures as pf
    from mosaic_fixtures import meta, quadrants
    from test_cli_mesh_mosaic import write_tiles
    from tin_engine.io.models import DemTile

    whole = DemTile(meta=meta(rows=pf.N, cols=pf.N, dx=pf.H, dy=pf.H), array=pf.hills())
    tiles = quadrants(whole, row_cut=CUT, col_cut=CUT, overlap=1)
    if raised:
        ne = tiles["ne.tif"].array.copy()
        ne[:CUT, 0] += np.float32(raised)
        tiles["ne.tif"] = DemTile(meta=tiles["ne.tif"].meta, array=ne)
    write_tiles(where / "tiles", tiles)
    domain = pf.write_domain(where / "square.geojson", pf.SQUARE)
    return ("--dem", str(where / "tiles"), "--domain", str(domain), "--tolerance", "1")


def reprojected_scene(where: Path) -> tuple[str, ...]:
    """DC10's: a 60 x 60 float32 geographic tile, `--out-crs EPSG:31983`."""
    from cli_helpers import geojson
    from geographic_fixtures import geographic_tile_tiff, project_ring
    from test_cli_mesh_geographic import DOMAIN_RC, TARGET, lonlat_ring, tile_array

    where.mkdir(parents=True, exist_ok=True)
    dem = where / "anadem_like.tif"
    dem.write_bytes(geographic_tile_tiff(tile_array()).getvalue())
    ring = project_ring("EPSG:4326", "EPSG:4674", lonlat_ring(DOMAIN_RC))
    domain = geojson(where / "domain.geojson", ring, crs="EPSG:4674")
    return ("--dem", str(dem), "--domain", str(domain), "--out-crs", TARGET, "--tolerance", "1.0")
