"""MP1-MP3 and DS1: each piece decodes, holds and reports only its own part of the DEM.

`docs/increments/23-basin-scale.md`, "23c-3: each piece's memory its own",
D1, D3, D5, D6 and "Red tests". A cut run plans the DEM from headers, then
loads one piece (`dem_input.load_piece`), refines it, writes its file and
drops it before the next. The spies are `piece_spies`'; none holds a piece.

- MP1 (projected, the four-tile scene): every decode, moved onto the
  lattice, lies inside the window of the piece being loaded, and covers it.
- MP2 (reprojected, DC10's scene): every source decode lies inside the run's
  source canvas and is strictly smaller than it; the run exits 0, so no
  check-point block read outside what the piece holds (`HeldWindows`).
- MP3: at each `load_piece` after the first, after `gc.collect()`, every
  earlier piece's tile, source canvas and written mesh arrays are dead; and
  when a piece's `write` returns, every check-point block yielded so far is.
  The finished pieces' files are on disk before the next piece loads (D6).
- DS1 (D6): a piece's overlap report is its own window's, read from the
  `seams` of each `PieceDem`. The report is a run-record input, which no
  piece file carries (`run_record.IN_FILE`), so the spy is where it shows.

The run's source canvas (MP2) is read from an uncut run of the same scene,
whose one decode is the canvas; nothing here re-derives the plan.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest

import piece_spies as ps
import pieces_fixtures as pf

# ------------------------------------------------------------------ MP1


def lattice_offset(tile: str, part: dict[str, Any]) -> tuple[int, int]:
    """Tile `tile`'s node (0, 0) on the partition's lattice, from the scene's
    layout (`quadrants` at `piece_spies.CUT`) and the lattice origin."""
    (ox, oy), h = part["origin"], part["spacing"]
    row, col = round((oy - pf.TIE_Y) / h), round((pf.TIE_X - ox) / h)
    north, west = tile.startswith("n"), tile[1] == "w"
    return row + (0 if north else ps.CUT), col + (0 if west else ps.CUT)


def test_mp1_each_piece_decodes_only_its_window(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = ps.spy_load_piece(monkeypatch)
    decodes = ps.spy_decodes(monkeypatch, calls)
    done = pf.cut(tmp_path / "run", *ps.four_tile_scene(tmp_path), "--pieces", "4")
    entries, part = done.index["pieces"], done.index["partition"]
    assert len(calls.loads) == len(entries) == 4, "one load_piece per piece, in index order"
    outside = [d for d in decodes if d.load is None or d.whole]
    assert outside == [], "decoded outside any piece's load, or a whole tile"
    for k, entry in enumerate(entries):
        w = entry["window"]
        rows = range(w["row0"], w["row0"] + w["rows"])
        cols = range(w["col0"], w["col0"] + w["cols"])
        covered: set[tuple[int, int]] = set()
        for d in (d for d in decodes if d.load == k):
            r, c = lattice_offset(d.name, part)
            r0, c0 = r + d.row0, c + d.col0
            box = (r0, r0 + d.rows - 1, c0, c0 + d.cols - 1)
            assert box[0] in rows and box[1] in rows, (entry["id"], d, w)
            assert box[2] in cols and box[3] in cols, (entry["id"], d, w)
            covered |= {
                (i, j) for i in range(box[0], box[1] + 1) for j in range(box[2], box[3] + 1)
            }
        missing = {(i, j) for i in rows for j in cols} - covered
        assert missing == set(), f"piece {entry['id']}: {len(missing)} window nodes not decoded"


# ------------------------------------------------------------------ MP2


@pytest.fixture(scope="module")
def run_canvas(tmp_path_factory: pytest.TempPathFactory) -> ps.Decode:
    """The uncut reprojected run's one decode: the run's source canvas
    (53 x 55 nodes of the 60 x 60 tile, measured on 23c-2)."""
    where = tmp_path_factory.mktemp("canvas")
    with pytest.MonkeyPatch.context() as mp:
        decodes = ps.spy_decodes(mp, None)
        done = pf.run(where / "run", *ps.reprojected_scene(where))
    assert done.code == 0, done.output
    assert len(decodes) == 1 and not decodes[0].whole, decodes
    assert (decodes[0].rows, decodes[0].cols) == (53, 55), decodes[0]
    return decodes[0]


def test_mp2_each_piece_decodes_a_part_of_the_run_source_canvas(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, run_canvas: ps.Decode
) -> None:
    calls = ps.spy_load_piece(monkeypatch)
    decodes = ps.spy_decodes(monkeypatch, calls)
    done = pf.cut(tmp_path / "run", *ps.reprojected_scene(tmp_path), "--pieces", "4")
    assert len(calls.loads) == len(done.index["pieces"]) == 4
    assert [d for d in decodes if d.load is None or d.whole] == []
    assert sorted({d.load for d in decodes}) == [0, 1, 2, 3], "a piece decoded nothing"
    c = run_canvas
    for d in decodes:
        assert c.row0 <= d.row0 and d.row0 + d.rows <= c.row0 + c.rows, (d, c)
        assert c.col0 <= d.col0 and d.col0 + d.cols <= c.col0 + c.cols, (d, c)
        assert d.rows * d.cols < c.rows * c.cols, (d, c)


# ------------------------------------------------------------------ MP3


def watched_run(
    where: Path, monkeypatch: pytest.MonkeyPatch, live: ps.Liveness, *args: str
) -> tuple[pf.Cut, list[set[str]]]:
    """A cut run with every piece's arrays watched by `live`, and the piece
    files on disk at each load."""
    out_dir = where / "run" / "x.vtk.pieces"
    on_disk: list[set[str]] = []

    def before(k: int) -> None:
        on_disk.append({p.name for p in out_dir.glob("*.vtk")} if out_dir.exists() else set())
        live.before_load(k)

    ps.spy_load_piece(monkeypatch, before, live.after_load)
    ps.spy_write(monkeypatch, live.after_write)
    return pf.cut(where / "run", *args), on_disk


def written_before_each_load(done: pf.Cut, on_disk: list[set[str]]) -> list[int]:
    """The loads at which some earlier piece's file was not yet on disk."""
    files = [e["file"] for e in done.index["pieces"]]
    return [k for k, seen in enumerate(on_disk) if seen != {f for f in files[:k] if f}]


def test_mp3_one_piece_alive_at_a_time_projected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The square scene at `--pieces 16`: tiles and written meshes."""
    live = ps.Liveness()
    done, on_disk = watched_run(
        tmp_path, monkeypatch, live, *ps.square_scene(tmp_path), "--pieces", "16"
    )
    assert len(live.groups) == len(done.index["pieces"]) == 16
    assert all(any(w.startswith("written") for w, _ in g) for g in live.groups), "a write unseen"
    assert live.problems == []
    assert written_before_each_load(done, on_disk) == []


def test_mp3_one_piece_alive_at_a_time_reprojected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """DC10's scene at `BLOCK` 8 (49 blocks, 42 meeting the domain): tiles,
    source canvases, written meshes, and every check-point block's arrays."""
    from tin_engine import target_grid

    monkeypatch.setattr(target_grid, "BLOCK", 8)
    live = ps.Liveness()
    ps.watch_source_canvases(monkeypatch, live)
    ps.watch_blocks(monkeypatch, live)
    done, on_disk = watched_run(
        tmp_path, monkeypatch, live, *ps.reprojected_scene(tmp_path), "--pieces", "4"
    )
    assert len(live.groups) == len(done.index["pieces"]) == 4
    assert all(any(w == "source canvas" for w, _ in g) for g in live.groups), "a canvas unseen"
    assert len(live.blocks) > 2 * 4, "too few blocks to judge: BLOCK not patched?"
    assert live.problems == []
    assert written_before_each_load(done, on_disk) == []


# ------------------------------------------------------------------ DS1

SEAM = re.compile(r"^(\S+) and (\S+) disagree at (\d+) nodes?, by up to ([\d.e+-]+) m")


def test_ds1_a_piece_reports_its_own_overlaps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`ne`'s copy of the column it shares with `nw` raised by 1 m: a piece
    whose window holds part of that column names `ne` and `nw` with the
    number of those nodes its window holds; a piece whose window holds none
    names neither. Heights are 30-90 m, so the 1 m step is far above
    `SEAM_THRESHOLD` (1 mm)."""
    reports: list[list[str]] = []
    ps.spy_load_piece(
        monkeypatch, after=lambda k, got: reports.append([s.entry() for s in got.seams])
    )
    done = pf.cut(tmp_path / "run", *ps.four_tile_scene(tmp_path, raised=1.0), "--pieces", "4")
    part, entries = done.index["partition"], done.index["pieces"]
    assert len(reports) == len(entries) == 4
    (ox, oy), h = part["origin"], part["spacing"]
    row, col = round((oy - pf.TIE_Y) / h), round((pf.TIE_X - ox) / h)
    raised = {(row + r, col + ps.CUT) for r in range(ps.CUT)}  # on the lattice
    run = (part["row0"], part["row0"] + part["rows"], part["col0"], part["col0"] + part["cols"])
    raised = {(r, c) for r, c in raised if run[0] <= r < run[1] and run[2] <= c < run[3]}
    holders = 0
    for entry, said in zip(entries, reports, strict=True):
        w = entry["window"]
        held = sum(
            w["row0"] <= r < w["row0"] + w["rows"] and w["col0"] <= c < w["col0"] + w["cols"]
            for r, c in raised
        )
        found = [m.groups() for m in map(SEAM.match, said) if m is not None]
        assert len(found) == len(said), f"an entry not in the report's form: {said}"
        if held:
            holders += 1
            assert [(a, b, int(n)) for a, b, n, _ in found] == [("ne.tif", "nw.tif", held)], (
                entry["id"],
                said,
            )
            assert float(found[0][3]) == pytest.approx(1.0, abs=1e-4), said
        else:
            assert said == [], (entry["id"], said)
    assert 0 < holders < len(entries), "the scene must have pieces on both sides"
