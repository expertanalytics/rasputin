"""The four red tests "Settled after 23c's red step (65e3990)" asks for.

`docs/increments/23-basin-scale.md`, that block's items 2, 8 and 14 and its
question B:

- 2: a cut run writes nothing at `--out` and says on stderr where the pieces
  are; piece files are `<j>-<i>-<k>.<ext>` in `<out>.pieces/`; a rerun over
  an existing directory deletes its `index.json` first and writes the new one
  last, so a rerun that fails (a truncated DEM) leaves no `index.json`.
- 8: `--pieces` or `--memory-budget` without `--domain` is a usage error
  naming `--domain`, and no file is written.
- 14: a partition line through a water ring given by `--features`: the ring
  does not divide pieces, its crossings with the seam are vertices of both
  pieces (off-node, the seam vertices whose z is bilinear), and conformity
  and both section 3D oracles hold over the union.
- B: a cell all NoData keeps its index entry with `file` and `sha256` null
  and zero `triangles` and `vertices`; the other pieces are written and pass
  the oracles. A cut of an all-NoData DEM is refused (exit 2) with no index.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon

import pieces_fixtures as pf
from feature_fixtures import Feat, write_geojson
from geotiff_fixtures import micro_tiff
from test_cli_mesh_dem import USAGE

SENTINEL = -32767.0


def args(dem: Path, domain: Path | None, *extra: str) -> tuple[str, ...]:
    where = ("--domain", str(domain)) if domain is not None else ()
    return ("--dem", str(dem), *where, "--tolerance", "1", *extra)


def scene(tmp_path: Path) -> tuple[Path, Path, np.ndarray]:
    array = pf.hills()
    return (
        pf.write_dem(tmp_path / "dem" / "dem.tif", array),
        pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE),
        array,
    )


def lines_of(done: pf.Cut) -> tuple[float, float]:
    xs, ys = pf.seam_lines(done.pieces)
    assert len(xs) == 1 and len(ys) == 1, (xs, ys)
    return xs.pop(), ys.pop()


def not_for_want_of_the_option(output: str) -> None:
    assert "No such option" not in output, "refused for the wrong reason"
    assert "Traceback" not in output


# ------------------------------------------------------------------ item 2


class TestTheOutputDirectory:
    def test_nothing_at_out_and_stderr_names_the_pieces(self, tmp_path: Path) -> None:
        dem, domain, _ = scene(tmp_path)
        done = pf.cut(tmp_path / "run", *args(dem, domain, "--pieces", "4"))
        assert not done.out.exists(), "23c writes nothing at --out in a cut run"
        assert done.dir is not None
        assert done.dir.name == done.out.name + ".pieces"
        assert done.dir.name in "".join(done.output.split())

    def test_piece_files_are_named_by_their_ids(self, tmp_path: Path) -> None:
        dem, domain, _ = scene(tmp_path)
        done = pf.cut(tmp_path / "run", *args(dem, domain, "--pieces", "4"))
        for entry in done.index["pieces"]:
            j, i, k = entry["id"]
            assert entry["file"] == f"{j}-{i}-{k}.vtk", entry

    def test_ply_pieces_are_named_by_their_ids_in_their_own_directory(self, tmp_path: Path) -> None:
        dem, domain, _ = scene(tmp_path)
        done = pf.run(tmp_path / "run", *args(dem, domain, "--pieces", "4"), out_name="x.ply")
        assert done.code == 0, done.output
        assert done.dir is not None and done.dir.name == "x.ply.pieces"
        index = json.loads((done.dir / "index.json").read_text())
        for entry in index["pieces"]:
            j, i, k = entry["id"]
            assert entry["file"] == f"{j}-{i}-{k}.ply", entry
            assert (done.dir / entry["file"]).exists()

    def test_a_rerun_that_fails_leaves_no_index(self, tmp_path: Path) -> None:
        dem, domain, array = scene(tmp_path)
        where = tmp_path / "run"
        first = pf.cut(where, *args(dem, domain, "--pieces", "4"))
        assert first.dir is not None and (first.dir / "index.json").exists()

        data = micro_tiff(array, scale=(pf.H, pf.H, 0.0)).getvalue()
        cut = tmp_path / "cut" / "dem.tif"
        cut.parent.mkdir()
        cut.write_bytes(data[: len(data) - array.nbytes // 2])
        again = pf.run(where, *args(cut, domain, "--pieces", "4"))
        assert again.code == USAGE, again.output
        not_for_want_of_the_option(again.output)
        assert not (first.dir / "index.json").exists(), "the old index outlived a failed rerun"
        assert not any(p.name.startswith("index") for p in first.dir.iterdir())


# ------------------------------------------------------------------ item 8


@pytest.mark.parametrize("option", [("--pieces", "4"), ("--memory-budget", "4096")])
def test_cutting_needs_a_domain(tmp_path: Path, option: tuple[str, str]) -> None:
    dem, _, _ = scene(tmp_path)
    done = pf.run(tmp_path / "run", *args(dem, None, *option))
    assert done.code == USAGE, done.output
    not_for_want_of_the_option(done.output)
    assert "--domain" in done.output
    assert not done.out.exists()
    assert done.dir is None
    assert [p.name for p in done.out.parent.iterdir()] == []


# ------------------------------------------------------------------ item 14


def test_a_partition_line_through_a_water_ring(tmp_path: Path) -> None:
    dem, domain, array = scene(tmp_path)
    x, y = lines_of(pf.cut(tmp_path / "first", *args(dem, domain, "--pieces", "4")))
    lake = Polygon(
        [
            (x - 40.3, y + 90.7),
            (x - 38.9, y + 30.1),
            (x + 41.1, y + 29.3),
            (x + 39.7, y + 91.9),
        ]
    )
    seam = LineString([(x, pf.TIE_Y), (x, pf.TIE_Y - 400.0)])
    crossings = shapely.get_coordinates(lake.exterior.intersection(seam))
    assert len(crossings) == 2
    assert not any(float((cy - pf.TIE_Y) / pf.H).is_integer() for _, cy in crossings), "on a node"
    lakes = write_geojson(tmp_path / "lake.geojson", [Feat("l", lake, {"property": "water"})])
    done = pf.cut(tmp_path / "run", *args(dem, domain, "--pieces", "4", "--features", str(lakes)))
    assert lines_of(done) == (x, y), "the lake moved the partition"
    assert sorted(p.id for p in done.pieces) == [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)]
    left, right = (next(p for p in done.pieces if p.id == (0, i, 0)) for i in (0, 1))
    for cx, cy in crossings:
        found = []
        for piece in (left, right):
            near = np.flatnonzero(
                np.hypot(piece.points[:, 0] - cx, piece.points[:, 1] - cy) <= pf.SNAP
            )
            assert len(near) == 1, (piece.id, (cx, cy), near)
            found.append(pf.xyz(piece.points[near[0]]))
            water = pf.bit_named(piece.vtk, "water")
            assert any(
                m & water and near[0] in (a, b)
                for (a, b), m in zip(piece.edges, piece.masks, strict=True)
            ), (piece.id, "the ring has no edge at the crossing")
        assert found[0] == found[1], "the crossing differs between the two pieces"
        assert found[0][0] == x
    assert pf.problems(done, array, Polygon(pf.SQUARE), 1.0) == []


# ------------------------------------------------------------------ question B


def read_index(done: pf.Cut) -> dict[str, Any]:
    assert done.dir is not None, done.output
    return json.loads((done.dir / "index.json").read_text())


def test_a_cell_all_nodata_keeps_its_entry_and_writes_no_file(tmp_path: Path) -> None:
    dem, domain, array = scene(tmp_path)
    x, y = lines_of(pf.cut(tmp_path / "first", *args(dem, domain, "--pieces", "4")))
    col, row = int((x - pf.TIE_X) / pf.H), int((pf.TIE_Y - y) / pf.H)
    holed = array.copy()
    holed[: row + 1, col:] = SENTINEL  # the top right cell, its seams included
    path = tmp_path / "holed" / "dem.tif"
    path.parent.mkdir()
    path.write_bytes(micro_tiff(holed, scale=(pf.H, pf.H, 0.0), nodata="-32767").getvalue())

    done = pf.run(tmp_path / "run", *args(path, domain, "--pieces", "4"))
    assert done.code == 0, done.output
    index = read_index(done)
    entries = {tuple(e["id"]): e for e in index["pieces"]}
    assert set(entries) == {(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)}
    empty = entries[(0, 1, 0)]
    assert empty["file"] is None and empty["sha256"] is None, empty
    assert (empty["counts"]["triangles"], empty["counts"]["vertices"]) == (0, 0), empty
    assert done.dir is not None
    assert not any(p.name.startswith("0-1-0") for p in done.dir.iterdir())
    pieces = []
    for key, entry in entries.items():
        if key == (0, 1, 0):
            continue
        file = done.dir / entry["file"]
        assert entry["sha256"] == pf.sha256(file)
        pieces.append(pf.read_piece(key, file))  # type: ignore[arg-type]
        assert entry["counts"]["triangles"] == len(pieces[-1].triangles) > 0
    assert pf.conformity_problems(pieces) == []
    u = pf.union(pieces)
    assert pf.overlapping_pairs(u) == []
    assert pf.delaunay_violations(u.points[:, :2], u.triangles, u.constrained) == []
    valid = holed != SENTINEL
    located, over = pf.located_over(u, holed, valid, 1.0)
    assert located > 500 and over == 0


def test_a_cut_of_an_all_nodata_dem_is_refused_with_no_index(tmp_path: Path) -> None:
    _, domain, array = scene(tmp_path)
    void = np.full_like(array, SENTINEL)
    path = tmp_path / "void" / "dem.tif"
    path.parent.mkdir()
    path.write_bytes(micro_tiff(void, scale=(pf.H, pf.H, 0.0), nodata="-32767").getvalue())
    done = pf.run(tmp_path / "run", *args(path, domain, "--pieces", "4"))
    assert done.code == USAGE, done.output
    not_for_want_of_the_option(done.output)
    assert done.dir is None or not (done.dir / "index.json").exists()
    assert not done.out.exists()


# ------------------------------------------------------------------ the lattice


def test_the_index_places_every_partition_line(tmp_path: Path) -> None:
    """Point 1 of "Three points from 23c-1's green (3ffae13)": the partition
    record's `row0`, `col0` (the window's first node, possibly negative),
    `spacing` (h) and `origin` (lattice node (0, 0) in the index's `crs`) put
    line `i` at `x = origin_x + (col0 + i dx) h` and line `j` at
    `y = origin_y - (row0 + j dy) h`: exactly the seams the piece files hold."""
    dem, domain, _ = scene(tmp_path)
    done = pf.cut(tmp_path / "run", *args(dem, domain, "--pieces", "9"))
    part = done.index["partition"]
    for key in ("row0", "col0"):
        assert isinstance(part[key], int) and not isinstance(part[key], bool), part
    assert part["spacing"] == pf.H
    ox, oy = part["origin"]
    h, nx, ny, dx, dy = part["spacing"], part["nx"], part["ny"], part["dx"], part["dy"]
    assert (nx, ny) == (3, 3)
    # The window's first node is the lattice node it names.
    assert float((ox + part["col0"] * h - pf.TIE_X) / pf.H).is_integer()
    assert float((pf.TIE_Y - (oy - part["row0"] * h)) / pf.H).is_integer()
    xs = {ox + (part["col0"] + i * dx) * h for i in range(1, nx)}
    ys = {oy - (part["row0"] + j * dy) * h for j in range(1, ny)}
    assert pf.seam_lines(done.pieces) == (xs, ys)


# ------------------------------------------------------------------ point 12


@pytest.mark.parametrize(
    ("option", "target", "out_name"),
    [
        ("--record", "run.json", "x.vtk"),
        ("--stats", "stats.md", "x.vtk"),
        ("--out-edges", "edges.ply", "x.ply"),
    ],
)
def test_options_describing_one_mesh_are_refused_on_a_cut_run(
    tmp_path: Path, option: str, target: str, out_name: str
) -> None:
    """ "23c-2 stopped at dbbe1c6", point 12 and 3 (d): `--record`, like
    `--stats` and `--out-edges`, describes one mesh, and a cut run has none
    until 23d's stitched file. Exit 2 naming the option, and nothing
    written: no `--out`, no pieces directory, no file at the option's path."""
    dem, domain, _ = scene(tmp_path)
    where = tmp_path / "run"
    done = pf.run(
        where, *args(dem, domain, "--pieces", "4", option, str(where / target)), out_name=out_name
    )
    assert done.code == USAGE, done.output
    not_for_want_of_the_option(done.output)
    assert option in "".join(done.output.split())
    assert list(where.iterdir()) == []


@pytest.mark.parametrize(
    ("option", "target", "out_name"),
    [
        ("--record", "run.json", "x.vtk"),
        ("--stats", "stats.md", "x.vtk"),
        ("--out-edges", "edges.ply", "x.ply"),
    ],
)
def test_the_refusal_comes_before_any_piece_is_refined(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, option: str, target: str, out_name: str
) -> None:
    """Ola's ruling of 2026-10-04 (@reviewer round 3's suggestion): once the
    partition is known to have more than one cell, the options that
    describe one mesh are refused before any piece is refined. Spies on
    the seam pass, `refine`, `refine_points` and `basin_run.cut_run`,
    wherever `tin_engine` holds them, must see no call."""
    import sys

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

        monkeypatch.setattr(owner, name, spy)
        for module in [m for key, m in list(sys.modules.items()) if key.startswith("tin_engine")]:
            if getattr(module, name, None) is original:
                monkeypatch.setattr(module, name, spy)
    dem, domain, _ = scene(tmp_path)
    where = tmp_path / "run"
    done = pf.run(
        where, *args(dem, domain, "--pieces", "4", option, str(where / target)), out_name=out_name
    )
    assert done.code == USAGE, done.output
    not_for_want_of_the_option(done.output)
    assert option in "".join(done.output.split())
    assert list(where.iterdir()) == []
    assert calls == [], f"refined before the refusal: {sorted(set(calls))}"
