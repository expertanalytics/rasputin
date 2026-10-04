"""tester.md section 3A on a cut run: NoData and NaN on a seam, a narrow
corridor crossed by a seam, and tolerance 0.

`docs/increments/23-basin-scale.md`, "Degeneracy policy": NoData on a seam
(a check point whose stencil touches NoData is skipped; a seam vertex
without a height is invalid and the triangles around it are carved in each
piece), a seam edge shorter than a cell (the corridor between a hole and the
outline, 0.4 m wide where the horizontal seam crosses it: one check point at
most, N15), and "at T = 0 every node is a vertex".

Where NoData carves the mesh, DC3's area and boundary checks do not apply
(the carved holes are neither), so the union is checked for what still
holds: each seam edge in exactly two files bit for bit, no edge in three
triangles, no overlap, constrained Delaunay with the seams, and section 3D's
tolerance oracle in its own words, every valid node in a triangle with three
valid vertices (`pieces_fixtures.located_over`).

HOW THIS FILE GOES RED: `--pieces` is not an option of `mesh`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon

import pieces_fixtures as pf
from geotiff_fixtures import micro_tiff
from pieces_fixtures import rel

SENTINEL = -32767.0


def args(dem: Path, domain: Path, tolerance: str = "1", *extra: str) -> tuple[str, ...]:
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", tolerance, *extra)


def learn_lines(where: Path, dem: Path, domain: Path) -> tuple[float, float]:
    first = pf.cut(where, *args(dem, domain), "--pieces", "4")
    xs, ys = pf.seam_lines(first.pieces)
    assert len(xs) == 1 and len(ys) == 1, (xs, ys)
    return xs.pop(), ys.pop()


@pytest.mark.parametrize("kind", ["sentinel", "nan"])
def test_nodata_on_and_beside_a_seam(tmp_path: Path, kind: str) -> None:
    """Voids on the vertical seam, beside it, and at the cell corner."""
    array = pf.hills()
    domain = pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE)
    clean = pf.write_dem(tmp_path / "clean" / "dem.tif", array)
    x, y = learn_lines(tmp_path / "first", clean, domain)
    col, row = int((x - pf.TIE_X) / pf.H), int((pf.TIE_Y - y) / pf.H)
    holed = array.copy()
    void = SENTINEL if kind == "sentinel" else np.float32(np.nan)
    holed[row - 12 : row - 8, col] = void  # on the seam
    holed[row + 6, col + 1] = void  # beside it
    holed[row, col] = void  # the corner
    path = tmp_path / "holed" / "dem.tif"
    path.parent.mkdir()
    nodata = "-32767" if kind == "sentinel" else None
    path.write_bytes(micro_tiff(holed, scale=(pf.H, pf.H, 0.0), nodata=nodata).getvalue())

    done = pf.cut(tmp_path / "run", *args(path, domain), "--pieces", "4")
    assert pf.seam_lines(done.pieces) == ({x}, {y})
    assert pf.conformity_problems(done.pieces) == []
    u = pf.union(done.pieces)
    assert u.z_conflicts == []
    assert np.isfinite(u.points).all()
    assert all(len(held) <= 2 for held in pf.sides(u.triangles).values())
    assert pf.overlapping_pairs(u) == []
    assert pf.delaunay_violations(u.points[:, :2], u.triangles, u.constrained) == []
    valid = np.isfinite(holed) & (holed != SENTINEL)
    located, over = pf.located_over(u, holed, valid, 1.0)
    assert located > 900 and over == 0
    # No vertex sits on a void node: the carving happened in each piece.
    xs_v, ys_v, _ = pf.node_xy(holed)
    voids = {
        (float(a), float(b))
        for a, b in zip(xs_v[~valid.ravel()], ys_v[~valid.ravel()], strict=True)
    }
    assert not voids & {(float(p[0]), float(p[1])) for p in u.points}


def test_a_seam_across_a_narrow_corridor(tmp_path: Path) -> None:
    """A hole 0.4 m from the outline's left side, crossed by the horizontal
    seam: the seam edge across the corridor is shorter than a cell."""
    array = pf.hills()
    hole = [rel(30.9, -150.3), rel(29.2, -260.3), rel(149.3, -260.1), rel(150.7, -149.9)]
    outer = Polygon(pf.SQUARE)
    assert outer.contains(Polygon(hole))
    assert 0.0 < outer.exterior.distance(Polygon(hole)) < 0.5
    dem = pf.write_dem(tmp_path / "dem.tif", array)
    domain = pf.write_domain(tmp_path / "corridor.geojson", pf.SQUARE, (hole,))
    done = pf.cut(tmp_path / "run", *args(dem, domain), "--pieces", "4")
    _, ys = pf.seam_lines(done.pieces)
    (y,) = ys
    assert Polygon(hole).bounds[1] < y < Polygon(hole).bounds[3], "the seam misses the hole"
    short = [
        (a, b)
        for piece in done.pieces
        for a, b in pf.seam_edges(piece)
        if a[1] == b[1] == y and abs(a[0] - b[0]) < pf.H
    ]
    assert short, "no seam edge shorter than a cell"
    assert pf.problems(done, array, Polygon(pf.SQUARE, [hole]), 1.0) == []


def test_tolerance_zero_makes_every_node_in_the_domain_a_vertex(tmp_path: Path) -> None:
    array = pf.hills()
    dem = pf.write_dem(tmp_path / "dem.tif", array)
    domain = pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE)
    done = pf.cut(tmp_path / "run", *args(dem, domain, "0"), "--pieces", "4")
    assert done.index["partition"]["bytes_per_node"] == 637
    assert pf.problems(done, array, Polygon(pf.SQUARE), 0.0) == []
    u = pf.union(done.pieces)
    x, y, _ = pf.node_xy(array)
    inside = shapely.contains_xy(Polygon(pf.SQUARE), x, y)
    vertices = {(float(p[0]), float(p[1])) for p in u.points}
    missing = [
        (a, b) for a, b, keep in zip(x, y, inside, strict=True) if keep and (a, b) not in vertices
    ]
    assert missing == []
