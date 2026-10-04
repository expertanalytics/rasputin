"""DC10, the reprojected path cut into four pieces.

`docs/increments/23-basin-scale.md`, "Tests @tester can write red", 23c,
DC10, and "Tolerance, the final check and the constraint check points per
piece": a synthetic geographic tile (15c-2's G6 scene: EPSG:4326 source, a
domain in EPSG:4674, `--out-crs EPSG:31983`) cut with `--pieces 4`. 15c's
independent final check (`geographic_fixtures.source_check`, G6's oracle)
finds 0 source nodes over tolerance on the union of the piece files, and no
piece counts a source node on a seam (`on_frozen` 0: projected nodes landing
exactly on a lattice line are not expected).

The partition is on the target grid's lattice (J6: node `(R, K)` at
`(K h, -R h)`), so every seam lies on `x = K h` or `y = -R h`, whole metres
at h = 30.

Conformity, a valid union and constrained Delaunay with the seams are
checked here too (section 3D: both oracles on every path). The union's area
is the run's noded outline's to 1e-9, and that outline is the domain moved
by pyproj to within the noder's snap (DC3 as corrected under "23c-2 stopped
at dbbe1c6", 1).

PINNED HERE (handback, choice 11): each index entry's `counts` carries
`on_frozen`, the piece's phase 2 count of source nodes on a frozen edge.

Split under "23c-2 stopped at dbbe1c6", 2: DC10's cases run at
`--start-min-angle 0`; at the default start quality, increment 20's pass
can insert a NoData node, which trim then carves into a hole (the uncut path
has the same defect). Those two cases are `xfail(strict=True)` naming the
fix, so they report themselves once it lands.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import Polygon

import pieces_fixtures as pf
from geographic_fixtures import (
    ANADEM_STEP,
    LAT0,
    LON0,
    NODATA,
    geographic_tile_tiff,
    nodes,
    project,
    project_ring,
    source_check,
)
from test_cli_mesh_geographic import DOMAIN_RC, TARGET, lonlat_ring, tile_array, write_geojson

H = 30.0
TOLERANCE = 1.0


def cut_reprojected(where: Path, *extra: str) -> pf.Cut:
    dem = where / "anadem_like.tif"
    dem.write_bytes(geographic_tile_tiff(tile_array()).getvalue())
    ring = project_ring("EPSG:4326", "EPSG:4674", lonlat_ring(DOMAIN_RC))
    domain = write_geojson(where / "domain.geojson", ring, "EPSG:4674")
    return pf.cut(
        where / "run",
        *("--dem", str(dem), "--domain", str(domain)),
        *("--out-crs", TARGET, "--tolerance", repr(TOLERANCE), "--pieces", "4", *extra),
    )


@pytest.fixture(scope="module")
def done(tmp_path_factory: pytest.TempPathFactory) -> pf.Cut:
    """DC10's subject, the cut on the reprojected path, with the start-quality
    pass off ("23c-2 stopped at dbbe1c6", 2)."""
    return cut_reprojected(tmp_path_factory.mktemp("reprojected"), "--start-min-angle", "0")


@pytest.fixture(scope="module")
def with_quality(tmp_path_factory: pytest.TempPathFactory) -> pf.Cut:
    """The same at the default start quality (25 degrees)."""
    return cut_reprojected(tmp_path_factory.mktemp("reprojected-quality"))


def target_domain() -> Polygon:
    ring = project_ring("EPSG:4326", "EPSG:4674", lonlat_ring(DOMAIN_RC))
    return Polygon(project_ring("EPSG:4674", TARGET, ring))


def test_four_pieces_in_the_target_crs(done: pf.Cut) -> None:
    assert done.index["crs"] == TARGET
    assert sorted(p.id for p in done.pieces) == [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)]
    for piece in done.pieces:
        assert piece.vtk.field_data["crs"].values[0] == TARGET


def test_seams_lie_on_the_target_lattice(done: pf.Cut) -> None:
    xs, ys = pf.seam_lines(done.pieces)
    assert xs and ys
    assert all(float(x / H).is_integer() for x in xs), xs
    assert all(float(y / H).is_integer() for y in ys), ys
    for piece in done.pieces:
        for a, b in pf.seam_edges(piece):
            assert a[0] == b[0] or a[1] == b[1], (piece.id, a, b)


def test_no_source_node_over_tolerance_on_the_union(done: pf.Cut) -> None:
    """G6's independent final check, over the union of the piece files."""
    u = pf.union(done.pieces)
    array = tile_array().ravel().astype(np.float64)
    lon, lat = nodes(LON0, LAT0, ANADEM_STEP, ANADEM_STEP, tile_array().shape)
    valid = array != NODATA
    xy = project("EPSG:4326", TARGET, lon[valid], lat[valid])
    found = source_check(u.points, u.triangles, xy, array[valid], target_domain(), TOLERANCE, H)
    assert found.nodes > 1000
    assert (found.over_interior, found.over_strip) == (0, 0), found


def test_no_piece_counts_a_source_node_on_a_seam(done: pf.Cut) -> None:
    assert done.dir is not None
    index = json.loads((done.dir / "index.json").read_text())
    counts = [entry["counts"]["on_frozen"] for entry in index["pieces"]]
    assert counts == [0] * len(index["pieces"])


def test_conformity_a_valid_union_and_constrained_delaunay(done: pf.Cut) -> None:
    assert pf.conformity_problems(done.pieces) == []
    u = pf.union(done.pieces)
    assert u.z_conflicts == []
    assert pf.validity_problems(u, target_domain(), done.noded) == []
    assert pf.overlapping_pairs(u) == []
    assert pf.delaunay_violations(u.points[:, :2], u.triangles, u.constrained) == []


# ------------------------------------------------------------------ at the default start quality

QUALITY_FIX = (
    "increment 20's start-quality pass inserts NoData nodes and trim leaves holes; "
    'fixed by "Fix: the start-quality pass skips NoData nodes" in 20-start-quality.md '
    "(improve skips void nodes, skipped_void)"
)


@pytest.mark.xfail(strict=True, reason=QUALITY_FIX)
def test_at_the_default_start_quality_no_source_node_over_tolerance(with_quality: pf.Cut) -> None:
    u = pf.union(with_quality.pieces)
    array = tile_array().ravel().astype(np.float64)
    lon, lat = nodes(LON0, LAT0, ANADEM_STEP, ANADEM_STEP, tile_array().shape)
    valid = array != NODATA
    xy = project("EPSG:4326", TARGET, lon[valid], lat[valid])
    found = source_check(u.points, u.triangles, xy, array[valid], target_domain(), TOLERANCE, H)
    assert (found.over_interior, found.over_strip) == (0, 0), found


@pytest.mark.xfail(strict=True, reason=QUALITY_FIX)
def test_at_the_default_start_quality_a_valid_union(with_quality: pf.Cut) -> None:
    u = pf.union(with_quality.pieces)
    assert pf.validity_problems(u, target_domain(), with_quality.noded) == []
