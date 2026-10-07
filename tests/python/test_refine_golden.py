"""T3 of increment 18: ``refine``'s output is unchanged by the row-span scan.

`docs/increments/18-row-span-scan.md`, T3 and R4. Under Ola's C1 (a) every
per-triangle scan result is bit-identical to increment 17's, so the refined
mesh is too, for every thread count. The digests below were RECORDED FROM
INCREMENT 17'S SCAN, at `e0578e5` (the design commit, whose production tree is
increment 17's), before any production change on this branch. No commit may
update them to agree with new code: a changed digest is a changed mesh.

The digest is SHA-256 over the raw bytes of every array ``refine`` returns --
vertices, z, valid, triangles, edges, masks, in that order, each prefixed with
its dtype and shape -- plus the four counters. The inputs are built exactly as
``rasputin mesh --dem <tile> [--domain quarter.geojson] --tolerance 1`` builds
them, through the CLI's own helpers.

Increment 20 (`docs/increments/20-start-quality.md`, Q7 and R11) re-runs them
with the start-quality pass off: the binding's ``min_angle_deg=0`` given
explicitly, and ``rasputin mesh --start-min-angle 0`` itself, whose ``refine``
call is captured and digested. Both must still give increment 17's digests;
the CLI's default (25) must not. Since increment 20b the CLI also turns
constraint feet on by default, so the pre-20 output needs both
``--start-min-angle 0`` and ``--no-constraint-feet``; the binding's default
(``constraint_feet=False``) already leaves them off.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import tin_engine.cli as cli
from cli_driver import geojson, invoke
from geotiff_fixtures import KARTVERKET, needs_codecs
from test_cli_mesh_domain import quarter_circle
from tin_engine import _core
from tin_engine._core import ChainRole
from tin_engine.chains import start_chains
from tin_engine.cli import DEFAULT_SNAP_SPACING, ROLES, _constraint_arrays, _engine
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.grid_domain import refine_start_stride, subsample
from tin_engine.io.domain_file import read_domain
from tin_engine.io.geotiff import decode_dem
from tin_engine.raster import to_core

TOLERANCE = 1.0

# Recorded from increment 17's build; see the module docstring.
GOLDEN = {
    "tile": "adb663588f182be18b5f4c6cbd2ba5ae00d3063b7024976439851b84537c3334",
    "quarter_circle": "05cab2249cd80a9f70e294463ee6bbb1f154948747928784b471d3db07e4f102",
}


def digest(out: _core.RefineOutcome) -> str:
    h = hashlib.sha256()
    for a in (out.vertices, out.z, out.valid, out.triangles, out.edges, out.masks):
        arr = np.ascontiguousarray(a)
        h.update(f"{arr.dtype.str}{arr.shape}".encode())
        h.update(arr.tobytes())
    h.update(f"{out.rounds} {out.inserted} {out.flips} {out.uncovered} {out.carved}".encode())
    h.update(np.float64(out.max_error).tobytes())
    return h.hexdigest()


def refined(
    case: str, tmp_path: Path, threads: int, min_angle_deg: float | None = None
) -> _core.RefineOutcome:
    with KARTVERKET.open("rb") as stream:
        tile = decode_dem(stream)
    if case == "tile":
        xy, ring = subsample(tile.meta, refine_start_stride(tile.meta))
        chains = [(ring, ChainRole.Outer, 0)]
    else:
        path = geojson(tmp_path / "quarter.geojson", quarter_circle())
        # 16b R2: the domain half of the chains moved from `cli._domain_chains`
        # to `chains.start_chains` (test amendment in 16b-1/2's red step,
        # `e99c8ea`); the digest pins that the move changed nothing.
        started = start_chains(read_domain(path), (), DEFAULT_VOCABULARY)
        xy = started.vertices
        chains = [([int(i) for i in c], ROLES[role], int(m)) for c, role, m in started.chains]
    run = _engine(xy, chains, True, DEFAULT_SNAP_SPACING)
    assert run.mesh is not None and run.noded is not None, run.message
    edges, masks = _constraint_arrays(run.mesh, run.noded)
    if min_angle_deg is None:
        out = _core.refine(
            to_core(tile), run.mesh, edges, masks, tolerance=TOLERANCE, threads=threads
        )
    else:
        out = _core.refine(
            to_core(tile),
            run.mesh,
            edges,
            masks,
            tolerance=TOLERANCE,
            threads=threads,
            min_angle_deg=min_angle_deg,
        )
    assert out.ok(), out.message
    return out


@needs_codecs
@pytest.mark.parametrize("threads", [1, 3, 0])
@pytest.mark.parametrize("case", sorted(GOLDEN))
def test_refine_matches_the_increment_17_digest(case: str, threads: int, tmp_path: Path) -> None:
    assert digest(refined(case, tmp_path, threads)) == GOLDEN[case]


@needs_codecs
@pytest.mark.parametrize("threads", [1, 0])
@pytest.mark.parametrize("case", sorted(GOLDEN))
def test_min_angle_0_matches_the_increment_17_digest(
    case: str, threads: int, tmp_path: Path
) -> None:
    assert digest(refined(case, tmp_path, threads, min_angle_deg=0.0)) == GOLDEN[case]


def _cli_outcome(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> _core.RefineOutcome:
    """The RefineOutcome ``rasputin mesh`` itself computes for ``case``."""
    seen: list[_core.RefineOutcome] = []
    real = cli.refine

    def spy(*args: Any, **kwargs: Any) -> _core.RefineOutcome:
        out = real(*args, **kwargs)
        seen.append(out)
        return out

    monkeypatch.setattr(cli, "refine", spy)
    args = ["--dem", str(KARTVERKET), "--tolerance", "1", "--out", str(tmp_path / "x.vtk")]
    if case == "quarter_circle":
        args += ["--domain", str(geojson(tmp_path / "quarter.geojson", quarter_circle()))]
    code, output = invoke("mesh", *args, *extra)
    assert code == 0, output
    (out,) = seen
    return out


@needs_codecs
@pytest.mark.parametrize("case", sorted(GOLDEN))
def test_the_cli_with_start_min_angle_0_matches_the_digest(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Pre-20 output needs both of the CLI's post-17 passes off: the start-quality
    # pass (increment 20) and constraint feet (increment 20b, on by default per
    # its R9; at 1 m they add 6 vertices to the quarter circle).
    out = _cli_outcome(
        case, tmp_path, monkeypatch, "--start-min-angle", "0", "--no-constraint-feet"
    )
    assert digest(out) == GOLDEN[case]


@needs_codecs
def test_the_cli_default_changes_the_domain_mesh(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The converse, so the test above cannot pass with a pass that never runs."""
    out = _cli_outcome("quarter_circle", tmp_path, monkeypatch)
    assert digest(out) != GOLDEN["quarter_circle"]


# ---------------------------------------------------------------- increment 27

#: S7 of `docs/increments/27-node-sampling.md`: ``rasputin mesh --dem
#: KARTVERKET --tolerance 1`` with default flags (start quality 25, constraint
#: feet on), through ``_cli_outcome``'s spy. RECORDED FROM THE PRE-CHANGE
#: PRODUCTION CODE in increment 27's red commit (parent 7eb0b79), before
#: ``RasterGeometry::node_at`` existed. 27 changes ``raster::bilinear`` only at
#: a point that is a node by refine's own test, and refine never calls
#: ``bilinear`` there ("Which paths change"), so no commit may update these.
GOLDEN_DEFAULT_FLAGS = {
    "tile": "a8e8720d37147f3e18d1702362c94418aec87a16cbbd2491d83315f94627060f",
    "quarter_circle": "e2d575076603dbb48c78fe036397a35a4286fef35e7c606ecd4df5ae90b5fd32",
}

#: S8: SHA-256 of the whole ``.vtk`` that ``rasputin mesh --dem KARTVERKET``
#: writes at its default stride. The file carries no rasputin version (only
#: the format's own ``# vtk DataFile Version 4.2`` header) and ``dem_source``
#: is the file's basename, so the whole file is hashed, FieldData included:
#: ``nodata_vertices_removed`` (397) is pinned with the mesh. Recorded from the
#: pre-change code, like S7. Every stride vertex is a node, and on this tile
#: every valid node already samples to its own value bit for bit and no valid
#: node is refused ("What changes in the stride output, exactly"), so the file
#: must not change.
GOLDEN_STRIDE_VTK = "2e0bb5dabc97970babeb0fc4e1f8a035c6aa1e6e6ec20ea316da67ef3b01e633"


@needs_codecs
@pytest.mark.parametrize("case", sorted(GOLDEN_DEFAULT_FLAGS))
def test_the_cli_default_flags_digest_is_unchanged_by_node_sampling(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert digest(_cli_outcome(case, tmp_path, monkeypatch)) == GOLDEN_DEFAULT_FLAGS[case]


@needs_codecs
def test_the_kartverket_stride_vtk_is_unchanged_by_node_sampling(tmp_path: Path) -> None:
    out = tmp_path / "x.vtk"
    # `--ascii` since increment 31 made binary the default: the hash was
    # recorded from the text file, and a hash of it still pins the same mesh.
    code, output = invoke("mesh", "--dem", str(KARTVERKET), "--ascii", "--out", str(out))
    assert code == 0, output
    assert hashlib.sha256(out.read_bytes()).hexdigest() == GOLDEN_STRIDE_VTK
