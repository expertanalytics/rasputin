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
from geotiff_fixtures import KARTVERKET, needs_codecs
from test_cli_mesh_dem import invoke
from test_cli_mesh_domain import geojson, quarter_circle
from tin_engine import _core
from tin_engine._core import ChainRole
from tin_engine.cli import DEFAULT_SNAP_SPACING, _constraint_arrays, _domain_chains, _engine
from tin_engine.domain import read_domain
from tin_engine.grid_domain import refine_start_stride, subsample
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
        xy, chains, _ = _domain_chains(read_domain(path, tile.meta), path.name)
    run = _engine(xy, chains, True, DEFAULT_SNAP_SPACING)
    assert run.mesh is not None and run.noded is not None, run.message
    edges, masks = _constraint_arrays(run.mesh, run.noded)
    if min_angle_deg is None:
        out = _core.refine(
            to_core(tile), run.mesh, edges, masks, tolerance=TOLERANCE, threads=threads
        )
    else:
        out = _core.refine(
            to_core(tile), run.mesh, edges, masks, tolerance=TOLERANCE, threads=threads,
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
    code, output = invoke(*args, *extra)
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
