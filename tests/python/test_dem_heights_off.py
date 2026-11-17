"""Without `--dem-heights` every file is byte for byte master's, and the
land-cover overlap keeps its bytes after the shared clip is refactored.

`docs/increments/137-dem-stats.md`, invariants 1 and 6, test 12 of section
5. The digests were recorded by `@tester` at commit 98fe80c3 (public master
before the height shares), in a worktree whose code was that commit's, by
the functions below, three runs alike. The `.gpkg` digest is
`height_fixtures.gpkg_digest`'s, which leaves out the run rows the host
fills (`command`, `threads`, `rasputin_version`). The cut, Velhas and
`class_areas` digests were re-recorded the same way at f4ba684c, when the
core began building with `-ffp-contract=off`; the steps digests did not move.

A recorded digest pins a platform as well as a behaviour: every digest here
is keyed by (`sys.platform`, `platform.machine()`) and a platform without
one skips, saying so, as `R11_SHA256` in `test_cli_catchment_geographic.py`
does. The design asks for `class_areas`'s arrays on every platform; they are
keyed too, because they were recorded on one platform only. Where they
skip, invariant 6 rests on the land-cover suite, which runs everywhere,
and invariant 1 on `TestOffWritesNothing` in `test_cli_mesh_dem_heights.py`.
"""

from __future__ import annotations

import platform
import sys
from pathlib import Path

import numpy as np
import pytest

from class_fixtures import LAT0, STEP
from cli_helpers import invoke
from fraction_fixtures import geod_row_areas, mixed_grid, square_triangulation
from geotiff_fixtures import micro_tiff
from height_fixtures import VELHAS, gpkg_digest, sha256
from piece_spies import square_scene
from pieces_fixtures import cut
from tin_engine import _core

HOST = (sys.platform, platform.machine())
DARWIN_ARM64 = ("darwin", "arm64")

#: Test 1's command, `--stride 10 --tolerance 10`, no flag.
STEPS_DIGESTS = {
    DARWIN_ARM64: {
        "steps.vtk": "78246e74b4bd60f16d115b533b91d886b9c40e649cd66c8f3b4822d57328fbc5",
        "steps.ply": "e2dc70ecd727ee88f023abaed9bb5f52515cee483d9ea27f663068042a049923",
        "steps.gpkg": "dc39e992a286dc131feb15bb34e91cca8d2828d8c528d2b248fed4fff50a8e7f",
    },
}
#: `square_scene` (test 10), `--pieces 4`, no flag: each piece file.
CUT_DIGESTS = {
    DARWIN_ARM64: {
        "0-0-0.vtk": "a26b997df80411d4f5fa9937743e8f783eadf07018f95ddb96c11c48267b294b",
        "0-1-0.vtk": "f0ac2466c2a1661308fda0a5f5bf1cb2dc127a852228ad568354efade0c48488",
        "1-0-0.vtk": "8a6df047a0c9c41f7155a37b3c8c47002ba2ec483ca619719e0947a6e4a838eb",
        "1-1-0.vtk": "b010d0bdca1e158881102c5285fa8195886fff41dfbe02408fe188e05cbb5ca2",
    },
}
#: Velhas (test 11), `--out-crs EPSG:31983 --tolerance 5`, no flag.
VELHAS_DIGEST = {
    DARWIN_ARM64: "21304ae625b979716795105a297bcc09a316aeecd75fc57677a47c39d44a8e52",
}
#: `class_areas` on `test_land_cover_overlap.py`'s `TestThreads` input, one
#: thread: the offsets, classes and areas arrays' bytes.
CLASS_AREAS_DIGESTS = {
    DARWIN_ARM64: (
        "684020a1dcffa7c5a26809bf69a63da102d344db15387c24ee5d25298b4e545b",
        "86f38dd12bdfac6fc13dbe95ac0c63ed1047b4405b3f9532c72b81e930036cd7",
        "95d1ff77edd298e6263d345956513544b96911db9d436f06e4513394c2e0cb46",
    ),
}
#: `TestThreads`' input, restated: three seeds of 400 points on a 20 x 20 grid.
EXTENT = 20
SEEDS = (26, 27, 28)


def recorded[T](table: dict[tuple[str, str], T]) -> T:
    if HOST not in table:
        pytest.skip(f"no digest recorded at 98fe80c3 for {HOST}; recorded for {sorted(table)}")
    return table[HOST]


def write_steps(where: Path) -> Path:
    a = np.zeros((11, 11), np.float32)
    a[:, 1:3] = 2.0
    path = where / "steps.tif"
    path.write_bytes(micro_tiff(a, scale=(10.0, 10.0, 0.0)).getvalue())
    return path


def steps_digest(where: Path, out: str) -> str:
    target = where / out
    code, output = invoke(
        "mesh",
        *("--dem", str(write_steps(where)), "--stride", "10", "--tolerance", "10"),
        *("--out", str(target)),
    )
    assert code == 0, output
    return gpkg_digest(target) if target.suffix == ".gpkg" else sha256(target.read_bytes())


def cut_digests(where: Path) -> dict[str, str]:
    done = cut(where, *square_scene(where), "--pieces", "4")
    return {p.path.name: sha256(p.path.read_bytes()) for p in done.pieces}


def velhas_digest(where: Path) -> str:
    target = where / "velhas.vtk"
    code, output = invoke(
        "mesh",
        *("--dem", str(VELHAS / "anadem_velhas.tif")),
        *("--domain", str(VELHAS / "catchment.geojson")),
        *("--out-crs", "EPSG:31983", "--tolerance", "5", "--out", str(target)),
    )
    assert code == 0, output
    return sha256(target.read_bytes())


def class_areas_digests() -> tuple[str, str, str]:
    corners = np.concatenate([square_triangulation(s, EXTENT, points=400) for s in SEEDS])
    grid = mixed_grid(EXTENT, EXTENT)
    weights = geod_row_areas(LAT0, STEP, EXTENT)
    offsets, classes, areas = _core.class_areas(
        np.asarray(corners, np.float64), grid, weights, threads=1
    )
    assert (offsets.dtype, classes.dtype, areas.dtype) == (np.uint32, np.uint8, np.float64)
    a, b, c = (sha256(x.tobytes()) for x in (offsets, classes, areas))
    return a, b, c


@pytest.mark.parametrize("out", ["steps.vtk", "steps.ply", "steps.gpkg"])
def test_steps_without_the_flag_are_master_bytes(tmp_path: Path, out: str) -> None:
    want = recorded(STEPS_DIGESTS)[out]
    assert steps_digest(tmp_path, out) == want


def test_the_cut_run_without_the_flag_is_master_bytes(tmp_path: Path) -> None:
    want = recorded(CUT_DIGESTS)
    assert cut_digests(tmp_path) == want


def test_velhas_without_the_flag_is_master_bytes(tmp_path: Path) -> None:
    want = recorded(VELHAS_DIGEST)
    assert velhas_digest(tmp_path) == want


def test_class_areas_keeps_its_bytes() -> None:
    """Invariant 6: the visitor does today's arithmetic in today's order."""
    want = recorded(CLASS_AREAS_DIGESTS)
    assert class_areas_digests() == want
