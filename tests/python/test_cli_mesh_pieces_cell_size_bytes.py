"""CS5: on lattices binary holds exactly, a cut writes the bytes it wrote before (K6).

`docs/increments/23-basin-scale.md`, "23c-4: cutting a DEM whose cell size is
not exact in binary", "Red tests", CS5, the step's invariant-critical test.
A double restores the rules from before the step in place of
`basin_run.seam_box` (`floor` and `ceil` of the ends' fractional indices,
minus and plus one) and `basin_run.corner_heights` (no corner fixed, so every
seam end takes its height from a seam edge), and the real run and the
double's write byte-identical piece files and `index.json`: the square scene
at 10 m and at 9.75 m (a power-of-two fraction), and DC10's reprojected
scene, whose cell is chosen.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

import piece_spies as ps
import pieces_fixtures as pf
from test_cli_mesh_pieces_cell_size import today_box

SLOPE = ("--tolerance-slope", "0.3", "5", "20")


def install_todays_rules(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """`seam_box` and `corner_heights` replaced by the double; the names of
    the ones called, so a double that never ran fails the test."""
    from tin_engine import basin_run

    called: list[str] = []

    def box(cut: Any, snap: float, a: Any, b: Any) -> tuple[int, int, int, int]:
        called.append("seam_box")
        return today_box(cut, snap, a, b)

    def corners(*args: Any, **kwargs: Any) -> dict[int, float | None]:
        called.append("corner_heights")
        return {}

    ps.patch_everywhere(monkeypatch, basin_run, "seam_box", box)
    ps.patch_everywhere(monkeypatch, basin_run, "corner_heights", corners)
    return called


def both(where: Path, monkeypatch: pytest.MonkeyPatch, *args: str) -> tuple[dict[str, bytes], ...]:
    real = pf.cut(where / "real", *args)
    with monkeypatch.context() as mp:
        called = install_todays_rules(mp)
        before = pf.cut(where / "todays-rules", *args)
    assert {"seam_box", "corner_heights"} <= set(called), called
    assert real.dir is not None and before.dir is not None
    return ps.piece_files(real.dir), ps.piece_files(before.dir)


def square(where: Path, h: float) -> tuple[str, ...]:
    dem = pf.write_dem(where / "dem.tif", pf.hills(), h=h)
    domain = pf.write_domain(where / "square.geojson", pf.SQUARE)
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1")


@pytest.mark.parametrize(
    "h, pieces, slope",
    [(10.0, "4", False), (10.0, "9", False), (10.0, "16", False),
     (10.0, "4", True), (10.0, "16", True), (9.75, "16", False)],
    ids=["10m-4", "10m-9", "10m-16", "10m-4-slope", "10m-16-slope", "9.75m-16"],
)  # fmt: skip
def test_cs5_square_scene_bytes_are_kept(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, h: float, pieces: str, slope: bool
) -> None:
    args = (*square(tmp_path, h), "--pieces", pieces, *(SLOPE if slope else ()))
    real, before = both(tmp_path, monkeypatch, *args)
    assert "index.json" in real
    assert len(real) == int(pieces) + 1, sorted(real)
    assert real == before


@pytest.mark.parametrize("quality", ["default", "off"])
def test_cs5_reprojected_scene_bytes_are_kept(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, quality: str
) -> None:
    angle = ("--start-min-angle", "0") if quality == "off" else ()
    args = (*ps.reprojected_scene(tmp_path), "--pieces", "4", *angle)
    real, before = both(tmp_path, monkeypatch, *args)
    assert "index.json" in real
    assert real == before
