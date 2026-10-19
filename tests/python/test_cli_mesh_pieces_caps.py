"""DC11 and MP4: the cap does not cut, the budget does not cap, the cap is per piece.

`docs/increments/23-basin-scale.md`, "23c-2 meets master", R1-R3 and
`@tester`'s step 4, and "23c-3: each piece's memory its own", D4 and MP4.
Two sizes on one command, each with one job: `--memory-budget` is the
machine's cap on every DEM window a run assembles (increment 38, section
4.3), which on a cut run is each piece's window, and `--piece-budget`
decides the cut and nothing else.

- (a) A cut run writes the same piece files under `--memory-budget` at the
  largest piece window's size and at ten times it, and with
  `resources.machine_memory` patched small and large while no
  `--memory-budget` is given: the cap never changes a cut (K5). MP4 (a): a
  cap one byte under the run's whole window, which every piece fits, runs
  and writes the same files.
- (b) MP4 (b): a cut run whose largest piece window is over
  `--memory-budget` exits 2 naming `--memory-budget` and that piece's
  figure, writes nothing, and refines no piece (D4).
- (c) An uncut run writes the same bytes with `--piece-budget` given or not,
  as long as the budget does not cut it (K1).

The largest piece window's size is the run's own, read from the refusal the
cap gives at one byte (`the request is R x C nodes, S bytes`), not
re-derived here; the run's whole window is the index's `rows x cols`. The
scene is `pieces_fixtures`' 41 x 41-node float32 DEM; with `--pieces 4` its
run window is 37 x 37 nodes (5,476 bytes) and its largest piece window
21 x 21 (1,764 bytes), measured on 23c-2. "Small" machine memory is twice
the largest piece window (half of it, the default cap, is that window
exactly) and "large" is 2**50 bytes.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import piece_spies as ps
import pieces_fixtures as pf
from cli_helpers import USAGE

REQUEST = re.compile(r"the request is (\d+) x (\d+) nodes, (\d+) bytes at float32")


@pytest.fixture(scope="module")
def scene(tmp_path_factory: pytest.TempPathFactory) -> tuple[str, ...]:
    where = tmp_path_factory.mktemp("caps")
    dem = pf.write_dem(where / "dem.tif", pf.hills())
    domain = pf.write_domain(where / "square.geojson", pf.SQUARE)
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1")


@pytest.fixture(scope="module")
def window(tmp_path_factory: pytest.TempPathFactory, scene: tuple[str, ...]) -> int:
    """The bytes of the largest piece window, from the run's own refusal at 1 byte."""
    done = pf.run(tmp_path_factory.mktemp("probe"), *scene, "--pieces", "4", "--memory-budget", "1")
    found = REQUEST.search(" ".join(done.output.split()))
    assert done.code == USAGE and found is not None, done.output
    rows, cols, size = (int(g) for g in found.groups())
    assert size == rows * cols * 4
    assert (rows, cols, size) == (21, 21, 1764), "not the largest piece window (D4)"
    return size


def files(done: pf.Cut) -> dict[str, bytes]:
    assert done.dir is not None
    return {p.name: p.read_bytes() for p in sorted(done.dir.iterdir())}


@pytest.fixture(scope="module")
def reference(tmp_path_factory: pytest.TempPathFactory, scene: tuple[str, ...]) -> dict[str, bytes]:
    """The cut with neither size given."""
    done = pf.cut(tmp_path_factory.mktemp("reference"), *scene, "--pieces", "4")
    assert len(done.pieces) == 4
    return files(done)


# ------------------------------------------------------------------ (a) the cap does not cut


@pytest.mark.parametrize("times", [1, 10])
def test_the_cap_at_the_window_and_ten_times_it_cuts_the_same(
    tmp_path: Path, scene: tuple[str, ...], window: int, reference: dict[str, bytes], times: int
) -> None:
    done = pf.cut(tmp_path, *scene, "--pieces", "4", "--memory-budget", str(times * window))
    assert files(done) == reference


@pytest.mark.parametrize("machine", ["small", "large"])
def test_the_machines_memory_cuts_the_same(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    scene: tuple[str, ...],
    window: int,
    reference: dict[str, bytes],
    machine: str,
) -> None:
    import tin_engine.resources as resources

    memory = 2 * window if machine == "small" else 2**50
    monkeypatch.setattr(resources, "machine_memory", lambda *_: memory)
    assert resources.Resources.fallback().memory_budget == memory // 2
    done = pf.cut(tmp_path, *scene, "--pieces", "4")
    assert files(done) == reference


# ------------------------------------------------------------------ (b) the cap refuses


def test_mp4a_a_cap_under_the_run_window_that_every_piece_fits_runs(
    tmp_path: Path, scene: tuple[str, ...], reference: dict[str, bytes]
) -> None:
    probe = pf.cut(tmp_path / "probe", *scene, "--pieces", "4")
    part = probe.index["partition"]
    run_window = part["rows"] * part["cols"] * 4
    assert run_window == 5476
    done = pf.cut(tmp_path / "run", *scene, "--pieces", "4", "--memory-budget", str(run_window - 1))
    assert files(done) == reference


def test_mp4b_a_piece_window_over_the_cap_is_refused_before_any_piece_is_refined(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: tuple[str, ...], window: int
) -> None:
    calls = ps.refine_calls(monkeypatch)
    where = tmp_path / "run"
    done = pf.run(where, *scene, "--pieces", "4", "--memory-budget", str(window - 1))
    assert done.code == USAGE, done.output
    flat = " ".join(done.output.split())
    assert "Traceback" not in flat
    assert "--memory-budget" in flat
    assert "21 x 21 nodes, 1764 bytes at float32" in flat
    assert f"{window} bytes at float32, over the memory budget of {window - 1} bytes" in flat
    assert list(where.iterdir()) == []
    assert calls == [], f"refined before the refusal: {sorted(set(calls))}"


# ------------------------------------------------------------------ (c) the budget does not cap


@pytest.mark.parametrize("budget", ["16G", "1T", str(2**40 + 1)])
def test_an_uncut_run_is_unchanged_by_the_piece_budget(
    tmp_path: Path, scene: tuple[str, ...], budget: str
) -> None:
    today = pf.run(tmp_path / "today", *scene, "--binary")
    given = pf.run(tmp_path / "given", *scene, "--binary", "--piece-budget", budget)
    assert (today.code, given.code) == (0, 0), (today.output, given.output)
    assert (today.dir, given.dir) == (None, None)
    assert given.out.read_bytes() == today.out.read_bytes()
