"""`--piece-budget` through the CLI: a size in master's one grammar.

`docs/increments/23-basin-scale.md`, "23c-2 meets master", R2 and D1: the
piece size is `--piece-budget`, parsed by `resources.parse_bytes` (a whole
number with K, M, G or T, 1024-based; `test_resources.py` covers the grammar
itself). A cut run given `16G` records `piece_budget` 17179869184 in the index
and writes the same piece files as the run given the bytes; a refusal is a
usage error naming `--piece-budget` that writes nothing, `16GB` included
(the two-letter units are not master's grammar); `--help` shows the default
as `16G`.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

import pieces_fixtures as pf
from cli_helpers import USAGE
from tin_engine.cli import app


def scene(tmp_path: Path) -> tuple[str, ...]:
    dem = pf.write_dem(tmp_path / "dem" / "dem.tif", pf.hills())
    domain = pf.write_domain(tmp_path / "sq.geojson", pf.SQUARE)
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", "1", "--pieces", "4")


def files(done: pf.Cut) -> dict[str, bytes]:
    assert done.dir is not None
    return {p.name: p.read_bytes() for p in sorted(done.dir.iterdir()) if p.suffix == ".vtk"}


def test_16g_is_the_byte_value(tmp_path: Path) -> None:
    args = scene(tmp_path)
    by_unit = pf.cut(tmp_path / "unit", *args, "--piece-budget", "16G")
    by_bytes = pf.cut(tmp_path / "bytes", *args, "--piece-budget", "17179869184")
    assert by_unit.index["partition"]["piece_budget"] == 17179869184
    assert by_bytes.index["partition"]["piece_budget"] == 17179869184
    assert "memory_budget" not in by_unit.index["partition"]
    assert files(by_unit) == files(by_bytes)
    assert len(files(by_unit)) == 4


@pytest.mark.parametrize("value", ["1.5G", "16 G", "16GB"])
def test_refusals_write_nothing(tmp_path: Path, value: str) -> None:
    args = scene(tmp_path)
    done = pf.run(tmp_path / "run", *args, "--piece-budget", value)
    assert done.code == USAGE, done.output
    assert "No such option" not in done.output and "Traceback" not in done.output
    flat = "".join(done.output.split())
    assert "--piece-budget" in flat
    assert "isnotanumberofbytes" in flat, done.output
    assert list((tmp_path / "run").iterdir()) == []


def test_the_help_shows_16g() -> None:
    text = CliRunner().invoke(app, ["mesh", "--help"], env={"COLUMNS": "200"}).output
    flat = " ".join(text.split())
    assert "--piece-budget" in flat
    assert "Default 16G, the same on every machine" in flat
