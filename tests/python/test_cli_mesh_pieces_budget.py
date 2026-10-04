"""B17 (b) through the CLI: `--memory-budget` with a unit.

`docs/increments/23-basin-scale.md`, "23c-2 stopped at dbbe1c6", point 11: a
cut run given `16G` records 17179869184 bytes in the index and writes the
same piece files as the run given the bytes; a refusal is a usage error
naming `--memory-budget` that writes nothing; `--help` shows the default as
`16G` and says units are powers of 1024.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

import pieces_fixtures as pf
from test_cli_mesh_dem import USAGE
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
    by_unit = pf.cut(tmp_path / "unit", *args, "--memory-budget", "16G")
    by_bytes = pf.cut(tmp_path / "bytes", *args, "--memory-budget", "17179869184")
    assert by_unit.index["partition"]["memory_budget"] == 17179869184
    assert by_bytes.index["partition"]["memory_budget"] == 17179869184
    assert files(by_unit) == files(by_bytes)
    assert len(files(by_unit)) == 4


@pytest.mark.parametrize(("value", "says"), [("1.5G", "1536M"), ("16 G", "--memory-budget")])
def test_refusals_write_nothing(tmp_path: Path, value: str, says: str) -> None:
    args = scene(tmp_path)
    done = pf.run(tmp_path / "run", *args, "--memory-budget", value)
    assert done.code == USAGE, done.output
    assert "No such option" not in done.output and "Traceback" not in done.output
    flat = "".join(done.output.split())
    assert "--memory-budget" in flat
    assert "".join(says.split()) in flat
    assert list((tmp_path / "run").iterdir()) == []


def test_the_help_shows_16g_and_powers_of_1024() -> None:
    text = CliRunner().invoke(app, ["mesh", "--help"], env={"COLUMNS": "200"}).output
    flat = " ".join(text.split())
    assert "16G" in flat
    assert "powers of 1024" in flat
