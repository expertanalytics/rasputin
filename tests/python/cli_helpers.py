"""One way to drive `rasputin` in the tests: the runner, the output, the inputs.

`docs/increments/python-audit.md`, section 9 (PR T1). Test support only;
nothing in `src_python` uses it. The CLI suites take these from here rather
than from each other, so no test module is a helper library for another.

`invoke` takes the subcommand as its first argument and refuses a name that is
not one: a call written for a mesh-only helper, such as `invoke("--dem", ...)`,
would otherwise run `rasputin --dem ...`, exit 2 with "No such option", and
pass a test that checks only for exit code 2.

The two fixture factories return a pytest fixture; a suite binds it at module
level to the name its tests request (`bumpy = rough_dem(16)`). Bound inside a
class, the fixture would receive the instance as `tmp_path`.
"""

from __future__ import annotations

import io
import json
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import typer
from typer.testing import CliRunner

from geotiff_fixtures import TIE_X, TIE_Y, micro_tiff
from tin_engine.cli import app
from vtkread import VtkFile, read_vtk

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

USAGE = 2
ANSI = re.compile(r"\x1b\[[0-9;]*m")
BOX = re.compile(r"[─-╿]")

#: The subcommands, read from the app rather than restated.
COMMANDS = frozenset(typer.main.get_command(app).commands)


def plain(text: str) -> str:
    """Output as a reader sees it: no colour, no box rule, no line wrapping.

    Typer renders errors inside a Rich panel whose width is the terminal's, so
    a message asserted verbatim would fail on a narrow one and pass on a wide
    one -- a flaky test by construction. Collapsing the box and the whitespace
    makes every assertion independent of where Rich chose to wrap, while
    leaving words and numbers intact, which is all any of them look at.
    """
    return " ".join(BOX.sub(" ", ANSI.sub("", text)).split())


def squashed(text: str) -> str:
    """No whitespace at all: Rich may break a long message anywhere inside its panel."""
    return "".join(text.split())


def invoke(command: str, *args: str) -> tuple[int, str]:
    """`rasputin <command> <args>`: the exit code and the `plain` output."""
    assert command in COMMANDS, f"{command!r} is not a subcommand; they are {sorted(COMMANDS)}"
    result = runner.invoke(app, [command, *args])
    return result.exit_code, plain(result.output)


def ran(command: str, *args: str) -> str:
    """`invoke`, required to succeed; the output."""
    code, output = invoke(command, *args)
    assert code == 0, output
    return output


def refused(
    tmp_path: Path, command: str, *args: str, says: tuple[str, ...], squash: bool = False
) -> str:
    """Refused with exit 2 by the command, saying each of `says`, writing nothing.

    With `squash`, each word is looked for with all whitespace removed from
    both sides, since Rich may break a long word inside its panel.
    """
    target = tmp_path / "refused.vtk"
    code, output = invoke(command, *args, "--out", str(target))
    assert code == USAGE, output
    assert "Traceback" not in output
    assert "No such option" not in output, "refused for the wrong reason"
    for word in says:
        if squash:
            assert squashed(word) in squashed(output), f"{word!r} not in {output!r}"
        else:
            assert word in output, f"{word!r} not in {output!r}"
    assert not target.exists()
    return output


def mesh_to_vtk(tmp_path: Path, *args: str, out: str = "x.vtk") -> VtkFile:
    """Mesh with ``--stats`` beside the file (``x.md``); the mesh read back."""
    target = tmp_path / out
    ran("mesh", *args, "--out", str(target), "--stats", str(target.with_suffix(".md")))
    return read_vtk(target.read_bytes())


def write_tiff(path: Path, stream: io.BytesIO) -> Path:
    path.write_bytes(stream.getvalue())
    return path


# ---------------------------------------------------------------- the domain

UTM33 = "urn:ogc:def:crs:EPSG::25833"
Ring = list[tuple[float, float]]

# micro_tiff's grid: nodes at x = TIE_X + 10 col, y = TIE_Y - 5 row, EPSG:25833,
# PixelIsPoint. 17 rows x 21 cols: x 500 000 .. 500 200, y 6 599 920 .. 6 600 000.
ROWS, COLS = 17, 21

# A square and a hole, every vertex off-node (x not a multiple of 10, y not of 5).
SQUARE: Ring = [
    (TIE_X + 12.3, TIE_Y - 73.3),
    (TIE_X + 187.7, TIE_Y - 72.9),
    (TIE_X + 186.1, TIE_Y - 6.7),
    (TIE_X + 13.9, TIE_Y - 7.1),
]
HOLE: Ring = [
    (TIE_X + 71.1, TIE_Y - 52.7),
    (TIE_X + 72.3, TIE_Y - 28.3),
    (TIE_X + 121.9, TIE_Y - 28.9),
    (TIE_X + 120.7, TIE_Y - 51.1),
]


def geojson(path: Path, outer: Ring, holes: tuple[Ring, ...] = (), crs: str | None = UTM33) -> Path:
    doc: dict[str, Any] = {
        "type": "Polygon",
        "coordinates": [[*r, r[0]] for r in (outer, *holes)],
    }
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(doc))
    return path


# ---------------------------------------------------------------- fixture factories


def rough_dem(seed: int) -> Callable[[Path], Path]:
    """A fixture: ``bumpy.tif``, ROWS x COLS uniform in 0 to 50 m from ``seed``."""

    @pytest.fixture
    def fixture(tmp_path: Path) -> Path:
        array = np.random.default_rng(seed).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
        return write_tiff(tmp_path / "bumpy.tif", micro_tiff(array))

    return fixture


def polygon_file(outer: Ring, holes: tuple[Ring, ...] = ()) -> Callable[[Path], Path]:
    """A fixture: ``square.geojson``, the polygon ``outer`` with ``holes``, in UTM33."""

    @pytest.fixture
    def fixture(tmp_path: Path) -> Path:
        return geojson(tmp_path / "square.geojson", outer, holes)

    return fixture
