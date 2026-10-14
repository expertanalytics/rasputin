"""`--z-unit` is a choice typer itself parses: a bad value is a usage error.

`docs/increments/48-declared-imports.md`, *6. The tests that pin it*, T1, T2
and T4. Every command that reads a DEM takes `--z-unit`; a spelling outside
the three is refused like any other bad option value (exit 2, one message
naming the option, the value and the choices), `--help` shows the choices,
and the spellings the option takes are the ones the unit table holds.

Every case is an argument list and no file is read: the option's value is
converted before click looks for missing required options, so the refusal
comes first even with nothing else given.
"""

from __future__ import annotations

import typing

import pytest

from cli_helpers import USAGE, help_text, plain, runner, squashed
from tin_engine.cli import app
from tin_engine.io import models

#: The commands that take `--z-unit`.
DEM_COMMANDS = ("mesh", "catchment", "station-catchments")
SPELLINGS = ("m", "ft", "us-ft")


@pytest.mark.parametrize("command", DEM_COMMANDS)
def test_an_unknown_unit_is_a_usage_error(command: str) -> None:
    result = runner.invoke(app, [command, "--z-unit", "furlong"])
    said = squashed(plain(result.output))
    assert result.exit_code == USAGE, said
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    for word in ("--z-unit", "furlong", *(f"'{s}'" for s in SPELLINGS)):
        assert word in said, (word, said)


@pytest.mark.parametrize("command", DEM_COMMANDS)
def test_help_shows_the_three_spellings(command: str) -> None:
    said = squashed(help_text(command))
    assert "m|ft|us-ft" in said
    assert "<function>" not in said


def test_the_option_takes_exactly_the_tables_spellings() -> None:
    table = tuple(row[2] for row in models.UNITS.values())
    assert typing.get_args(models.ZSpelling) == table == SPELLINGS
