"""B17 (b): `--memory-budget` takes a size with a unit, every unit a power of 1024.

`docs/increments/23-basin-scale.md`, "23c-2 stopped at dbbe1c6", point 11
(Ola's ruling of 2026-10-04: "accept human readable formatting for
--memory-budget"). The grammar: a whole number of digits, then an optional
unit straight after it, no space; units `K`, `M`, `G`, `T`, each optionally
followed by `B` or `iB`, all powers of 1024, case-insensitive; a bare number,
or one followed by `B` alone, is bytes. `tin_engine.decompose.parse_memory_budget`
is pure and raises `ValueError`; the CLI turns that into the usage error
(`test_cli_mesh_pieces_budget.py`).

PINNED HERE, where the point is silent (handback): the digits are ASCII
(fullwidth digits are refused, as `int()` alone would accept them); leading
zeros are digits like any other (`016G` is 16 GiB); whitespace before or after
the value is refused as a space inside is (ruled after the red step); a zero
budget is refused
here (no window fits in 0 bytes), while "below one node's bytes", which
needs the tolerance, stays `partition`'s and the CLI's refusal.
"""

from __future__ import annotations

import importlib
from types import ModuleType

import pytest

KIB, MIB, GIB, TIB = 2**10, 2**20, 2**30, 2**40


@pytest.fixture
def dec() -> ModuleType:
    return importlib.import_module("tin_engine.decompose")


ACCEPTED = [
    ("16G", 16 * GIB),
    ("16g", 16 * GIB),
    ("16GB", 16 * GIB),
    ("16gb", 16 * GIB),
    ("16GiB", 16 * GIB),
    ("16gib", 16 * GIB),
    ("16GIB", 16 * GIB),
    ("17179869184", 2**34),
    ("512M", 2**29),
    ("512MiB", 2**29),
    ("1T", 2**40),
    ("1tb", 2**40),
    ("1024K", 2**20),
    ("1024KB", 2**20),
    ("1024kib", 2**20),
    ("100B", 100),
    ("100b", 100),
    ("100", 100),
    ("016G", 16 * GIB),
    ("3G", 3 * GIB),
]


@pytest.mark.parametrize(("text", "expected"), ACCEPTED)
def test_accepted(dec: ModuleType, text: str, expected: int) -> None:
    value = dec.parse_memory_budget(text)
    assert value == expected
    assert type(value) is int


def test_every_unit_is_a_power_of_1024(dec: ModuleType) -> None:
    for unit, power in (("K", 1), ("M", 2), ("G", 3), ("T", 4)):
        for suffix in ("", "B", "iB"):
            assert dec.parse_memory_budget(f"7{unit}{suffix}") == 7 * 1024**power


REFUSED = [
    pytest.param("", id="empty"),
    pytest.param("+16G", id="plus"),
    pytest.param("-16G", id="minus"),
    pytest.param("-1", id="negative-bytes"),
    pytest.param("1.5G", id="decimal"),
    pytest.param("16.0", id="decimal-bytes"),
    pytest.param("16 G", id="inner-space"),
    pytest.param("16 GB", id="inner-space-gb"),
    pytest.param(" 16G", id="leading-space"),
    pytest.param("16G ", id="trailing-space"),
    pytest.param("\t16G", id="leading-tab"),
    pytest.param("G", id="unit-without-number"),
    pytest.param("GiB", id="unit-without-number-gib"),
    pytest.param("16X", id="unknown-unit"),
    pytest.param("16P", id="peta"),
    pytest.param("16PB", id="petabytes"),
    pytest.param("16iB", id="ib-without-unit"),
    pytest.param("16BB", id="double-b"),
    pytest.param("16Gi", id="gi-without-b"),
    pytest.param("0x10", id="hex"),
    pytest.param("1e9", id="exponent"),
    pytest.param("16_000", id="underscore"),
    pytest.param("\uff11\uff16G", id="fullwidth-digits"),
    pytest.param("0", id="zero"),
    pytest.param("0G", id="zero-gib"),
]


@pytest.mark.parametrize("text", REFUSED)
def test_refused(dec: ModuleType, text: str) -> None:
    with pytest.raises(ValueError):
        dec.parse_memory_budget(text)


def test_a_decimal_suggests_the_whole_number_in_a_smaller_unit(dec: ModuleType) -> None:
    with pytest.raises(ValueError, match="1536M"):
        dec.parse_memory_budget("1.5G")
