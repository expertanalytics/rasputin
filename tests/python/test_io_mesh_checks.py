"""`tin_engine.io.mesh_checks`: the mesh writers' two shared checks.

`docs/increments/python-audit-pr-d.md`, section 3 (the design) and section 9,
red test 1. Before PR D, `io/ply.py` and `io/vtk_legacy.py` each carried a
copy of the text gate and of the int32 check, and the two text gates had
drifted apart in wording. This suite pins the one copy: what it accepts, what
it refuses, and the exact sentence of each refusal, since both writers and
`rasputin mesh --crs` pass that sentence on to a person.

Committed red at `8dcfa2a`: the module did not exist. It is fetched inside a
fixture, as `test_run_record.py` fetches its module, so its absence failed
each test here with `ModuleNotFoundError` and left the rest of the session
collecting (a collection error would have stopped the whole run). It landed
in `7346a0e`. Not invariant-critical: no mutation round (section 9, "Lean").

Pinned here, beyond the design's wording (section 9, "Pinned by the red step
(`8dcfa2a`), ruled", lists each): a value
with both a control character and a non-ASCII one is refused as a control
character; the refusal names the first offending character; a C1 control
(U+0080 to U+009F) is refused as non-ASCII, since it is not ASCII; space and
`~` are the edges of what passes; `check_int32` judges an unsigned array by
its values, not its dtype.
"""

from __future__ import annotations

import importlib
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest

INT32_MIN = -(2**31)
INT32_MAX = 2**31 - 1

#: Every character the design names a control character, by kind: NUL, the
#: tab, the two line terminators a PLY or VTK reader splits on, ESC (a
#: terminal escape), and DEL. U+001F is the last below the space.
CONTROLS = ["\x00", "\t", "\n", "\r", "\x1b", "\x1f", "\x7f"]


@pytest.fixture
def checks() -> ModuleType:
    return importlib.import_module("tin_engine.io.mesh_checks")


class TestCheckedAscii:
    def test_plain_text_is_its_ascii_bytes(self, checks: ModuleType) -> None:
        assert checks.checked_ascii("crs EPSG:25833", "comment") == b"crs EPSG:25833"

    def test_the_empty_text_is_empty_bytes(self, checks: ModuleType) -> None:
        assert checks.checked_ascii("", "string") == b""

    def test_space_and_tilde_bound_what_passes(self, checks: ModuleType) -> None:
        # U+0020 is the first character above the controls and U+007E the
        # last below DEL: the two edges of the accepted range.
        assert checks.checked_ascii(" ~", "string") == b" ~"

    @pytest.mark.parametrize("what", ["comment", "string"])
    @pytest.mark.parametrize("bad", CONTROLS, ids=[f"U+{ord(c):04X}" for c in CONTROLS])
    def test_a_control_character_is_refused_naming_what_and_the_character(
        self, checks: ModuleType, bad: str, what: str
    ) -> None:
        with pytest.raises(ValueError) as refusal:
            checks.checked_ascii(f"crs x{bad}forged", what)
        assert str(refusal.value) == f"a {what} may not contain a control character; got {bad!r}"

    def test_the_first_control_character_is_the_one_named(self, checks: ModuleType) -> None:
        with pytest.raises(ValueError) as refusal:
            checks.checked_ascii("a\rb\nc", "comment")
        assert str(refusal.value) == "a comment may not contain a control character; got '\\r'"

    @pytest.mark.parametrize("what", ["comment", "string"])
    def test_a_non_ascii_character_is_refused_naming_it_and_the_value(
        self, checks: ModuleType, what: str
    ) -> None:
        value = "ETRS89 60\N{DEGREE SIGN}N"
        with pytest.raises(ValueError) as refusal:
            checks.checked_ascii(value, what)
        assert str(refusal.value) == f"a {what} must be ASCII; got '\N{DEGREE SIGN}' in {value!r}"

    def test_the_first_non_ascii_character_is_the_one_named(self, checks: ModuleType) -> None:
        value = "d\N{LATIN SMALL LETTER E WITH ACUTE}m \N{LATIN SMALL LETTER O WITH STROKE}"
        with pytest.raises(ValueError) as refusal:
            checks.checked_ascii(value, "string")
        assert str(refusal.value) == f"a string must be ASCII; got '\xe9' in {value!r}"

    def test_a_c1_control_is_refused_as_non_ascii(self, checks: ModuleType) -> None:
        # U+0085 (NEL) ends a line for a Unicode-aware reader. It is outside
        # ASCII, so the ASCII refusal is the one that meets it.
        with pytest.raises(ValueError) as refusal:
            checks.checked_ascii("x\x85y", "comment")
        assert str(refusal.value) == "a comment must be ASCII; got '\\x85' in 'x\\x85y'"

    def test_a_control_character_is_named_before_a_non_ascii_one(self, checks: ModuleType) -> None:
        with pytest.raises(ValueError) as refusal:
            checks.checked_ascii("\N{DEGREE SIGN}\r", "comment")
        assert str(refusal.value) == "a comment may not contain a control character; got '\\r'"


class TestCheckInt32:
    @pytest.mark.parametrize(
        "values",
        [
            np.array([INT32_MAX, INT32_MIN], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([0, INT32_MAX], dtype=np.uint64),
        ],
        ids=["both-bounds", "empty", "unsigned-at-the-top"],
    )
    def test_values_inside_int32_pass(
        self, checks: ModuleType, values: npt.NDArray[np.integer]
    ) -> None:
        checks.check_int32(values, "face_codes")

    @pytest.mark.parametrize("name", ["face_codes", "triangle_codes"])
    @pytest.mark.parametrize(
        "values",
        [
            np.array([311, INT32_MAX + 1], dtype=np.int64),
            np.array([INT32_MIN - 1, 311], dtype=np.int64),
            np.array([INT32_MAX + 1], dtype=np.uint64),
        ],
        ids=["one-past-the-top", "one-past-the-bottom", "unsigned-one-past-the-top"],
    )
    def test_a_value_outside_int32_is_refused_naming_the_array(
        self, checks: ModuleType, values: npt.NDArray[np.integer], name: str
    ) -> None:
        with pytest.raises(ValueError) as refusal:
            checks.check_int32(values, name)
        assert str(refusal.value) == f"{name} must fit int32"
