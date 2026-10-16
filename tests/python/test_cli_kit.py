"""`cli_kit.destination`: the output-path boundary every writing command shares.

Called directly, with no command and no compiled core, so each refusal is
pinned at the function a plugin also calls. Every refusal asserts what it says,
and that it names the option it was given as ``hint``.

Through `rasputin mesh` only the "outside" and the symlink refusals are reached
(`test_cli_mesh_plain_output.py`, `test_cli_mesh_stats.py`); the rest of the
boundary is pinned here.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import typer

from tin_engine.cli_kit import destination


@pytest.fixture
def permitted(tmp_path: Path) -> Path:
    """An explicitly permitted parent, one level below the test's own directory."""
    parent = tmp_path / "permitted"
    parent.mkdir()
    return parent


def refusal(out: Path, out_parent: Path | None, hint: str = "--out") -> typer.BadParameter:
    with pytest.raises(typer.BadParameter) as caught:
        destination(out, out_parent, "mesh", hint)
    return caught.value


class TestTheBoundary:
    def test_a_path_inside_the_permitted_parent_is_returned_resolved(self, permitted: Path) -> None:
        given = permitted / "sub" / ".." / "inside.ply"
        (permitted / "sub").mkdir()
        assert destination(given, permitted, "mesh") == (permitted / "inside.ply").resolve()

    def test_traversal_out_of_the_permitted_parent_is_refused(
        self, permitted: Path, tmp_path: Path
    ) -> None:
        error = refusal(permitted / ".." / "escaped.ply", permitted)
        assert "outside" in error.message
        assert not (tmp_path / "escaped.ply").exists()

    def test_an_absolute_path_elsewhere_is_refused(self, permitted: Path, tmp_path: Path) -> None:
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        error = refusal(elsewhere / "x.ply", permitted)
        assert "outside" in error.message

    def test_a_symlink_pointing_back_inside_is_refused(self, permitted: Path) -> None:
        # Containment alone accepts this link, since it resolves inside the
        # permitted parent: only the symlink check refuses it.
        real = permitted / "real.ply"
        real.write_text("untouched", encoding="utf-8")
        link = permitted / "link.ply"
        link.symlink_to(real)
        error = refusal(link, permitted)
        assert "symlink" in error.message
        assert real.read_text(encoding="utf-8") == "untouched"

    def test_a_missing_directory_is_a_refusal_naming_it(self, tmp_path: Path) -> None:
        missing = tmp_path / "nope"
        error = refusal(missing / "x.ply", None)
        assert str(missing) in error.message
        assert "not an existing directory" in error.message

    @pytest.mark.parametrize("hint", ["--out", "--out-edges"])
    def test_the_refusal_names_the_option_it_was_given(self, hint: str, tmp_path: Path) -> None:
        error = refusal(tmp_path / "nope" / "x.ply", None, hint)
        assert error.param_hint == hint
