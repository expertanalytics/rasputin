"""`cli_kit.destination`: the output-path boundary every writing command shares.

Called directly, with no command and no compiled core, so each refusal is
pinned at the function a plugin also calls. Every refusal asserts what it says,
and that it names the option it was given as ``hint``.

Through `rasputin mesh` only the "outside" and the symlink refusals are reached
(`test_cli_mesh_plain_output.py`, `test_cli_mesh_stats.py`); the rest of the
boundary is pinned here, with the missing output folders it makes
(`docs/increments/57-vector-sources.md`, section 4.6, test T8) and `palette
--out`, which goes through it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import typer

from cli_helpers import USAGE, invoke
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

    @pytest.mark.parametrize("hint", ["--out", "--out-edges"])
    def test_the_refusal_names_the_option_it_was_given(self, hint: str, tmp_path: Path) -> None:
        (tmp_path / "file").write_text("in the way", encoding="utf-8")
        error = refusal(tmp_path / "file" / "sub" / "x.ply", None, hint)
        assert error.param_hint == hint


class TestOutputFolders:
    """Missing folders of the output path are made, after every check has
    passed and only under the permitted parent; a refused path makes nothing."""

    def test_a_missing_chain_is_made_with_no_permitted_parent(self, tmp_path: Path) -> None:
        given = tmp_path / "a" / "b" / "c" / "x.ply"
        assert destination(given, None, "mesh") == given.resolve()
        assert (tmp_path / "a" / "b" / "c").is_dir()
        assert not given.exists()  # the folders only: the writer makes the file

    def test_a_missing_chain_is_made_under_the_permitted_parent(self, permitted: Path) -> None:
        given = permitted / "run1" / "deeper" / "x.ply"
        assert destination(given, permitted, "mesh") == given.resolve()
        assert (permitted / "run1" / "deeper").is_dir()

    def test_an_existing_folder_is_left_as_it_is(self, permitted: Path) -> None:
        (permitted / "kept.txt").write_text("kept", encoding="utf-8")
        assert (
            destination(permitted / "x.ply", permitted, "mesh") == (permitted / "x.ply").resolve()
        )
        assert sorted(p.name for p in permitted.iterdir()) == ["kept.txt"]

    def test_a_path_outside_the_permitted_parent_makes_nothing(
        self, permitted: Path, tmp_path: Path
    ) -> None:
        before = sorted(tmp_path.rglob("*"))
        error = refusal(permitted / ".." / "escaped" / "deeper" / "x.ply", permitted)
        assert "outside" in error.message
        assert sorted(tmp_path.rglob("*")) == before

    def test_a_missing_permitted_parent_is_refused_naming_it(self, tmp_path: Path) -> None:
        missing = tmp_path / "no-such-parent"
        before = sorted(tmp_path.rglob("*"))
        error = refusal(missing / "sub" / "x.ply", missing)
        assert "--out-parent" in f"{error.param_hint} {error.message}", error
        assert str(missing) in error.message, error.message
        assert sorted(tmp_path.rglob("*")) == before

    def test_a_file_in_the_way_is_refused_naming_the_folder(self, tmp_path: Path) -> None:
        blocker = tmp_path / "blocker"
        blocker.write_text("a file, not a folder", encoding="utf-8")
        error = refusal(blocker / "sub" / "x.ply", None)
        assert str(blocker) in error.message, error.message
        assert error.param_hint == "--out"
        assert blocker.read_text(encoding="utf-8") == "a file, not a folder"

    def test_a_symlinked_out_is_still_refused_before_any_folder_is_made(
        self, permitted: Path
    ) -> None:
        link = permitted / "link.ply"
        link.symlink_to(permitted / "nowhere" / "real.ply")
        error = refusal(link, permitted)
        assert "symlink" in error.message
        assert not (permitted / "nowhere").exists()


def test_palette_writes_into_a_folder_it_makes(tmp_path: Path) -> None:
    out = tmp_path / "absent" / "deeper" / "p.json"
    code, output = invoke("palette", "corine", "--out", str(out))
    assert code == 0, output
    preset = json.loads(out.read_text(encoding="utf-8"))
    assert preset, preset


def test_palette_refuses_a_symlinked_out(tmp_path: Path) -> None:
    """Through the shared boundary: a link is never followed."""
    real = tmp_path / "real.json"
    real.write_text("untouched", encoding="utf-8")
    link = tmp_path / "link.json"
    link.symlink_to(real)
    code, output = invoke("palette", "corine", "--out", str(link))
    assert code == USAGE, output
    assert "symlink" in output, output
    assert real.read_text(encoding="utf-8") == "untouched"
