"""`cli_kit.destination`: the output-path boundary every writing command shares.

Called directly, with no command and no compiled core, so each behaviour is
pinned at the function a plugin also calls. Every refusal asserts what it says,
that it names the option it was given as ``hint``, and that nothing was made.

`destination(out, hint)` refuses a symlinked final part, resolves, and makes
the missing folders above the file (`docs/increments/57-vector-sources.md`,
section 4.6, test T8a). It fences nothing: the fence against writing
everywhere is the MCP server's output root (`test_mcp_vector_sources.py`).
`palette --out` goes through it.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest
import typer

from cli_helpers import USAGE, invoke
from tin_engine.cli_kit import destination


def refusal(out: Path, hint: str = "--out") -> typer.BadParameter:
    with pytest.raises(typer.BadParameter) as caught:
        destination(out, hint)
    return caught.value


def everything_under(directory: Path) -> list[Path]:
    return sorted(directory.rglob("*"))


def test_the_signature_is_out_and_hint_only() -> None:
    """`out_parent` and the unused `name` are gone; a plugin calls `destination(out, hint)`."""
    parameters = inspect.signature(destination).parameters
    assert list(parameters) == ["out", "hint"]
    assert parameters["hint"].default == "--out"


class TestTheBoundary:
    def test_a_dotdot_path_is_returned_resolved(self, tmp_path: Path) -> None:
        (tmp_path / "sub").mkdir()
        given = tmp_path / "sub" / ".." / "inside.ply"
        assert destination(given) == (tmp_path / "inside.ply").resolve()

    def test_a_symlink_is_refused_and_its_target_untouched(self, tmp_path: Path) -> None:
        real = tmp_path / "real.ply"
        real.write_text("untouched", encoding="utf-8")
        link = tmp_path / "link.ply"
        link.symlink_to(real)
        before = everything_under(tmp_path)
        error = refusal(link)
        assert "symlink" in error.message
        assert error.param_hint == "--out"
        assert real.read_text(encoding="utf-8") == "untouched"
        assert everything_under(tmp_path) == before

    @pytest.mark.parametrize("hint", ["--out", "--out-edges"])
    def test_the_refusal_names_the_option_it_was_given(self, hint: str, tmp_path: Path) -> None:
        (tmp_path / "file").write_text("in the way", encoding="utf-8")
        error = refusal(tmp_path / "file" / "sub" / "x.ply", hint)
        assert error.param_hint == hint


class TestOutputFolders:
    """Missing folders above the output file are made; a refused path makes nothing."""

    def test_a_missing_chain_is_made_and_the_file_is_not(self, tmp_path: Path) -> None:
        given = tmp_path / "a" / "b" / "c" / "x.ply"
        assert destination(given) == given.resolve()
        assert (tmp_path / "a" / "b" / "c").is_dir()
        assert not given.exists()  # the folders only: the writer makes the file

    def test_an_existing_folder_is_left_as_it_is(self, tmp_path: Path) -> None:
        folder = tmp_path / "kept"
        folder.mkdir()
        (folder / "kept.txt").write_text("kept", encoding="utf-8")
        assert destination(folder / "x.ply") == (folder / "x.ply").resolve()
        assert sorted(p.name for p in folder.iterdir()) == ["kept.txt"]
        assert (folder / "kept.txt").read_text(encoding="utf-8") == "kept"

    def test_dotdot_makes_the_resolved_folder_not_the_spelled_one(self, tmp_path: Path) -> None:
        given = tmp_path / "a" / ".." / "b" / "x.ply"
        assert destination(given) == (tmp_path / "b" / "x.ply").resolve()
        assert (tmp_path / "b").is_dir()
        assert not (tmp_path / "a").exists()

    def test_a_file_in_the_way_is_refused_naming_the_folder(self, tmp_path: Path) -> None:
        blocker = tmp_path / "blocker"
        blocker.write_text("a file, not a folder", encoding="utf-8")
        before = everything_under(tmp_path)
        error = refusal(blocker / "sub" / "x.ply")
        assert str(blocker.resolve()) in error.message, error.message
        assert error.param_hint == "--out"
        assert blocker.read_text(encoding="utf-8") == "a file, not a folder"
        assert everything_under(tmp_path) == before

    def test_a_symlinked_out_is_refused_before_any_folder_is_made(self, tmp_path: Path) -> None:
        link = tmp_path / "link.ply"
        link.symlink_to(tmp_path / "nowhere" / "deeper" / "real.ply")
        before = everything_under(tmp_path)
        error = refusal(link)
        assert "symlink" in error.message
        assert not (tmp_path / "nowhere").exists()
        assert everything_under(tmp_path) == before


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
