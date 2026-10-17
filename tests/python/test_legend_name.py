"""One name for a land-cover legend: the old "class map" spellings are gone.

`docs/increments/58-legend.md`, section 6.1 (the guard) and section 6.2's
three pins that the old names are refused, not aliased. Ola ruled out
aliases, so every old spelling is either history or a bug.

The guard lists the tracked files (`git ls-files`, run from the repository
root above this file) and fails naming path and line for every tracked path
and every line that matches `OLD`, case-insensitively. Exactly three
exceptions, each a rule, no allow-list of lines:

1. paths under `HISTORY` and `ROADMAP.md`: the history records keep the
   names of their time;
2. `HISTORY_FILE_NAME` (a history record's file name, cited from code and
   tests) is removed from a path or a line before matching;
3. this file, which has to spell the old names to refuse them.

Files with a NUL byte in them are not text and are not read; their paths
are still matched. Nor are the prose files of `tools/ci_changes.py`
(`is_prose`) read: a pull request that changes only prose skips this suite,
and `conftest.py` fails a session that reads one. Their paths are matched;
their words are left to review.

The pins below spell the old names on purpose, and so live here (rule 3):
the old `mesh` options are refused by Typer as unknown options, a sources
file with the old key is refused whole, and the MCP `mesh` tool takes only
the new argument names.
"""

from __future__ import annotations

import importlib.util
import re
import subprocess
from pathlib import Path
from types import ModuleType

from cli_helpers import USAGE, invoke
from test_catalogue_data import PER_YEAR, USER_ID, changed, refusal, write
from tin_engine import sources
from tin_engine.cli import app
from tin_engine.mcp import schema as mcp_schema

#: The old spellings, matched case-insensitively.
OLD = re.compile(
    r"class[ _-]?maps?|features[-_]map\b|land[-_]cover[-_]map\b|\bcmap\b", re.IGNORECASE
)

#: Rule 1: history records, by path prefix, and the ROADMAP.
HISTORY = ("docs/increments/", "docs/benchmarks/", "docs/retrospectives/")
ROADMAP = "ROADMAP.md"

#: Rule 2: removed from a path or a line before matching.
HISTORY_FILE_NAME = "41-class-maps-as-data"

ROOT = Path(__file__).resolve().parents[2]
SELF = Path(__file__).resolve().relative_to(ROOT).as_posix()


def _ci_changes() -> ModuleType:
    """`tools/ci_changes.py`, which decides what a prose file is."""
    spec = importlib.util.spec_from_file_location(
        "legend_ci_changes", ROOT / "tools" / "ci_changes.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CI_CHANGES = _ci_changes()


def tracked() -> list[str]:
    """The tracked paths, relative to the repository root."""
    listed = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True
    ).stdout.decode("utf-8")
    return [path for path in listed.split("\0") if path]


def counted(path: str) -> bool:
    """Rules 1 and 3: whether `path` is checked at all."""
    return not (path.startswith(HISTORY) or path in (ROADMAP, SELF))


def old_names(path: str, text: str) -> list[str]:
    """`path:line: text` for every line of `text` with an old spelling;
    `path: (file name)` first when the path has one."""
    found = [f"{path}: (file name)"] if OLD.search(path.replace(HISTORY_FILE_NAME, "")) else []
    for number, line in enumerate(text.splitlines(), start=1):
        if OLD.search(line.replace(HISTORY_FILE_NAME, "")):
            found.append(f"{path}:{number}: {line.strip()}")
    return found


def read_text(path: str) -> str:
    """The file's text, or "" for a prose file, a file that is not text, or
    one that is not there (a tracked path deleted in the working tree)."""
    if CI_CHANGES.is_prose(path) or not (ROOT / path).is_file():
        return ""
    data = (ROOT / path).read_bytes()
    return "" if b"\0" in data else data.decode("utf-8", errors="replace")


def test_no_old_spelling_outside_the_history_records() -> None:
    paths = tracked()
    assert paths, f"git ls-files listed no file under {ROOT}"
    found = [hit for path in paths if counted(path) for hit in old_names(path, read_text(path))]
    assert found == [], f"{len(found)} old spellings of 'legend':\n" + "\n".join(found)


class TestTheRules:
    """Each rule of the guard, on a made-up path and line."""

    def test_a_line_and_a_path_are_both_matched(self) -> None:
        assert old_names("src_python/x.py", "a = 1\nload_class_map(x)\n") == [
            "src_python/x.py:2: load_class_map(x)"
        ]
        assert old_names("data/class_maps/a.toml", "") == ["data/class_maps/a.toml: (file name)"]

    def test_every_old_spelling_is_matched(self) -> None:
        for line in ("Class map", "class-map", "CLASS_MAPS", "ClassMap", "--features-map",
                     "features_map", "--land-cover-map", "land_cover_map", "cmap = 1"):  # fmt: skip
            assert old_names("x.py", line), line

    def test_the_new_names_and_their_neighbours_are_not_matched(self) -> None:
        for line in ("legend", "--features-legend", "land_cover_legend", "colormap",
 "features_mapper"):  # fmt: skip
            assert old_names("x.py", line) == [], line

    def test_the_history_records_and_this_file_are_not_counted(self) -> None:
        for path in ("docs/increments/58-legend.md", "docs/benchmarks/x/stats.txt",
                     "docs/retrospectives/x.md", "ROADMAP.md", SELF):  # fmt: skip
            assert not counted(path), path
        assert counted("docs/GLOSSARY.md") and counted("README.md")

    def test_a_prose_file_is_not_read_but_a_code_file_is(self) -> None:
        assert read_text("docs/GLOSSARY.md") == ""
        assert "def old_names" in read_text(SELF)

    def test_the_history_file_name_is_removed_before_matching(self) -> None:
        line = "see `docs/increments/41-class-maps-as-data.md`, section 7"
        assert old_names("x.py", line) == []
        assert old_names("x.py", line + "; a class map") != []


# ------------------------------------------------- no alias for the old names


def test_the_old_mesh_options_are_refused() -> None:
    for option, value in (("--features-map", "corine"), ("--land-cover-map", "mapbiomas-c11")):
        code, output = invoke("mesh", option, value)
        assert code == USAGE, output
        assert f"No such option: {option}" in output, output


def test_a_sources_file_with_the_old_key_is_refused(tmp_path: Path) -> None:
    user = write(
        tmp_path / "my-sources.toml", {f"source.{USER_ID}": changed(PER_YEAR, class_map="corine")}
    )
    refusal(sources, user, "class_map")


def test_the_mcp_mesh_tool_takes_only_the_new_names() -> None:
    (mesh,) = (d for d in mcp_schema.derive(app) if d.tool_name == "mesh")
    properties = set(mesh.tool().input_schema["properties"])
    assert {"features_legend", "land_cover_legend"} <= properties, sorted(properties)
    assert properties.isdisjoint({"features_map", "land_cover_map"}), sorted(properties)
