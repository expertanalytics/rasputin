r"""Readers of what `rasputin mesh` writes for a reader: increment 25's fields.

`docs/increments/25-plain-output.md`. Test support only; nothing in
`src_python` uses it. One named field from the mesh file (`file_field`), one
`--stats` row of the Inputs and Result sections by its name (`stats_row`),
and a `.ply` header's comments as fields (`ply_fields`). The design (D6) keeps
`file_field` and `stats_row` importable from `test_cli_mesh_refine.py`, which
re-exports them; they live here so `test_cli_mesh_dem.py`, which that module
imports, can use them without a cycle.

The `--stats` layout is the design's: a section `## Inputs` and a section
`## Result`, each one table whose rows are `| wording | value | name |`
(`stats_rows`' tuple, D1, printed in that order; D3 rule 1: the name is the
third column). A name may be written in backticks. A `|` inside a cell is
escaped `\|` and unescaped here (D4).
"""

from __future__ import annotations

import re

from vtkread import VtkFile


def file_field(vtk: VtkFile, name: str) -> str:
    """One named field of the mesh file (increment 25, D2), exactly as written."""
    assert name in vtk.field_data, f"no field {name!r}; the file has {sorted(vtk.field_data)}"
    (text,) = vtk.field_data[name].values
    assert isinstance(text, str)
    return text


_CELL = re.compile(r"(?<!\\)\|")


def _cells(line: str) -> list[str]:
    """A Markdown table row's cells, split on unescaped ``|`` and unescaped."""
    inner = line.strip()[1:-1]
    return [c.strip().replace("\\|", "|") for c in _CELL.split(inner)]


#: ``dem_seams`` when every overlap agrees (increment 25, "The fields, for
#: Ola"); a disagreeing pair keeps ``mosaic.Seam.entry``'s text.
SEAMS_AGREE = "none: the tiles agree where they overlap"

#: Comments the PLY writer adds itself, describing the file's own arrays
#: (D2, "Fields that stay"): not run results, so not compared with ``.vtk``.
PLY_WRITER_COMMENTS = ("feature_bit", "feature_vocabulary")


def ply_fields(comments: tuple[str, ...]) -> dict[str, str]:
    """A ``.ply`` header's ``name value`` comments as a dict, the writer's own
    excluded. A name given twice is a failure: the file has one value per name."""
    fields: dict[str, str] = {}
    for comment in comments:
        name, _, value = comment.partition(" ")
        if name in PLY_WRITER_COMMENTS:
            continue
        assert name not in fields, f"comment {name!r} twice in {comments}"
        fields[name] = value
    return fields


def stats_rows_named(report: str) -> list[tuple[str, str]]:
    """Every (name, value) row of the ``--stats`` Inputs and Result sections,
    in report order (D4). A row is ``| wording | value | name |``; the header
    row of each table is skipped, backticks around the name are dropped, and
    a value's ``\\|`` is unescaped."""
    rows: list[tuple[str, str]] = []
    for heading in ("Inputs", "Result"):
        match = re.search(rf"^## {heading}\n(.*?)(?=^## |\Z)", report, re.M | re.S)
        if match is None:
            continue
        lines = [ln for ln in match.group(1).splitlines() if ln.startswith("|")]
        for i, line in enumerate(lines):
            if line.startswith("|---") or (i + 1 < len(lines) and lines[i + 1].startswith("|---")):
                continue
            cells = _cells(line)
            assert len(cells) == 3, f"an Inputs or Result row has three cells: {line!r}"
            rows.append((cells[2].strip("`"), cells[1]))
    return rows


def stats_row(report: str, name: str) -> str:
    """The value ``--stats`` prints for the field ``name`` (its third column)."""
    found = [value for key, value in stats_rows_named(report) if key == name]
    assert found, f"no --stats row {name!r} in\n{report}"
    assert len(found) == 1, f"--stats has {len(found)} rows named {name!r}"
    return found[0]


def stats_names(report: str) -> list[str]:
    return [name for name, _ in stats_rows_named(report)]


def start_stride(report: str) -> int:
    """The stride of a stride start, from ``start_mesh`` (``every 40th DEM node``)."""
    text = stats_row(report, "start_mesh")
    match = re.fullmatch(r"every (?:(\d+)(?:st|nd|rd|th) )?DEM node", text)
    assert match is not None, f"start_mesh is not a stride start: {text!r}"
    return int(match.group(1) or 1)


def sizes_row(report: str, item: str) -> str:
    """The count ``--stats`` prints in its Sizes table for ``item``."""
    match = re.search(r"^## Sizes\n(.*?)(?=^## |\Z)", report, re.M | re.S)
    assert match is not None, f"no Sizes section in\n{report}"
    for line in match.group(1).splitlines():
        if line.startswith("|") and not line.startswith("|---"):
            cells = _cells(line)
            if cells[0] == item:
                return cells[1]
    raise AssertionError(f"no Sizes row {item!r} in\n{report}")
