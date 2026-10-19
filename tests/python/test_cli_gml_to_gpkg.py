"""`rasputin gml-to-gpkg SOURCE --out OUT [--crs TEXT]`: the conversion as a command.

`docs/increments/73-gml-vector-source.md`, sections 4.3 and 9, test 9: one
line on stderr in plain words, ``OUT`` on stdout, ``OUT`` written as
``OUT.part`` and renamed only when finished, and every refusal a usage error
(exit 2) that names SOURCE and leaves no file.

Each test first requires the command to exist, so a missing command fails as
that and never passes as an exit-2 refusal.
"""

from __future__ import annotations

import sqlite3
from contextlib import closing
from pathlib import Path

import pytest
from click.testing import Result

from cli_helpers import COMMANDS, USAGE, plain, runner
from gml32_fixtures import AR50_SAMPLE, document, feature, polygon
from tin_engine.cli import app
from tin_engine.io.geopackage import layer_info, layer_summary
from tin_engine.io.repository import open_geopackage

COMMAND = "gml-to-gpkg"
LINE = (
    "6 features of ArealressursFlate in EPSG:25832, 0 without a surface; columns gml_id,"
    " informasjon, oppdateringsdato, arealtype, skogbonitet, treslag, jordbruk,"
    " vegetasjonsdekke, kartstandard"
)


def run(*args: str) -> Result:
    assert COMMAND in COMMANDS, f"{COMMAND!r} is not a subcommand; they are {sorted(COMMANDS)}"
    return runner.invoke(app, [COMMAND, *args])


def features_in(path: Path) -> int:
    with closing(open_geopackage(path)) as conn:
        return layer_summary(conn, layer_info(conn, None)).features


def test_converts_fixture_a(tmp_path: Path) -> None:
    out = tmp_path / "x.gpkg"
    result = run(str(AR50_SAMPLE), "--out", str(out))
    assert result.exit_code == 0, result.output
    assert result.stdout.strip() == str(out)
    stderr = plain(result.stderr)
    assert "x.gpkg:" in stderr and LINE in stderr, stderr
    assert features_in(out) == 6
    assert not (tmp_path / "x.gpkg.part").exists()


def test_an_existing_output_is_replaced(tmp_path: Path) -> None:
    out = tmp_path / "x.gpkg"
    out.write_bytes(b"not a geopackage")
    result = run(str(AR50_SAMPLE), "--out", str(out))
    assert result.exit_code == 0, result.output
    assert features_in(out) == 6


def test_a_stale_part_file_is_removed_first(tmp_path: Path) -> None:
    out = tmp_path / "x.gpkg"
    sqlite3.connect(tmp_path / "x.gpkg.part").execute("CREATE TABLE stale (a)").connection.close()
    result = run(str(AR50_SAMPLE), "--out", str(out))
    assert result.exit_code == 0, result.output
    with closing(open_geopackage(out)) as conn:
        names = {n for (n,) in conn.execute("SELECT name FROM sqlite_master")}
    assert "stale" not in names
    assert not (tmp_path / "x.gpkg.part").exists()


def test_the_crs_option_names_a_crs_the_file_lacks(tmp_path: Path) -> None:
    source = tmp_path / "bare.gml"
    source.write_bytes(document(feature("f1", polygon(srs_name=None))).getvalue())
    out = tmp_path / "x.gpkg"
    result = run(str(source), "--out", str(out), "--crs", "EPSG:25833")
    assert result.exit_code == 0, result.output
    with closing(open_geopackage(out)) as conn:
        assert layer_info(conn, None).srs_id == 25833
    assert "in EPSG:25833" in plain(result.stderr)


@pytest.mark.parametrize("name", ["x.sqlite", "x.gpkg.part", "x"])
def test_an_output_not_ending_in_gpkg_is_a_usage_error(tmp_path: Path, name: str) -> None:
    out = tmp_path / name
    result = run(str(AR50_SAMPLE), "--out", str(out))
    assert result.exit_code == USAGE, result.output
    assert ".gpkg" in plain(result.output)
    assert "No such" not in result.output, "refused for the wrong reason"
    assert sorted(p.name for p in tmp_path.iterdir()) == []


def test_a_refused_document_leaves_no_file(tmp_path: Path) -> None:
    source = tmp_path / "curve.gml"
    line = (
        '<gml:LineString srsName="EPSG:25832"><gml:posList>0 0 1 1</gml:posList></gml:LineString>'
    )
    source.write_bytes(document(feature("f.curve", line)).getvalue())
    out = tmp_path / "x.gpkg"
    result = run(str(source), "--out", str(out))
    assert result.exit_code == USAGE, result.output
    text = plain(result.output)
    assert "curve.gml" in text and "f.curve" in text, text
    assert "Traceback" not in result.output
    assert not out.exists()
    assert not (tmp_path / "x.gpkg.part").exists()


def test_a_refusal_after_a_written_chunk_leaves_no_file(tmp_path: Path) -> None:
    """The bad feature comes after more than a chunk has been inserted."""
    good = [feature(f"f{k}", polygon()) for k in range(1_500)]
    source = tmp_path / "late.gml"
    source.write_bytes(document(*good, feature("f.late", polygon("0 0 1 0 0 0"))).getvalue())
    out = tmp_path / "x.gpkg"
    result = run(str(source), "--out", str(out))
    assert result.exit_code == USAGE, result.output
    assert "f.late" in plain(result.output)
    assert not out.exists()
    assert not (tmp_path / "x.gpkg.part").exists()


def test_a_missing_source_is_a_usage_error(tmp_path: Path) -> None:
    out = tmp_path / "x.gpkg"
    result = run(str(tmp_path / "absent.gml"), "--out", str(out))
    assert result.exit_code == USAGE, result.output
    assert "absent.gml" in plain(result.output)
    assert not out.exists()


def test_help_says_what_it_does() -> None:
    result = run("--help")
    assert result.exit_code == 0, result.output
    assert "Convert a GML 3.2 file of surfaces to a GeoPackage" in plain(result.output)
