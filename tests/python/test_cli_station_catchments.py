"""`rasputin station-catchments`: a catchment per station, a table, a summary (increment 29, PR 4).

`docs/increments/29-nve-reference-catchments.md`, "The batch" (the command,
its options and what it writes), "Ola's rulings" (the stderr line of copies
dropped and `river_copies_dropped`) and "The red suites", PR 4's
`test_cli_station_catchments.py`. The inputs are `batch_fixtures`' five
stations, rivers and references, written as GeoJSON, and the terrain as four
GeoTIFF tiles; no network.

Pinned here beyond the design:

- `DIR/<station>.geojson` is written for every station with a catchment,
  none for a refused one; its properties hold `station`, `name`, `series`,
  the placement and the sensitivity (`placed_on`, `causes`, `swing`,
  `reach_fork`) beside 22's (`nodes`, `fine_area_m2`, ...).
- `DIR/results.csv` has a header row and one row per station in file order;
  its columns are `StationResult`'s fields, with `station_class` headed
  `class`; a tuple (`causes`, `grid_tiles`) is written joined by `;`, and
  None as an empty cell.
- `DIR/summary.json` is the batch's summary (`reference.Summary`) plus
  `river_copies_dropped`.
- stderr has one line per station holding its number, its name and its
  class (the word `refused` for a refusal), and the summary after them.
- A stations or rivers file without a `crs` member is refused, naming the
  option and the missing CRS, and nothing is written.
- A reference file whose CRS is not the river file's is refused naming
  `--reference` and both CRSs by their EPSG codes, and nothing is written.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

import batch_fixtures as bf
from test_cli_mesh import plain
from test_cli_mesh_mosaic import write_tiles
from tin_engine.cli import app
from tin_engine.domain import read_domain

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

PLACED = [s for s in bf.FIVE if s is not bf.LANGT]
CLASSES = {
    bf.TREFF.station: "match",
    bf.BOM.station: "miss",
    bf.SAMLOP.station: "uncertain",
    bf.SLUTT.station: "uncertain",
    bf.LANGT.station: "refused",
}


def invoke(*args: str) -> tuple[int, str]:
    result = runner.invoke(app, list(args))
    return result.exit_code, result.output


@pytest.fixture(scope="module")
def data(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("inputs")
    write_tiles(root / "dem", bf.basin_tiles())
    return {
        "dem": root / "dem",
        "stations": bf.write_stations(root / "stations.geojson"),
        "rivers": bf.write_rivers(root / "rivers.geojson"),
        "reference": bf.write_references(root / "reference.geojson"),
    }


def args(data: dict[str, Path], out: Path, *, reference: bool = True) -> list[str]:
    a = ["station-catchments", "--dem", str(data["dem"]), "--stations", str(data["stations"]),
         "--rivers", str(data["rivers"]), "--out-dir", str(out)]  # fmt: skip
    if reference:
        a += ["--reference", str(data["reference"])]
    return a


@pytest.fixture(scope="module")
def full(data: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, str]:
    out = tmp_path_factory.mktemp("full") / "out"
    code, output = invoke(*args(data, out))
    assert code == 0, output
    return out, output


def table(out: Path) -> list[dict[str, str]]:
    with (out / "results.csv").open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
# What it writes
# ---------------------------------------------------------------------------


def test_a_catchment_file_per_placed_station_and_none_for_the_refused(
    full: tuple[Path, str],
) -> None:
    out, _ = full
    names = sorted(p.name for p in out.glob("*.geojson"))
    assert names == sorted(f"{s.station}.geojson" for s in PLACED)
    assert {p.name for p in out.iterdir()} == {*names, "results.csv", "summary.json"}


def test_the_catchment_file_is_a_domain_with_the_station_and_its_gauge(
    full: tuple[Path, str],
) -> None:
    out, _ = full
    path = out / f"{bf.TREFF.station}.geojson"
    assert read_domain(path).polygon.is_valid
    props = json.loads(path.read_text(encoding="utf-8"))["features"][0]["properties"]
    assert props["station"] == bf.TREFF.station
    assert props["name"] == bf.TREFF.name
    assert props["series"] == ["1001.0"]
    assert props["placed_on"] == "number"
    assert props["causes"] == []
    assert props["reach_fork"] is False
    assert props["nodes"] == int(bf.flood_mask(bf.TREFF).sum())
    assert props["swing"] == pytest.approx(bf.swing_on_raw(bf.TREFF), rel=1e-9)


def test_mesh_reads_a_station_catchment_as_its_domain(
    full: tuple[Path, str], data: dict[str, Path], tmp_path: Path
) -> None:
    out, _ = full
    vtk = tmp_path / "m.vtk"
    code, output = invoke(
        "mesh", "--dem", str(data["dem"]), "--domain", str(out / f"{bf.TREFF.station}.geojson"),
        "--tolerance", "5", "--out", str(vtk),
    )  # fmt: skip
    assert code == 0, output
    assert vtk.is_file() and vtk.stat().st_size > 0


def test_the_table_has_a_row_per_station_in_file_order(full: tuple[Path, str]) -> None:
    rows = table(full[0])
    assert [r["station"] for r in rows] == [s.station for s in bf.FIVE]
    assert {r["station"]: r["class"] for r in rows} == CLASSES
    by = {r["station"]: r for r in rows}
    assert by[bf.TREFF.station]["match_by"] == "overlap"
    assert by[bf.LANGT.station]["refusal_cause"] == "no_river"
    assert by[bf.LANGT.station]["placed_on"] == ""
    assert by[bf.SLUTT.station]["causes"] == "downstream_unread"
    assert by[bf.TREFF.station]["name"] == bf.TREFF.name


def test_the_summary_json(full: tuple[Path, str]) -> None:
    s = json.loads((full[0] / "summary.json").read_text(encoding="utf-8"))
    assert s["classes"] == {"match": 1, "close": 0, "miss": 1, "uncertain": 2, "refused": 1}
    assert s["river_copies_dropped"] == 1


# ---------------------------------------------------------------------------
# What it says
# ---------------------------------------------------------------------------


def station_line(output: str, station: str) -> str:
    lines = [
        ln for ln in output.splitlines() if re.search(rf"(?<![\d.]){re.escape(station)}\b", ln)
    ]
    assert lines, f"no line names {station}:\n{output}"
    return lines[0]


def test_one_stderr_line_per_station_with_its_name_and_class(full: tuple[Path, str]) -> None:
    _, output = full
    for spec in bf.FIVE:
        line = station_line(output, spec.station)
        assert spec.name in line, line
        assert CLASSES[spec.station] in line, line


def test_the_summary_comes_after_the_station_lines(full: tuple[Path, str]) -> None:
    _, output = full
    last = station_line(output, bf.LANGT.station)
    after = output[output.index(last) + len(last) :]
    assert re.search(r"\buncertain\b", after), output
    assert re.search(r"\bmatch\b", after), output


def test_the_copies_dropped_line(full: tuple[Path, str]) -> None:
    _, output = full
    assert re.search(r"\brivers: 4 segments read, 1 exact copy dropped\b", output), output


# ---------------------------------------------------------------------------
# Options and refusals
# ---------------------------------------------------------------------------


def test_only_runs_the_named_stations(data: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "out"
    code, output = invoke(*args(data, out), "--only", bf.SLUTT.station, "--only", bf.TREFF.station)
    assert code == 0, output
    assert [r["station"] for r in table(out)] == [bf.TREFF.station, bf.SLUTT.station]
    assert sorted(p.name for p in out.glob("*.geojson")) == sorted(
        [f"{bf.TREFF.station}.geojson", f"{bf.SLUTT.station}.geojson"]
    )


def test_without_reference_catchments_and_no_scored_classes(
    data: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "out"
    code, output = invoke(*args(data, out, reference=False))
    assert code == 0, output
    rows = {r["station"]: r for r in table(out)}
    assert rows[bf.TREFF.station]["class"] == ""
    assert rows[bf.BOM.station]["class"] == ""
    assert rows[bf.SAMLOP.station]["class"] == "uncertain"
    assert rows[bf.LANGT.station]["class"] == "refused"
    assert {r["class"] for r in rows.values()} <= {"", "uncertain", "refused"}
    assert len(list(out.glob("*.geojson"))) == len(PLACED)


@pytest.mark.parametrize("which", ["stations", "rivers"])
def test_a_file_without_crs_is_refused_and_nothing_is_written(
    data: dict[str, Path], tmp_path: Path, which: str
) -> None:
    bare = dict(data)
    if which == "stations":
        bare["stations"] = bf.write_stations(tmp_path / "s.geojson", crs=None)
    else:
        bare["rivers"] = bf.write_rivers(tmp_path / "r.geojson", crs=None)
    out = tmp_path / "out"
    code, output = invoke(*args(bare, out))
    words = plain(output)
    assert code != 0, words
    assert "No such command" not in words and "No such option" not in words, words
    assert f"--{which}" in words, words
    assert re.search(r"(?i)\bcrs\b", words), words
    assert not (out / "results.csv").exists()
    assert not out.exists() or not list(out.glob("*.geojson"))


def test_a_reference_in_another_crs_is_refused_and_nothing_is_written(
    data: dict[str, Path], tmp_path: Path
) -> None:
    """The rivers are in EPSG:25833 and `run_batch` takes the references in
    that CRS, so a reference file in EPSG:32633 is refused, naming the option
    and both CRSs; no polygon is reprojected ("PR 4's red step")."""
    other = dict(data)
    other["reference"] = bf.write_references(tmp_path / "ref.geojson", crs="EPSG:32633")
    out = tmp_path / "out"
    code, output = invoke(*args(other, out))
    words = plain(output)
    assert code != 0, words
    assert "No such command" not in words and "No such option" not in words, words
    assert "--reference" in words, words
    assert "32633" in words and "25833" in words, words
    assert not (out / "results.csv").exists()
    assert not (out / "summary.json").exists()
    assert not out.exists() or not list(out.glob("*.geojson"))


def test_the_dem_option_is_required(data: dict[str, Path], tmp_path: Path) -> None:
    a: list[Any] = args(data, tmp_path / "out")
    i = a.index("--dem")
    del a[i : i + 2]
    code, output = invoke(*a)
    assert code != 0
    assert "No such command" not in plain(output), plain(output)
