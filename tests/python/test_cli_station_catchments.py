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
- ("PR 4's green step", changes (a) and (c).) An `--only` naming a station
  not in the stations file is refused naming `--only` and the number, before
  `--out-dir` is created; a run stopped by an exception in the batch keeps
  `results.csv`'s header and every finished row, and writes no
  `summary.json`.
- ("PR 4's code review, round 1", changes (d) to (f).) A river file whose
  CRS is not the DEM's is refused naming `--rivers` and both EPSG codes,
  with `delineate` never called and `--out-dir` not created. A write that
  fails, making `--out-dir` itself included (an existing regular file, or a
  parent that cannot be written), is a refusal (no exception escapes to
  the runner) naming `--out-dir` and the path, as given or resolved (macOS
  resolves `/var` to `/private/var`), without `--dem` or "cannot read". A
  station whose `name` is missing or empty gives a stderr line starting `<station>: <class>`.
"""

from __future__ import annotations

import csv
import json
import os
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

import batch_fixtures as bf
import tin_engine.catchment_batch as catchment_batch
from cli_driver import plain, runner
from nve_fixtures import collection, write
from test_cli_mesh_mosaic import write_tiles
from tin_engine.cli import app
from tin_engine.domain import read_domain

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


def test_an_unknown_only_is_refused_naming_the_option_before_out_dir_is_made(
    data: dict[str, Path], tmp_path: Path
) -> None:
    """Change (a) of "PR 4's green step": the command checks `--only` against the
    stations read, before `--out-dir` is created, as a refusal of `--only`
    (not a plain `Error:` line after the directory exists)."""
    out = tmp_path / "out"
    code, output = invoke(*args(data, out), "--only", bf.TREFF.station, "--only", "9.9.9")
    words = plain(output)
    assert code != 0, words
    assert "No such command" not in words and "No such option" not in words, words
    assert "--only" in words, words
    assert "9.9.9" in words, words
    assert not out.exists()


# ---------------------------------------------------------------------------
# A run stopped by a bug
# ---------------------------------------------------------------------------


def test_a_run_stopped_by_a_bug_keeps_the_finished_rows(
    data: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Change (c) of "PR 4's green step": the directory sink writes `results.csv`'s
    header when it is made and flushes each row as it arrives, so an exception
    at the second station leaves the first station's row; `summary.json` is
    written only at the end. `delineate` is replaced by a stub that runs the
    real one for the first station and raises `RuntimeError` at the second
    (a bug, not a refusal the batch classes)."""
    real = catchment_batch.delineate
    calls: list[int] = []

    def first_then_bug(*a: Any, **k: Any) -> Any:
        calls.append(1)
        if len(calls) > 1:
            raise RuntimeError("a planted bug at the second station")
        return real(*a, **k)

    monkeypatch.setattr(catchment_batch, "delineate", first_then_bug)
    out = tmp_path / "out"
    result = runner.invoke(app, args(data, out))
    assert result.exit_code != 0, result.output
    assert isinstance(result.exception, RuntimeError), result.output
    assert len(calls) == 2
    with (out / "results.csv").open(encoding="utf-8", newline="") as f:
        lines = list(csv.reader(f))
    assert lines[0][:2] == ["station", "name"]
    assert "class" in lines[0]
    assert [ln[0] for ln in lines[1:]] == [bf.TREFF.station]
    assert lines[1][lines[0].index("class")] == CLASSES[bf.TREFF.station]
    assert not (out / "summary.json").exists()


# ---------------------------------------------------------------------------
# PR 4's code review, round 1: changes (d), (e) and (f)
# ---------------------------------------------------------------------------


def test_a_river_file_not_in_the_dems_crs_is_refused_before_any_station_runs(
    data: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Change (d): the rivers and the references in EPSG:32633 over the DEM's
    EPSG:25833. The references are in the river file's CRS, so their own check
    passes and the river file's is the one that fires: a refusal of `--rivers`
    naming both CRSs, before any station is delineated and before `--out-dir`
    exists (before change (d), each station was a `refused` row with cause
    `other`)."""
    other = dict(data)
    other["rivers"] = bf.write_rivers(tmp_path / "r.geojson", crs="EPSG:32633")
    other["reference"] = bf.write_references(tmp_path / "ref.geojson", crs="EPSG:32633")
    calls: list[int] = []
    real = catchment_batch.delineate

    def counted(*a: Any, **k: Any) -> Any:
        calls.append(1)
        return real(*a, **k)

    monkeypatch.setattr(catchment_batch, "delineate", counted)
    out = tmp_path / "out"
    code, output = invoke(*args(other, out))
    words = plain(output)
    assert code != 0, words
    assert "No such command" not in words and "No such option" not in words, words
    assert "--rivers" in words, words
    assert "32633" in words and "25833" in words, words
    assert calls == []
    assert not out.exists()


def squashed(text: str) -> str:
    """`plain` with every space gone: Rich breaks a long path across panel
    lines wherever the terminal's width falls, so a path is looked for with
    the whitespace removed (no path here holds a space)."""
    return "".join(plain(text).split())


def an_existing_file(out: Path) -> Path:
    """`--out-dir` itself is a regular file, so making it fails."""
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("not a directory\n", encoding="utf-8")
    return out


def a_read_only_parent(out: Path) -> Path:
    """`--out-dir`'s parent exists but cannot be written, so making it fails.
    Root ignores permission bits, so the case is skipped when run as root."""
    if not hasattr(os, "geteuid") or os.geteuid() == 0:
        pytest.skip("permission bits are not enforced for this user")
    out.parent.mkdir(parents=True)
    out.parent.chmod(0o555)
    return out


def a_directory_at(name: str) -> Any:
    """A directory made beforehand where the command writes the file `name`."""

    def setup(out: Path) -> Path:
        (out / name).mkdir(parents=True)
        return out / name

    return setup


@pytest.fixture
def writable_again(tmp_path: Path) -> Iterator[None]:
    """Gives every directory under `tmp_path` its write bit back afterwards,
    so pytest can remove what `a_read_only_parent` locked."""
    yield
    for p in [tmp_path, *tmp_path.rglob("*")]:
        if p.is_dir():
            p.chmod(0o755)


@pytest.mark.parametrize(
    "blocked",
    [
        a_directory_at("results.csv"),
        a_directory_at(f"{bf.TREFF.station}.geojson"),
        a_directory_at("summary.json"),
        an_existing_file,
        a_read_only_parent,
    ],
    ids=["results_csv", "catchment_file", "summary_json", "out_dir_is_a_file",
         "out_dir_parent_read_only"],
)  # fmt: skip
def test_a_write_that_fails_names_out_dir_and_the_path(
    data: dict[str, Path], tmp_path: Path, blocked: Any, writable_again: None
) -> None:
    """Change (e): a write the command cannot make. A directory made
    beforehand where the command writes a file (results.csv when the sink is
    made, the first station's catchment file, summary.json at the end), or
    `--out-dir` itself that cannot be made: an existing regular file, or a
    directory in a parent that cannot be written. The refusal names
    `--out-dir` and the path, not `--dem` or "cannot read", and is a refusal,
    not a traceback. Before change (e), the catchment file's failure read as a
    refusal of `--dem` and the other two were tracebacks; before the change
    for `--out-dir` itself, making it raised `FileExistsError` or
    `PermissionError` as a traceback."""
    out = tmp_path / "parent" / "out"
    path = blocked(out)
    result = runner.invoke(app, args(data, out))
    words = plain(result.output)
    assert result.exit_code != 0, words
    assert result.exception is None or isinstance(result.exception, SystemExit), repr(
        result.exception
    )
    assert "Traceback" not in result.output
    assert "--out-dir" in words, words
    assert any(str(p) in squashed(result.output) for p in (path, path.resolve())), words
    assert "--dem" not in words, words
    assert "cannot read" not in words, words


@pytest.mark.parametrize("name", [None, ""], ids=["no_name", "empty_name"])
def test_a_station_with_no_name_prints_no_space_before_the_colon(
    data: dict[str, Path], tmp_path: Path, name: str | None
) -> None:
    """Change (f): the stderr line is `<station>: <class>` when the station
    has no name (no `name` property, or an empty one), not `<station> : ...`."""
    feature = bf.station_feature(bf.TREFF)
    if name is None:
        del feature["properties"]["name"]
    else:
        feature["properties"]["name"] = name
    unnamed = dict(data)
    unnamed["stations"] = write(tmp_path / "s.geojson", collection([feature]))
    out = tmp_path / "out"
    code, output = invoke(*args(unnamed, out))
    assert code == 0, output
    line = station_line(output, bf.TREFF.station)
    assert line.startswith(f"{bf.TREFF.station}: match"), line


# ---------------------------------------------------------------------------
# PR 4, lake gauges: `--lakes` ("Lake gauges", "The row and the summary")
# ---------------------------------------------------------------------------
#
# Pinned here beyond the design: the five columns are headed by their field
# names (`seeded_by`, `lake_rule`, `lake_number`, `lake_name`,
# `lake_distance_m`) right after `reach_fork`; "without `--lakes` ... the
# five new columns are empty" is read as the four `lake_*` columns empty and
# `seeded_by` `river` on every placed row (the design's own definition of
# `seeded_by`), empty on a row `place` refused. Before the change the
# command has no `--lakes` option and the table none of the five columns.

LAKE_HEADER = ["seeded_by", "lake_rule", "lake_number", "lake_name", "lake_distance_m"]


@pytest.fixture(scope="module")
def lakes_full(data: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, str]:
    """The five stations with THE LAKE of `batch_fixtures` round `SAMLOP`."""
    root = tmp_path_factory.mktemp("lakes")
    lakes = bf.write_lakes(root / "lakes.geojson")
    out = root / "out"
    code, output = invoke(*args(data, out), "--lakes", str(lakes))
    assert code == 0, output
    return out, output


def test_lakes_writes_the_five_columns_after_reach_fork(lakes_full: tuple[Path, str]) -> None:
    out, _ = lakes_full
    with (out / "results.csv").open(encoding="utf-8", newline="") as f:
        header = next(csv.reader(f))
    at = header.index("reach_fork")
    assert header[at + 1 : at + 6] == LAKE_HEADER


def test_the_lake_row_and_the_river_rows(lakes_full: tuple[Path, str]) -> None:
    rows = {r["station"]: r for r in table(lakes_full[0])}
    samlop = rows[bf.SAMLOP.station]
    assert [samlop[k] for k in LAKE_HEADER] == [
        "lake", "inside", str(bf.LAKE_NUMBER), bf.LAKE_NAME, "0.0",
    ]  # fmt: skip
    assert samlop["class"] not in ("uncertain", "refused"), samlop["class"]
    assert samlop["causes"] == "" and samlop["swing"] == ""
    for spec in (bf.TREFF, bf.BOM, bf.SLUTT):
        assert [rows[spec.station][k] for k in LAKE_HEADER] == ["river", "", "", "", ""]
    assert [rows[bf.LANGT.station][k] for k in LAKE_HEADER] == [""] * 5


def test_the_lake_rows_stderr_words(lakes_full: tuple[Path, str]) -> None:
    line = station_line(lakes_full[1], bf.SAMLOP.station)
    assert re.search(
        rf"\b(match|close|miss) \(seeded by the lake {re.escape(bf.LAKE_NAME)}\)", line
    ), line
    assert "seeded by" not in station_line(lakes_full[1], bf.TREFF.station)


def test_the_lake_rows_catchment_file_carries_the_five(lakes_full: tuple[Path, str]) -> None:
    path = lakes_full[0] / f"{bf.SAMLOP.station}.geojson"
    props = json.loads(path.read_text(encoding="utf-8"))["features"][0]["properties"]
    assert props["seeded_by"] == "lake" and props["lake_rule"] == "inside"
    assert props["lake_number"] == bf.LAKE_NUMBER and props["lake_name"] == bf.LAKE_NAME
    assert props["lake_distance_m"] == 0.0
    assert props["causes"] == [] and props["swing"] is None


def test_the_summary_json_has_by_seed(lakes_full: tuple[Path, str]) -> None:
    s = json.loads((lakes_full[0] / "summary.json").read_text(encoding="utf-8"))
    assert list(s["by_seed"]) == ["river", "lake"]
    assert s["by_seed"]["lake"]["stations"] == 1
    assert s["by_seed"]["river"]["stations"] == 3  # SAMLOP is the lake's; LANGT neither's


@pytest.mark.parametrize(
    ("number", "name", "words"),
    [(bf.LAKE_NUMBER, None, f"the lake {bf.LAKE_NUMBER}"), (None, None, "its lake")],
    ids=["number_no_name", "neither"],
)
def test_the_stderr_words_without_a_lake_name(
    data: dict[str, Path], tmp_path: Path, number: int | None, name: str | None, words: str
) -> None:
    lakes = bf.write_lakes(
        tmp_path / "lakes.geojson",
        [bf.lake_feature(bf.lake_polygon(), 1, vatnlnr=number, navn=name)],
    )
    out = tmp_path / "out"
    code, output = invoke(*args(data, out), "--only", bf.SAMLOP.station, "--lakes", str(lakes))
    assert code == 0, output
    line = station_line(output, bf.SAMLOP.station)
    assert f"(seeded by {words})" in line, line


def test_a_lakes_file_in_another_crs_is_refused_before_out_dir_is_made(
    data: dict[str, Path], tmp_path: Path
) -> None:
    lakes = bf.write_lakes(tmp_path / "lakes.geojson", crs="EPSG:32633")
    out = tmp_path / "out"
    code, output = invoke(*args(data, out), "--lakes", str(lakes))
    words = plain(output)
    assert code != 0, words
    assert "No such command" not in words and "No such option" not in words, words
    assert "--lakes" in words, words
    assert "32633" in words and "25833" in words, words
    assert not out.exists()


def test_a_lakes_file_without_crs_is_refused_naming_lakes(
    data: dict[str, Path], tmp_path: Path
) -> None:
    lakes = bf.write_lakes(tmp_path / "lakes.geojson", crs=None)
    out = tmp_path / "out"
    code, output = invoke(*args(data, out), "--lakes", str(lakes))
    words = plain(output)
    assert code != 0, words
    assert "No such option" not in words, words
    assert "--lakes" in words and re.search(r"(?i)\bcrs\b", words), words
    assert not out.exists()


def test_without_lakes_every_row_is_a_river_row(full: tuple[Path, str]) -> None:
    out, output = full
    rows = {r["station"]: r for r in table(out)}
    for spec in bf.FIVE:
        seeded = "" if spec is bf.LANGT else "river"
        assert [rows[spec.station][k] for k in LAKE_HEADER] == [seeded, "", "", "", ""]
    assert "seeded by" not in output
