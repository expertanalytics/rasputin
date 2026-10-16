"""PR F's safety net: what `catchment` and `station-catchments` write, as hashes.

`docs/increments/python-audit.md`, section 11 ("The safety net"). Builds the
existing suites' fixtures (`batch_fixtures`, `gauge_fixtures`) in a temporary
directory, runs eight commands through Typer's runner, and prints one line per
output: a SHA-256 of each file the command wrote and of its stderr, with the
exit code. Timings are masked: the `seconds` column of results.csv is
emptied, and every `<number> s` on stderr becomes `# s`; the temporary
directory's path becomes `<tmp>`. Nothing else is changed. Run it before and
after the change, from the repository root, with the interpreter whose
`tin_engine` is the tree to measure (its first line prints that path):

    PYTHONPATH=tests/python .venv/bin/python \
        docs/increments/python-audit-probes/catchment_bytes.py

The two outputs, less the first line, must be identical.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re
import sys
import tempfile
from pathlib import Path

from typer.testing import CliRunner

import batch_fixtures as bf
import gauge_fixtures as gf
import tin_engine
from mosaic_fixtures import quadrants
from nve_fixtures import collection, river, write
from test_cli_mesh_mosaic import write_tiles
from tin_engine.cli import app

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})


def masked(text: str, tmp: Path) -> str:
    return re.sub(r"\d+(\.\d+)? s\b", "# s", text.replace(str(tmp), "<tmp>"))


def table_bytes(path: Path) -> bytes:
    rows = list(csv.reader(path.open(encoding="utf-8", newline="")))
    k = rows[0].index("seconds")
    out = io.StringIO()
    csv.writer(out).writerows([*rows[:1], *([*r[:k], "", *r[k + 1 :]] for r in rows[1:])])
    return out.getvalue().encode()


def run(label: str, tmp: Path, out: Path, *args: str) -> None:
    result = runner.invoke(app, list(args))
    digest = hashlib.sha256(masked(result.output, tmp).encode()).hexdigest()
    print(f"{label} exit={result.exit_code} stderr+stdout {digest}")
    files = sorted(out.rglob("*")) if out.is_dir() else [out] if out.exists() else []
    for path in files:
        data = table_bytes(path) if path.name == "results.csv" else path.read_bytes()
        print(f"{label} {path.relative_to(tmp)} {hashlib.sha256(data).hexdigest()}")


def main() -> None:
    print(f"tin_engine from {Path(tin_engine.__file__).parent}")
    with tempfile.TemporaryDirectory() as name:
        tmp = Path(name)
        write_tiles(tmp / "basins", bf.basin_tiles())
        write_tiles(tmp / "mixed", bf.mixed_tiles())
        valley = gf.tile_of(gf.valley(dam=True))
        write_tiles(tmp / "valley", quadrants(valley, row_cut=120, col_cut=30, overlap=1))
        stations = bf.write_stations(tmp / "stations.geojson")
        grense = bf.write_stations(tmp / "grense.geojson", [bf.TREFF, bf.GRENSE])
        rivers = bf.write_rivers(tmp / "rivers.geojson")
        mixed_rivers = bf.write_rivers(tmp / "mixed_rivers.geojson", mixed=True)
        lake_rivers = write(tmp / "lake_rivers.geojson",
                            collection(bf.lake_line_rivers(), crs=gf.EPSG))  # fmt: skip
        reference = bf.write_references(tmp / "reference.geojson")
        lakes = bf.write_lakes(tmp / "lakes.geojson")
        line_lake = bf.write_lakes(
            tmp / "line_lake.geojson", [bf.lake_feature(bf.lake_line_polygon(), 2)]
        )
        line = list(gf.column_line(gf.CC + 0.4, 20, 230))
        valley_rivers = write(tmp / "valley_rivers.geojson", collection(
            [river(8841, line, elvid="2-11-1", elvenavn="Nea", vassdragsnr="002.A")], crs=gf.EPSG
        ))  # fmt: skip
        batch = ("station-catchments", "--dem", str(tmp / "basins"), "--stations", str(stations))
        run("batch", tmp, tmp / "o1", *batch, "--rivers", str(rivers), "--out-dir",
            str(tmp / "o1"), "--reference", str(reference))  # fmt: skip
        run("batch-lakes", tmp, tmp / "o2", *batch, "--rivers", str(rivers), "--out-dir",
            str(tmp / "o2"), "--lakes", str(lakes))  # fmt: skip
        run("batch-lake-line", tmp, tmp / "o3", *batch, "--rivers", str(lake_rivers),
            "--out-dir", str(tmp / "o3"), "--lakes", str(line_lake))  # fmt: skip
        run("batch-mixed", tmp, tmp / "o4", "station-catchments", "--dem", str(tmp / "mixed"),
            "--stations", str(grense), "--rivers", str(mixed_rivers), "--out-dir",
            str(tmp / "o4"))  # fmt: skip
        x, y = gf.lat(gf.CC + 2, 150)
        seed = ("--seed", repr(x), repr(y), "--seed-crs", gf.EPSG)
        valley_dem = ("catchment", "--dem", str(tmp / "valley"), *seed)
        run("single-river", tmp, tmp / "r.geojson", *valley_dem, "--rivers",
            str(valley_rivers), "--out", str(tmp / "r.geojson"))  # fmt: skip
        run("single-far", tmp, tmp / "f.geojson", *valley_dem, "--rivers", str(valley_rivers),
            "--map-radius", "15", "--out", str(tmp / "f.geojson"))  # fmt: skip
        run("single-outlet", tmp, tmp / "p.geojson", *valley_dem, "--out",
            str(tmp / "p.geojson"))  # fmt: skip
        lx, ly = bf.lake_polygon().centroid.coords[0]
        run("single-lake", tmp, tmp / "l.geojson", "catchment", "--dem", str(tmp / "basins"),
            "--seed", repr(lx), repr(ly), "--seed-crs", gf.EPSG, "--lakes", str(lakes),
            "--out", str(tmp / "l.geojson"))  # fmt: skip


if __name__ == "__main__":
    sys.exit(main())
