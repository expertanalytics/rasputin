"""Emulate increment 111 on Sula at 20 m without changing code.

111 builds the land-cover map once (CORINE with NVE's lakes preferred, the
refilled pieces joined to their neighbours) and takes both the edges and the
class fractions from it. Emulated in two steps:

1. Today's direct command, run in this process, with the step that builds
   the land cover's lines (`_Tally._preferred`) wrapped so that its output,
   the final map as the lines see it, is written to a GeoJSON in EPSG:25833.
2. That map given as the only land cover: `--features map --land-cover map`,
   no `--prefer`, every class an edge (Sol's rules file), the same clean-up.
3. The map as 111 builds it: CORINE read as the class fractions read it
   (`read_class_polygons`, which repairs self-crossing rings), over the read
   region, then the same cut, refill and join; meshed as in step 2.

Prints triangles, lines and how far the class fractions of the two runs
differ, triangle by triangle. Run:
.venv/bin/python -I docs/increments/111-probes/one_map_probe.py <work dir>
(needs the clc2018 and nve-lakes caches, dtm10-norway, and the Sula files
under ~/rasputin-output/mcp/ola/work/).
"""

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import numpy as np
import shapely
from typer.testing import CliRunner

from tin_engine import cli, feature_input
from tin_engine.feature_input import read_class_polygons, read_polygons, source_region
from tin_engine.io.domain_file import read_domain
from tin_engine.meshing_rules import load

HOME = Path.home()
WORK = HOME / "rasputin-output/mcp/ola/work"
DOMAIN = WORK / "sula-73-27-0/sula-official-catchment.geojson"
CLC = HOME / ".cache/rasputin/clc2018/features.gpkg"
RULES = WORK / "uploads/sula-corine-constraints.toml"
NVE = HOME / ".cache/rasputin/nve-lakes/features.gpkg"
CRS = "EPSG:25833"
RASPUTIN = str(Path(sys.executable).parent / "rasputin")


def common(out: Path) -> list[str]:
    return [
        "mesh", "--out", str(out), "--dem", "dtm10-norway", "--tolerance", "20",
        "--domain", str(DOMAIN), "--out-crs", CRS, "--features-rules", str(RULES),
        "--features-repair", "0.1", "--features-tolerance", "20",
        "--features-outline-snap", "0", "--threads", "4", "--memory-budget", "4G",
        "--z-unit", "m",
    ]  # fmt: skip


def read(path: Path) -> tuple[np.ndarray, np.ndarray, dict[int, dict[int, float]]]:
    """Triangle polygons, line geometries and each triangle's fractions."""
    con = sqlite3.connect(path)

    def geoms(table: str) -> np.ndarray:
        out = []
        for (b,) in con.execute(f"select geom from {table} order by fid"):
            envelope = (0, 32, 48, 48, 64)[(b[3] >> 1) & 7]
            out.append(shapely.from_wkb(bytes(b[8 + envelope :])))
        return np.array(out, dtype=object)

    shares: dict[int, dict[int, float]] = {}
    for t, c, f in con.execute("select triangle, code, fraction from land_cover_fractions"):
        shares.setdefault(t, {})[c] = f
    tri, lines = geoms("triangles"), geoms("lines")
    con.close()
    return tri, lines, shares


def write_map(path: Path, polygons: list[tuple[str, object]]) -> None:
    features = [
        {"type": "Feature", "properties": {"Code_18": c}, "geometry": g.__geo_interface__}
        for c, g in polygons
    ]
    named = {"type": "name", "properties": {"name": CRS}}
    path.write_text(json.dumps({"type": "FeatureCollection", "crs": named, "features": features}))


def one_map(work: Path, final: Path, name: str) -> Path:
    new = work / f"{name}.gpkg"
    subprocess.run(
        [RASPUTIN, *common(new), "--features", str(final), "--land-cover", str(final),
         "--land-cover-rules", "corine"],
        check=True, capture_output=True,
    )  # fmt: skip
    return new


def built_map() -> list[tuple[str, object]]:
    """Step 3: the map from the repaired CORINE polygons of the read region."""
    domain = read_domain(DOMAIN)  # in EPSG:25833, the DEM's CRS
    region = source_region(domain, CRS, CRS, cell=10.0)
    xy = shapely.get_coordinates(region.exterior)
    corine = load("corine")
    base, _ = read_class_polygons(CLC, None, corine, CRS, xy)
    lakes, _ = read_polygons(NVE, None, CRS, shapely.get_coordinates(domain.polygon.exterior))
    dem = feature_input.parse_crs(CRS)
    tally = feature_input._Tally(domain, dem, feature_input.DEFAULT_VOCABULARY)
    cover = [(i, 1, None, str(c), g) for i, (c, g) in enumerate(base)]
    pre = feature_input.Prefer(source=0, rules=corine, chosen=[(512, lakes)])
    out = tally._preferred(pre, cover, [c for c, _ in base], region, corine)
    return [(str(item[3]), item[4]) for item in out]


def main(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    captured: list[tuple[str, object]] = []
    inner = feature_input._Tally._preferred

    def wrapped(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        out = inner(self, *args, **kwargs)
        captured.extend((str(item[3]), item[4]) for item in out)  # merge key = Code_18
        return out

    feature_input._Tally._preferred = wrapped  # type: ignore[method-assign]
    today = work / "today.gpkg"
    args = [*common(today), "--features", str(CLC), "--land-cover", "clc2018",
            "--prefer", "lake=nve-lakes"]  # fmt: skip
    done = CliRunner().invoke(cli.app, args)
    assert done.exit_code == 0, done.output
    feature_input._Tally._preferred = inner  # type: ignore[method-assign]
    final = work / "final-map.geojson"
    write_map(final, captured)
    print(f"map: {len(captured)} polygons")
    t0, l0, s0 = read(today)
    print("step 2:")
    compare(t0, l0, s0, *read(one_map(work, final, "one-map")))
    built = work / "built-map.geojson"
    polygons = built_map()
    write_map(built, polygons)
    print(f"step 3, map: {len(polygons)} polygons")
    compare(t0, l0, s0, *read(one_map(work, built, "built-map")))


def compare(t0, l0, s0, t1, l1, s1) -> None:  # type: ignore[no-untyped-def]
    """Triangles, line vertices, and fractions on the triangles both share."""

    def corners(g: object) -> tuple[tuple[float, float], ...]:
        return tuple(sorted(map(tuple, np.round(shapely.get_coordinates(g)[:3], 6))))

    k0 = {corners(g): i for i, g in enumerate(t0)}
    k1 = {corners(g): i for i, g in enumerate(t1)}
    v0 = {tuple(p) for p in np.round(shapely.get_coordinates(l0), 6)}
    v1 = {tuple(p) for p in np.round(shapely.get_coordinates(l1), 6)}
    moved = [min(np.hypot(p[0] - q[0], p[1] - q[1]) for q in v0 - v1) for p in v1 - v0]
    print(f"triangles: today {len(t0)}, one map {len(t1)}, not shared {len(k0.keys() - k1.keys())}")
    print(f"line vertices: today {len(v0)}, one map {len(v1)}, not shared {len(v0 - v1)}, "
          f"largest move {max(moved, default=0.0):.4f} m")  # fmt: skip
    worst = 0.0
    for k, i in k0.items():
        if (j := k1.get(k)) is not None:
            a, b = s0.get(i + 1, {}), s1.get(j + 1, {})  # fractions number triangles from 1
            worst = max([worst, *(abs(a.get(c, 0.0) - b.get(c, 0.0)) for c in a.keys() | b.keys())])
    print(f"largest fraction difference on shared triangles: {worst:.3e}")

if __name__ == "__main__":
    main(Path(sys.argv[1]))
