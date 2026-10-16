"""Summarise the land-cover simplification probe (probe.sh).

    /Users/skavhaug/projects/rasputin/.venv/bin/python summarise.py [MESHDIR] > summary.md

Per run: start and output triangles, max and RMS height error (from
*_errors.json, Ola's mesh_error_stats.py), triangles under 1 degree counted
from the .vtk, worst angle and land-cover vertices after clean-up and moved
area (from --stats), clean-up seconds and total seconds, and every timing row
at 40 % or more. Then, per --features-tolerance, the area of each CORINE class
in the minimal mesh against the source clipped to the outline (EPSG:25832),
and the area labelled with another class than the source's (exact: the mesh's
class polygons intersected with the source's). Needs matplotlib-free deps
only: numpy, shapely, pyproj, and lands/scripts/vtk_to_gpkg.read_vtk.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import shape

sys.path.insert(0, "/Users/skavhaug/projects/rasputin_data/lands/scripts")
from vtk_to_gpkg import read_vtk  # noqa: E402

HERE = Path(__file__).parent
MESH = Path(sys.argv[1] if len(sys.argv) > 1 else "/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify")
OUTLINE = "/Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/fused_catchment_glo30.geojson"
CORINE = "/Users/skavhaug/projects/rasputin_data/germany_corine/isar_loisach_clc2018_3035.geojson"
FTS = (0, 10, 30, 100)
RUNS = ["minimal"] + [f"tol{t}_a{a}" for t in (50, 20, 10) for a in (25, 0)]


def row(text: str, label: str) -> str:
    m = re.search(rf"^\| {re.escape(label)} \| ([^|]+) \|", text, re.M)
    return m.group(1).strip() if m else ""


def angles(pts: np.ndarray, tris: np.ndarray) -> np.ndarray:
    """Smallest plan angle per triangle, degrees."""
    p = pts[:, :2][tris]
    out = []
    for i in range(3):
        u = p[:, (i + 1) % 3] - p[:, i]
        v = p[:, (i + 2) % 3] - p[:, i]
        c = (u * v).sum(1) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1))
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return np.min(out, axis=0)


def runs_table() -> None:
    print("| FT m | run | start tris | tris | max err m | RMS m | < 1° | worst | LC vertices after | moved m² | clean-up s | total s | phases ≥ 40 % |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for ft in FTS:
        for r in RUNS:
            n = f"ft{ft}_{r}"
            st = (HERE / f"{n}_stats.md").read_text()
            er = json.loads((HERE / f"{n}_errors.json").read_text())
            pts, tris, _, _ = read_vtk(MESH / f"{n}.vtk")
            a = angles(pts, tris)
            lc = row(st, "Land-cover vertices before and after clean-up")
            lc_after = re.search(r"(\d+) after clean-up", lc).group(1)
            worst = re.search(r"^\| minimum angle \|.*\| ([0-9.]+°) \|$", st, re.M).group(1)
            total = re.search(r"\*\*total\*\* \| \*\*([0-9.]+)\*\*", st).group(1)
            big = [f"{m.group(1)} {m.group(3)} %" for m in re.finditer(r"^\| ([^|*]+?) \| ([0-9.]+) \| ([0-9.]+) % \|$", st, re.M) if float(m.group(3)) >= 40]
            print(f"| {ft} | {r} | {row(st, 'start triangles')} | {len(tris)} | {er['max_abs_m']:.1f} | {er['rms_m']:.2f} | {int((a < 1).sum())} | {worst} | {lc_after} | {row(st, 'Land-cover area inside the outline that the outline rule gave to another polygon, m2')} | {row(st, 'features clip: clean-up')} | {total} | {'; '.join(big) or '-'} |")


def class_areas() -> None:
    dom = shape(json.load(open(OUTLINE))["features"][0]["geometry"])
    to = Transformer.from_crs(3035, 25832, always_xy=True)
    src: dict[int, list] = {}
    for f in json.load(open(CORINE))["features"]:
        g = shapely.transform(shape(f["geometry"]), lambda xy: np.column_stack(to.transform(xy[:, 0], xy[:, 1])))
        src.setdefault(int(f["properties"]["Code_18"]), []).append(g)
    ref = {c: shapely.intersection(shapely.union_all(gs), dom) for c, gs in src.items()}
    ref = {c: g for c, g in ref.items() if g.area > 0}
    mesh = {}
    for ft in FTS:
        pts, tris, code, _ = read_vtk(MESH / f"ft{ft}_minimal.vtk")
        polys = shapely.polygons(pts[:, :2][tris])
        mesh[ft] = {int(c): shapely.coverage_union_all(polys[code == c]) for c in np.unique(code)}
    print("\n| class | source km² | " + " | ".join(f"FT {ft}: km² (Δ %)" for ft in FTS) + " |")
    print("|---|---|" + "---|" * len(FTS))
    for c in sorted(ref):
        a0 = ref[c].area
        cells = []
        for ft in FTS:
            a = mesh[ft].get(c, shapely.Polygon()).area
            cells.append(f"{a / 1e6:.3f} ({100 * (a - a0) / a0:+.1f})")
        print(f"| {c} | {a0 / 1e6:.3f} | " + " | ".join(cells) + " |")
    total = dom.area
    print("\n| FT m | area labelled with another class than the source's, km² | share of the outline |")
    print("|---|---|---|")
    for ft in FTS:
        wrong = sum(g.area - shapely.intersection(g, ref.get(c, shapely.Polygon())).area for c, g in mesh[ft].items())
        print(f"| {ft} | {wrong / 1e6:.3f} | {100 * wrong / total:.2f} % |")
    print(f"\nOutline area {total / 1e6:.1f} km²; {len(ref)} classes inside.")


if __name__ == "__main__":
    runs_table()
    class_areas()
