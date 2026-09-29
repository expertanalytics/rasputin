"""Increment 16c acceptance on Bygdin (@perf, 2026-09-29): every table in
README.md comes from this script. Usage, from the repository root, .venv active:

    python docs/benchmarks/2026-09-29/bygdin-landcover/analyse.py <SCRATCH>

<SCRATCH> holds run.sh's meshes. Prints to stdout; run.sh's logs are read from
./logs. The CORINE reference is read here straight from the GeoPackage with
sqlite3 and shapely, not through tin_engine.feature_input, so it is independent
of the code under test.
"""

from __future__ import annotations

import re
import sqlite3
import statistics
import sys
from pathlib import Path

import numpy as np
import shapely
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOLegacy import vtkPolyDataReader

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "tools"))
import bench  # noqa: E402

from tin_engine.landcover import regions  # noqa: E402

GPKG = ROOT.parent / "rasputin_data" / "corine2018_dtm10_utm33.gpkg"
MARGIN = 2 * 1e-3  # cli: margin = 2 * snap spacing, default 1 mm
# 22's table (bygdin/README.md on increment22-autocatchment), percent.
TABLE_22 = {333: 42.34, 332: 21.32, 322: 17.00, 512: 16.46, 335: 2.44, 412: 0.33, 142: 0.10}


def gpkg_geometry(blob: bytes) -> shapely.Geometry:
    """A GeoPackage binary geometry: 8-byte header, envelope, then WKB."""
    flags = blob[3]
    envelope = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[(flags >> 1) & 0b111]
    return shapely.from_wkb(blob[8 + envelope :])


def corine(domain: shapely.Geometry) -> list[tuple[shapely.Geometry, int]]:
    """Every CORINE polygon whose envelope meets the domain's, unclipped."""
    x0, y0, x1, y1 = domain.bounds
    con = sqlite3.connect(f"file:{GPKG}?mode=ro", uri=True)
    rows = con.execute(
        "select c.geom, c.code_18 from corine2018 c join rtree_corine2018_geom r on c.fid = r.id "
        "where r.maxx >= ? and r.minx <= ? and r.maxy >= ? and r.miny <= ?",
        (x0, x1, y0, y1),
    ).fetchall()
    out = [(gpkg_geometry(g), int(code)) for g, code in rows]
    return [(g, c) for g, c in out if g.intersects(domain)]


def read(path: Path) -> dict:
    r = vtkPolyDataReader()
    r.SetFileName(str(path))
    r.ReadAllFieldsOn()
    r.Update()
    pd = r.GetOutput()
    arr = pd.GetCellData().GetArray("land_cover_code")
    n_lines, n_polys = pd.GetNumberOfLines(), pd.GetNumberOfPolys()
    lines = vtk_to_numpy(pd.GetLines().GetConnectivityArray()).reshape(-1, 2)
    tris = vtk_to_numpy(pd.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    fd = pd.GetFieldData().GetAbstractArray("land_cover_codes")
    return dict(
        cells=pd.GetNumberOfCells(), n_lines=n_lines, n_polys=n_polys,
        codes=None if arr is None else vtk_to_numpy(arr).astype(np.int64),
        points=vtk_to_numpy(pd.GetPoints().GetData()).astype(np.float64),
        tris=tris.astype(np.int64), lines=lines.astype(np.int64),
        text=None if fd is None else fd.GetValue(0),
        scalars=pd.GetCellData().GetScalars().GetName() if pd.GetCellData().GetScalars() else None,
    )  # fmt: skip


def rounds(tri: np.ndarray, edges: np.ndarray) -> tuple[int, np.ndarray]:
    """landcover.regions' loop, counting its hook rounds (outer iterations
    that hooked something). Checked equal to regions() by the caller."""
    n = int(max(tri.max(), edges.max())) + 1
    key = lambda p: np.minimum(p[:, 0], p[:, 1]) * n + np.maximum(p[:, 0], p[:, 1])  # noqa: E731
    sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    keys, owner = key(sides), np.tile(np.arange(len(tri)), 3)
    free = ~np.isin(keys, key(edges))
    order = np.argsort(keys[free], kind="stable")
    keys, owner = keys[free][order], owner[free][order]
    pair = np.flatnonzero(keys[1:] == keys[:-1])
    u, v = owner[pair], owner[pair + 1]
    parent, count = np.arange(len(tri)), 0
    while True:
        ru, rv = parent[u], parent[v]
        if np.array_equal(ru, rv):
            return count, parent
        count += 1
        np.minimum.at(parent, np.maximum(ru, rv), np.minimum(ru, rv))
        while not np.array_equal(parent, jumped := parent[parent]):
            parent = jumped


def oracle(xy, tri, polys) -> tuple[np.ndarray, np.ndarray]:
    """The design's oracle (R1): each triangle's centroid against the
    polygons, smallest area wins, ties to the smaller code; and the mask of
    triangles with inradius > 1.5 * margin, where it is exact."""
    a, b, c = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
    per = sum(np.hypot(*(q - p).T) for p, q in ((b, c), (c, a), (a, b)))
    cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    r = np.abs(cross) / per
    tree = shapely.STRtree([p for p, _ in polys])
    found, which = tree.query(shapely.points((a + b + c) / 3), predicate="intersects")
    area = np.array([p.area for p, _ in polys])
    code = np.array([k for _, k in polys])
    order = np.lexsort((code[which], area[which], found))
    found, which = found[order], which[order]
    head = np.unique(found, return_index=True)[1]
    out = np.zeros(len(tri), dtype=np.int64)
    out[found[head]] = code[which[head]]
    return out, r > 1.5 * MARGIN


def stats_phases(path: Path) -> dict[str, float]:
    text = path.read_text()
    out = {m[1]: float(m[2]) for m in re.finditer(r"^\| ([^|]+?) \| ([0-9.]+) \| [0-9.]+ % \|$", text, re.M)}
    total = re.search(r"^\| \*?\*?total\*?\*? \| \*?\*?([0-9.]+)", text, re.M)
    if total:
        out["total"] = float(total[1])
    return out


def wall(path: Path) -> tuple[float, int]:
    text = path.read_text()
    return float(re.search(r"([0-9.]+) real", text)[1]), int(re.search(r"(\d+)\s+maximum resident", text)[1])


def main(scratch: Path) -> None:
    domain = shapely.from_geojson((HERE / "bygdin_reduced_t20.geojson").read_text())
    domain = shapely.union_all([g for g in shapely.get_parts(domain)]) if domain.geom_type == "GeometryCollection" else domain
    polys = corine(domain)
    ref: dict[int, float] = {}
    for p, k in polys:
        ref[k] = ref.get(k, 0.0) + shapely.intersection(p, domain).area
    ref_total = sum(ref.values())
    print(f"CORINE polygons meeting the domain: {len(polys)}; clipped area {ref_total / 1e6:.6f} km2, domain {domain.area / 1e6:.6f} km2")

    for tol in ("10", "1"):
        tag = f"lc_t{tol}"
        print(f"\n== {tol} m ==")
        m = read(scratch / f"{tag}_run1.vtk")
        codes, tri, lines, pts = m["codes"], m["tris"], m["lines"], m["points"]
        print(f"vtkPolyDataReader: {m['cells']} cells = {m['n_lines']} lines + {m['n_polys']} triangles; active SCALARS {m['scalars']}")
        assert codes is not None, "land_cover_code missing"
        print(f"land_cover_code: {len(codes)} values (one per cell: {len(codes) == m['cells']}); "
              f"nonzero on lines: {int((codes[: m['n_lines']] != 0).sum())}")  # fmt: skip
        print(f"FieldData land_cover_codes: {m['text']!r}")
        tc = codes[m["n_lines"] :]
        xy = pts[:, :2]
        a, b, c = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
        area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
        total = area.sum()
        print(f"mesh area {total / 1e6:.6f} km2; triangles with code 0: {int((tc == 0).sum())}")
        print("| code | mesh km2 | mesh share % | CORINE clipped share % | diff (pp) | 22's table % |")
        worst = 0.0
        for k in sorted(set(ref) | set(np.unique(tc).tolist()), key=lambda k: -ref.get(k, 0)):
            ms = 100 * area[tc == k].sum() / total
            rs = 100 * ref.get(k, 0.0) / ref_total
            worst = max(worst, abs(ms - rs))
            print(f"| {k} | {area[tc == k].sum() / 1e6:.6f} | {ms:.5f} | {rs:.5f} | {ms - rs:+.6f} | {TABLE_22.get(k, '-')} |")
        print(f"largest share difference: {worst:.6f} pp; codes equal: {set(np.unique(tc).tolist()) == set(ref)}")
        # I1: across every interior non-constraint edge, the codes agree.
        n = len(pts)
        key = lambda p: np.minimum(p[:, 0], p[:, 1]) * n + np.maximum(p[:, 0], p[:, 1])  # noqa: E731
        sides = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
        k3, own = key(sides), np.tile(np.arange(len(tri)), 3)
        free = ~np.isin(k3, key(lines))
        o = np.argsort(k3[free], kind="stable")
        ks, ow = k3[free][o], own[free][o]
        pr = np.flatnonzero(ks[1:] == ks[:-1])
        print(f"I1: {len(pr)} interior unconstrained edges, {int((tc[ow[pr]] != tc[ow[pr + 1]]).sum())} with differing codes")
        # I2: the centroid oracle.
        want, exact = oracle(xy, tri, polys)
        print(f"I2: oracle checks {int(exact.sum())} of {len(tri)} triangles (r > {1.5 * MARGIN} m); "
              f"disagreements {int((want[exact] != tc[exact]).sum())}; "
              f"disagreements on the {int((~exact).sum())} unchecked: {int((want[~exact] != tc[~exact]).sum())}")  # fmt: skip
        nr, parent = rounds(tri, lines)
        same = np.array_equal(parent, regions(tri, lines))
        print(f"union-find: {nr} hook rounds; {len(np.unique(parent))} components; equal to landcover.regions: {same}")
        # The ASCII run: same codes, and bench.quality.
        asc = read(scratch / f"{tag}_ascii.vtk")
        print(f"ASCII run codes equal to binary run1: {np.array_equal(asc['codes'], codes)}")
        vm = bench.read_vtk_ascii(scratch / f"{tag}_ascii.vtk")
        err = (HERE / "logs" / f"{tag}_ascii.err").read_text()
        max_err = float(re.search(r"achieved max error ([0-9.e+-]+)", err)[1])
        q = bench.quality(vm.points, vm.triangles, vm.edges, float(tol), max_err)
        print(f"quality: worst angle {q.worst_angle:.5f} deg, max degree {q.max_degree}, "
              f"max error {max_err} (within: {q.within_tolerance}), Delaunay checked/ambiguous/violations "
              f"{q.delaunay_checked}/{q.delaunay_ambiguous}/{q.delaunay_violations}")  # fmt: skip
        # Times over the three binary runs.
        runs = [stats_phases(HERE / "logs" / f"{tag}_run{i}.stats.md") for i in (1, 2, 3)]
        walls = [wall(HERE / "logs" / f"{tag}_run{i}.err") for i in (1, 2, 3)]
        for ph in ("land cover", "refine", "features clip", "total"):
            vals = [r.get(ph, float("nan")) for r in runs]
            print(f"time {ph}: median {statistics.median(vals):.3f} s (runs {', '.join(f'{v:.3f}' for v in vals)})")
        print(f"wall: median {statistics.median(w for w, _ in walls):.2f} s (runs {', '.join(f'{w:.2f}' for w, _ in walls)}); "
              f"max RSS median {statistics.median(r for _, r in walls) / 1e6:.0f} MB")  # fmt: skip
        for i in (1, 2, 3):
            line = [ln for ln in (HERE / "logs" / f"{tag}_run{i}.err").read_text().splitlines() if ln.startswith("land cover")]
            print(f"run{i} stderr: {line[0]}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
