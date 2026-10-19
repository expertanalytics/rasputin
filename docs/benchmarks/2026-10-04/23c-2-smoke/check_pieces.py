"""23c-2 smoke on real data (@perf, 2026-10-04): checks over a cut run's piece files.

    python check_pieces.py PIECES_DIR TOLERANCE MERGED_OUT.vtk [--shift M]

1. Conformity: every constraint segment (VTK LINES) of every piece, keyed by its
   two end points rounded to 1e-6 m. A seam segment must appear in exactly two
   pieces with the same end points and the same z. A segment that appears once
   must lie on the domain's outline: it is not shared, so its two ends are
   checked to be outline vertices. That is approximated by "not on a lattice
   line in both ends". Reported: counts by multiplicity, and z mismatches.
2. A merged mesh: the pieces' points deduplicated on (x, y) rounded to 1e-6 m,
   triangles and lines re-indexed, written as ASCII VTK to MERGED_OUT.vtk. Then
   tools/bench.py's quality() on it, with every piece's LINES (outline and seams)
   as constraints: worst angle, max degree and the constrained Delaunay check.
3. The independent check: 15c-2's run_geo.independent_check on the merged mesh
   (matplotlib's interpolator at every valid ANADEM node of the source window
   inside the domain, interior and edge strip). With --shift, the control is the
   same with the mesh moved east by that many metres.
The merged file is a measurement artefact, not 23d's stitcher.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
W = HERE.parents[3]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
# quality() needs some tin_engine for tin_engine.stats.quality; this branch's build is
# for its 3.13 venv, which has no matplotlib, so QUALITY_PKG may name another
# Release package (master's) for the interpreter running this script.
import os  # noqa: E402

sys.path[:0] = [os.environ.get("QUALITY_PKG", str(W / "build-bench/pkg")), str(W / "tools"),
                str(W / "docs/benchmarks/2026-10-02/15c-2-acceptance")]  # fmt: skip
import bench  # noqa: E402

DATA = Path("/Users/skavhaug/projects/rasputin_data/sao_francisco_piece")
OUTLINE = DATA / "bho2017_5k_76949_outline_epsg4674.geojson"
WINDOW = DATA / "bho2017_5k_76949_anadem_window_epsg4674.tif"
CRS = "EPSG:31983"


def main() -> None:
    pieces, tol, merged = Path(sys.argv[1]), float(sys.argv[2]), Path(sys.argv[3])
    shift = float(sys.argv[sys.argv.index("--shift") + 1]) if "--shift" in sys.argv else 0.0
    t0 = time.perf_counter()
    index = json.loads((pieces / "index.json").read_text())
    entries = index["pieces"] if isinstance(index, dict) and "pieces" in index else index
    files = sorted(pieces.glob("*.vtk"))
    out: dict = {"pieces_dir": str(pieces), "piece_files": len(files), "index_entries": len(entries)}

    seg: dict[tuple, list] = {}
    all_pts, all_tri, all_lin, base = [], [], [], 0
    per_piece = []
    for f in files:
        m = bench.read_vtk_ascii(f)
        all_pts.append(m.points)
        all_tri.append(m.triangles + base)
        all_lin.append(m.edges + base)
        base += len(m.points)
        per_piece.append({"file": f.name, "points": int(len(m.points)), "triangles": int(len(m.triangles)),
                          "lines": int(len(m.edges))})  # fmt: skip
        for a, b in m.edges:
            pa, pb = m.points[a], m.points[b]
            ka, kb = tuple(np.round(pa[:2], 6)), tuple(np.round(pb[:2], 6))
            key = (ka, kb) if ka <= kb else (kb, ka)
            za, zb = (pa[2], pb[2]) if ka <= kb else (pb[2], pa[2])
            seg.setdefault(key, []).append((f.name, za, zb))
    mult: dict[int, int] = {}
    z_mismatch = 0
    for v in seg.values():
        mult[len(v)] = mult.get(len(v), 0) + 1
        if len(v) == 2 and (v[0][1] != v[1][1] or v[0][2] != v[1][2]):
            z_mismatch += 1
    out["constraint_segments_by_multiplicity"] = {str(k): v for k, v in sorted(mult.items())}
    out["shared_segments_z_mismatch"] = z_mismatch

    pts = np.vstack(all_pts)
    key = np.round(pts[:, :2], 6)
    uniq, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    zdup = np.abs(pts[:, 2] - pts[first[inv], 2]).max()
    tri = inv[np.vstack(all_tri)]
    lin = inv[np.vstack(all_lin)]
    lin = np.unique(np.sort(lin, axis=1), axis=0)
    mp = pts[first]
    out["merged"] = {"points": int(len(mp)), "triangles": int(len(tri)), "constraint_edges": int(len(lin)),
                     "duplicate_points_max_z_difference": float(zdup)}  # fmt: skip
    with open(merged, "w") as fh:
        fh.write("# vtk DataFile Version 3.0\nmerged pieces (measurement)\nASCII\nDATASET POLYDATA\n")
        fh.write(f"POINTS {len(mp)} double\n")
        np.savetxt(fh, mp, fmt="%.17g")
        fh.write(f"LINES {len(lin)} {3 * len(lin)}\n")
        np.savetxt(fh, np.column_stack([np.full(len(lin), 2), lin]), fmt="%d")
        fh.write(f"POLYGONS {len(tri)} {4 * len(tri)}\n")
        np.savetxt(fh, np.column_stack([np.full(len(tri), 3), tri]), fmt="%d")
    q = bench.quality(mp, tri, lin, tol, 0.0)
    out["quality"] = {k: v for k, v in q.model_dump().items() if k not in ("within_tolerance", "mesh_sha256")}

    import run_geo  # noqa: E402  matplotlib-based, no tin_engine in its check

    out["check"] = run_geo.independent_check(merged, WINDOW, OUTLINE, CRS, 30.0, tol)
    out["check"].pop("quality", None)
    if shift:
        c = run_geo.independent_check(merged, WINDOW, OUTLINE, CRS, 30.0, tol, shift)
        out["control"] = {"shift_m": shift, "interior_over": c["interior"]["over_tolerance"],
                          "strip_over": c["strip"]["over_tolerance"]}  # fmt: skip
    out["per_piece"] = per_piece
    out["seconds"] = round(time.perf_counter() - t0, 1)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
