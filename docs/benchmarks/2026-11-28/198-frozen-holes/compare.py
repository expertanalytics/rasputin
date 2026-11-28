"""Increment 198's acceptance (design section 8), untimed: off is identity, and the flag's figures.

Usage (from the worktree, with its venv):
    .venv/bin/python docs/benchmarks/2026-11-28/198-frozen-holes/compare.py MASTER_SRC WORK OUT_JSON

MASTER_SRC is master 949bc7d66's ``src_python`` (``git archive 949bc7d66 src_python``)
with the head's built ``_core`` copied into ``tin_engine/`` (the increment changes no
C++). WORK is a scratch folder for inputs and meshes (kept out of the repository).

Part 1, off: every case runs once with master and once with the head, no
``--frozen-holes``, each writing ``.vtk`` (binary), ``.gpkg`` and ``--record``; the
files are compared byte for byte, the record JSON with its ``command`` entry dropped.
Master is also run twice on one case, to show the outputs are reproducible at all.

Part 2, on: the design's fixture with the flag; the ring's vertex count, the
recorded entries, the error of every DEM node off and on the hole rings against
the written mesh, and the hole-edge error recomputed by dense sampling.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TREE = HERE.parents[3]
sys.path.insert(0, str(TREE / "tests/python"))

import test_cli_mesh_frozen_holes as fx  # noqa: E402
from vtkread import polygons_as_array, read_vtk  # noqa: E402

PY = sys.executable
HEAD = [str(Path(PY).parent / "rasputin")]
# The data folder: $RASPUTIN_DATA, else rasputin_data beside the main checkout
# (TREE is a worktree under <checkout>/.claude/worktrees/).
DATA = Path(os.environ.get("RASPUTIN_DATA", TREE.parents[3] / "rasputin_data"))


def bench_quick_cases() -> list[tuple[str, list[str]]]:
    import tomllib

    doc = tomllib.loads((TREE / "docs/benchmarks/quick/cases.toml").read_text())
    out = []
    for case in doc["case"]:
        args = [a.replace("$TREE", str(TREE)).replace("$RASPUTIN_DATA", str(DATA)) for a in case["args"]]
        for t in case["threads"]:
            extra = [] if t == 0 else ["--threads", str(t)]
            out.append((f"{case['name']}-t{t}", [*args[1:], *extra]))
    return out


def fixture_cases(inp: fx.Inputs) -> list[tuple[str, list[str]]]:
    tol = ["--tolerance", "0.1"]
    return [
        ("fx-square", ["--dem", str(inp.dem), "--domain", str(inp.square), *tol]),
        ("fx-two-holes", ["--dem", str(inp.dem), "--domain", str(inp.two_holes), *tol]),
        ("fx-off-grid", ["--dem", str(inp.dem), "--domain", str(inp.off_grid), *tol]),
        ("fx-no-hole", ["--dem", str(inp.dem), "--domain", str(inp.no_hole), *tol]),
        ("fx-square-budget", ["--dem", str(inp.dem), "--domain", str(inp.square), "--triangles", "3000"]),
        ("fx-square-line-through", ["--dem", str(inp.dem), "--domain", str(inp.square), *tol,
                                    "--lines", f"wall={inp.through_vertices}", "--edges", "wall"]),
        ("fx-square-line-between", ["--dem", str(inp.dem), "--domain", str(inp.square), *tol,
                                    "--lines", f"wall={inp.between_vertices}", "--edges", "wall"]),
        ("fx-nodata", ["--dem", str(inp.nodata_dem), "--domain", str(inp.square), *tol]),
    ]


def run(cmd: list[str], out: Path, args: list[str], ext: str) -> tuple[int, str, float]:
    full = [*cmd, "mesh", *args, "--out", str(out.with_suffix(ext)), "--record", str(out.with_suffix(ext + ".json"))]
    t = time.monotonic()
    p = subprocess.run(full, capture_output=True, text=True)
    return p.returncode, p.stderr[-400:], time.monotonic() - t


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def record_minus_command(path: Path) -> str | None:
    if not path.exists():
        return None
    doc = json.loads(path.read_text())
    doc.pop("command", None)
    return hashlib.sha256(json.dumps(doc, sort_keys=True).encode()).hexdigest()


def off_identity(trees: dict[str, list[str]], work: Path, cases: list[tuple[str, list[str]]]) -> list[dict]:
    """Both trees write to the same path (the gpkg stores the command, paths included), kept after under work/<tree>/."""
    (a, b) = trees
    rows = []
    for name, args in cases:
        row: dict = {"case": name}
        for ext in (".vtk", ".gpkg"):
            got = {}
            for tree, cmd in trees.items():
                out = work / "out" / name
                out.parent.mkdir(parents=True, exist_ok=True)
                code, err, secs = run(cmd, out, args, ext)
                got[tree] = (code, sha(out.with_suffix(ext)), record_minus_command(out.with_suffix(ext + ".json")))
                row[f"{tree}{ext}_exit"] = code
                row[f"{tree}{ext}_s"] = round(secs, 2)
                if code:
                    row[f"{tree}{ext}_stderr"] = err
                keep = work / tree
                keep.mkdir(parents=True, exist_ok=True)
                for f in (out.with_suffix(ext), out.with_suffix(ext + ".json")):
                    if f.exists():
                        f.replace(keep / f.name)
            row[f"{ext[1:]}_sha_{a}"] = got[a][1]
            row[f"{ext[1:]}_equal"] = got[a][:2] == got[b][:2]
            row[f"{ext[1:]}_record_equal"] = got[a][2] == got[b][2]
        print(json.dumps(row), flush=True)
        rows.append(row)
    return rows


def bilinear(z: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """The fixture DEM at local (x, y): nodes every 10 m, row 0 at y = 1,000."""
    c, r = x / 10.0, (1000.0 - y) / 10.0
    c0, r0 = np.clip(np.floor(c).astype(int), 0, 99), np.clip(np.floor(r).astype(int), 0, 99)
    fc, fr = c - c0, r - r0
    return ((1 - fr) * ((1 - fc) * z[r0, c0] + fc * z[r0, c0 + 1])
            + fr * ((1 - fc) * z[r0 + 1, c0] + fc * z[r0 + 1, c0 + 1]))


def on_segment(p: np.ndarray, a: np.ndarray, b: np.ndarray, eps: float) -> np.ndarray:
    ab = b - a
    t = np.clip(((p - a) @ ab) / (ab @ ab), 0.0, 1.0)
    return np.linalg.norm(p - (a + t[:, None] * ab), axis=1) <= eps


def frozen_figures(work: Path, inp: fx.Inputs) -> list[dict]:
    z = fx.heights().astype(np.float64)
    gx, gy = np.meshgrid(np.arange(101) * 10.0, 1000.0 - np.arange(101) * 10.0)
    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    zn = z.ravel()
    runs = [
        ("square", inp.square, [fx.SQUARE_HOLE], ["--tolerance", "0.1"]),
        ("two-holes", inp.two_holes, [fx.SQUARE_HOLE, fx.DIAMOND], ["--tolerance", "0.1"]),
        ("off-grid", inp.off_grid, [fx.OFF_GRID_HOLE], ["--tolerance", "0.1"]),
        ("square-budget", inp.square, [fx.SQUARE_HOLE], ["--triangles", "3000"]),
    ]
    rows = []
    for name, domain, holes, extra in runs:
        out = work / "on" / name
        out.parent.mkdir(parents=True, exist_ok=True)
        code, err, _ = run(HEAD, out, ["--dem", str(inp.dem), "--domain", str(domain), *extra, "--frozen-holes"], ".vtk")
        assert code == 0, err
        rec = json.loads(out.with_suffix(".vtk.json").read_text())
        vtk = read_vtk(out.with_suffix(".vtk").read_bytes())
        pts = vtk.points[:, :2] - np.array([fx.X0, fx.Y0])
        tri = polygons_as_array(vtk)
        on_ring = np.zeros(len(nodes), bool)
        ring_vertices, moved = 0, 0.0
        dense_max = 0.0
        for hole in holes:
            ring = np.array(hole) - np.array([fx.X0, fx.Y0])
            hit = np.zeros(len(pts), bool)
            for a, b in zip(ring, np.roll(ring, -1, axis=0)):
                hit |= on_segment(pts, a, b, fx.ON_RING)
                on_ring |= on_segment(nodes, a, b, 1e-6)
                # The edge's two written ends; the chord between them, against the DEM, every 1 mm.
                ia = int(np.argmin(np.linalg.norm(pts - a, axis=1)))
                ib = int(np.argmin(np.linalg.norm(pts - b, axis=1)))
                t = np.linspace(0.0, 1.0, int(np.linalg.norm(b - a) / 1e-3) + 1)
                p = pts[ia] + t[:, None] * (pts[ib] - pts[ia])
                chord = vtk.points[ia, 2] + t * (vtk.points[ib, 2] - vtk.points[ia, 2])
                dense_max = max(dense_max, float(np.abs(bilinear(z, p[:, 0], p[:, 1]) - chord).max()))
            ring_vertices += int(hit.sum())
            near = np.linalg.norm(pts[hit][:, None, :] - ring[None, :, :], axis=2).min(axis=1)
            moved = max(moved, float(near.max()))
        # Each DEM node's error against the plane of the written triangle holding it.
        err = np.full(len(nodes), np.nan)
        P = vtk.points
        for t3 in tri:
            a, b, c = P[t3]
            lo, hi = np.minimum(np.minimum(a, b), c)[:2] - 1e-6, np.maximum(np.maximum(a, b), c)[:2] + 1e-6
            q = nodes + np.array([fx.X0, fx.Y0])
            sel = np.flatnonzero((q[:, 0] >= lo[0]) & (q[:, 0] <= hi[0]) & (q[:, 1] >= lo[1]) & (q[:, 1] <= hi[1]))
            if not len(sel):
                continue
            m = np.array([[b[0] - a[0], c[0] - a[0]], [b[1] - a[1], c[1] - a[1]]])
            uv = np.linalg.solve(m, (q[sel] - a[:2]).T).T
            inside = (uv[:, 0] >= -1e-9) & (uv[:, 1] >= -1e-9) & (uv.sum(axis=1) <= 1 + 1e-9)
            s = sel[inside]
            w = uv[inside]
            zm = a[2] + w[:, 0] * (b[2] - a[2]) + w[:, 1] * (c[2] - a[2])
            err[s] = np.abs(zn[s] - zm)
        meshed = ~np.isnan(err)
        rows.append({
            "run": name, "args": extra, "triangles": len(tri),
            "hole_vertices_in": sum(len(h) for h in holes), "ring_vertices_out": ring_vertices,
            "farthest_ring_vertex_from_input_m": moved,
            "nodes_off_rings": int((meshed & ~on_ring).sum()), "nodes_on_rings": int((meshed & on_ring).sum()),
            "max_error_off_rings_m": float(err[meshed & ~on_ring].max()),
            "max_error_on_ring_nodes_m": float(err[meshed & on_ring].max()),
            "hole_edge_error_dense_1mm_m": dense_max,
            **{k: rec.get(k) for k in ("max_error_m", "frozen_holes", "frozen_edges_max_error_m",
                                         "hole_vertices_snapped_max_m")},
        })
        print(json.dumps(rows[-1]), flush=True)
    return rows


def main() -> None:
    master = [PY, str(HERE / "run_tree.py"), sys.argv[1]]
    work = Path(sys.argv[2])
    (work / "inputs").mkdir(parents=True, exist_ok=True)
    inp = fx.Inputs(work / "inputs")
    fixtures = fixture_cases(inp)
    result = {
        "head": subprocess.run(["git", "-C", str(TREE), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
        "master": "949bc7d66",
        "on": frozen_figures(work, inp),
        "repeat": off_identity({"master": master, "master2": master}, work / "repeat",
                               fixtures[:1] + bench_quick_cases()[:1]),
        "off": off_identity({"master": master, "head": HEAD}, work / "off", fixtures + bench_quick_cases()),
    }
    Path(sys.argv[3]).write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
