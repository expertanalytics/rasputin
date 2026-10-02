"""Tables and the basin extrapolation from run_sweep.py's and sample_boxes.py's JSON.

A measurement script, not production code; nothing imports it.

    python docs/benchmarks/2026-10-02/basin-piece-anadem/analyse.py PIECE_RESULTS PIECE_OUTLINE \
        BASIN_OUTLINE BASIN_BOXES_JSON PIECE_BOXES_JSON [SCALING_JSON]

Prints Markdown. Areas are geodesic on GRS80 (pyproj.Geod). The basin area
used is BHO's level-2 ottobasin 76 outline's; the OAS figure (636,920 km²) is
printed beside it. A box sample's interval is a percentile bootstrap of the
mean density (10,000 resamples, seed 1).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pyproj
from shapely.geometry import shape

OAS_KM2 = 636_920.0
CANVAS_GIB = 8.6  # Q9: the basin's bounding box at 30 m as float32


def area_km2(path: str) -> float:
    g = shape(json.loads(Path(path).read_text())["features"][0]["geometry"])
    return abs(pyproj.Geod(ellps="GRS80").geometry_area_perimeter(g)[0]) / 1e6


def med(xs: list[float]) -> float:
    return float(np.median(xs))


def piece_table(res: dict, a: float) -> dict:
    n = len(res["runs"][0]["timed"])
    print(f"## The piece: {res['piece']}, {a:,.1f} km² (geodesic)\n")
    print(f"Medians of {n} timed runs per tolerance; one more ASCII run for quality and the "
          "final check.\n")
    print("| tolerance | triangles | vertices | vertices / grid nodes | triangles / km² | refine "
          "| process wall | max RSS | threads |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    rows = []
    for r in res["runs"]:
        t, fc = r["timed"], r["final_check"]
        tri, ver = t[0]["triangles"], t[0]["vertices"]
        assert all(x["triangles"] == tri and x["vertices"] == ver for x in t), "runs differ"
        row = dict(tol=r["tolerance"], tri=tri, ver=ver, refine=med([x["refine_s"] for x in t]),
                   proc=med([x["proc_s"] for x in t]), rss=med([x["max_rss_bytes"] for x in t]))
        rows.append(row)
        print(f"| {r['tolerance']:g} m | {tri:,} | {ver:,} | "
              f"{ver / fc['control_resampled_grid']['nodes']:.1%} | {tri / a:,.1f} | "
              f"{row['refine']:.2f} s | {row['proc']:.2f} s | {row['rss'] / 2**30:.2f} GiB | "
              f"{t[0]['threads']} |")
    print("\n| tolerance | worst angle | max degree | Delaunay violations (checked, ambiguous) | "
          "within tolerance | control: grid nodes over tol | source nodes over tol, interior "
          "| source nodes over tol, boundary strip |\n|---:|---:|---:|---:|---|---:|---:|---:|")
    for r in res["runs"]:
        q, fc = r["quality"], r["final_check"]
        i, s = fc["source_dem_interior"], fc["source_dem_strip"]
        print(f"| {r['tolerance']:g} m | {q['worst_angle']:.3g}° | {q['max_degree']} | "
              f"{q['delaunay_violations']} ({q['delaunay_checked']:,}, {q['delaunay_ambiguous']:,}) | "
              f"{q['within_tolerance']} | {fc['control_resampled_grid']['over_tolerance']} of "
              f"{fc['control_resampled_grid']['nodes']:,} | {i['over_tolerance']:,} of "
              f"{i['nodes']:,} ({i['over_share']:.2%}; max {i['max']:.1f} m) | "
              f"{s['over_tolerance']:,} of {s['nodes']:,} ({s['over_share']:.2%}; "
              f"max {s['max']:.1f} m) |")
    power = {("AC Power" in r["pmset_before"] and "AC Power" in r["pmset_after"]) for r in res["runs"]}
    print(f"\nPower: AC before and after every tolerance block: {power == {True}}.\n")
    tri = np.array([r["tri"] for r in rows], dtype=float)
    rss = np.array([r["rss"] for r in rows], dtype=float)
    slope, base = np.polyfit(tri, rss, 1)
    print(f"Max RSS against triangles over the sweep, least squares: {base / 2**30:.2f} GiB + "
          f"{slope:.0f} B per triangle (max residual "
          f"{np.abs(rss - (base + slope * tri)).max() / 2**20:.0f} MiB).\n")
    return {"rows": rows, "rss_slope": slope, "area": a}


def boxes(path: str) -> tuple[dict, list[dict]]:
    d = json.loads(Path(path).read_text())
    return d, d["boxes"]


def density(bx: list[dict], t: str, side_km2: float, key: str = "triangles") -> np.ndarray:
    return np.array([b["tolerances"][t][key] / side_km2 for b in bx])


def interval(x: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    boot = rng.choice(x, (10_000, x.size)).mean(axis=1)
    lo, hi = np.quantile(boot, [0.025, 0.975])
    return float(lo), float(hi)


def basin_table(path: str, basin: float, piece: dict) -> None:
    d, bx = boxes(path)
    side = (d["side_m"] / 1000.0) ** 2
    rng = np.random.default_rng(1)
    print(f"## The basin, extrapolated from {len(bx)} random {d['side_m'] / 1000:g} km boxes "
          f"(seed {d['seed']})\n")
    print("Mean density over the boxes times the basin's area; the interval is the bootstrap "
          "95 % interval of the mean. Memory floor: the piece's per-triangle RSS slope "
          f"({piece['rss_slope']:.0f} B) times the basin's triangles, plus the {CANVAS_GIB} GiB "
          "canvas of Q9. Refine time: the piece's refine time per triangle (10 threads), "
          "times the basin's triangles, i.e. linear.\n")
    print("| tolerance | triangles / km², mean (95 %) | median | p10-p90 | basin triangles "
          "(95 %) | basin vertices | refine, linear | memory floor |\n"
          "|---:|---:|---:|---:|---:|---:|---:|---:|")
    by_tol = {r["tol"]: r for r in piece["rows"]}
    for t in bx[0]["tolerances"]:
        x, v = density(bx, t, side), density(bx, t, side, "vertices")
        lo, hi = interval(x, rng)
        tri = x.mean() * basin
        p = by_tol[float(t)]
        print(f"| {t} m | {x.mean():,.1f} ({lo:,.1f}-{hi:,.1f}) | {np.median(x):,.1f} | "
              f"{np.quantile(x, 0.1):,.1f}-{np.quantile(x, 0.9):,.1f} | {tri / 1e6:,.1f} M "
              f"({lo * basin / 1e6:,.1f}-{hi * basin / 1e6:,.1f} M) | "
              f"{v.mean() * basin / 1e6:,.1f} M | {p['refine'] / p['tri'] * tri / 60:,.1f} min | "
              f"{(piece['rss_slope'] * tri) / 2**30 + CANVAS_GIB:,.0f} GiB |")
    print("\nThe final check of Q6, measured: source-DEM nodes outside the tolerance after "
          "meshing the resampled grid, **interior only** (a box's 4-corner edges make its "
          "boundary strip unrepresentative of a real outline; the piece gives the strip). "
          "Basin count = mean per km² times the basin's area.\n")
    print("| tolerance | share of source nodes | per km², mean (95 %) | basin source nodes over "
          "tol | as a share of the basin's vertices |\n|---:|---:|---:|---:|---:|")
    for t in bx[0]["tolerances"]:
        over = np.array([b["tolerances"][t]["final_check"]["source_dem_interior"]["over_tolerance"]
                         for b in bx]) / side  # fmt: skip
        nodes = sum(b["tolerances"][t]["final_check"]["source_dem_interior"]["nodes"] for b in bx)
        lo, hi = interval(over, rng)
        verts = density(bx, t, side, "vertices").mean() * basin
        print(f"| {t} m | {over.sum() * side / nodes:.2%} | {over.mean():,.1f} ({lo:,.1f}-{hi:,.1f}) "
              f"| {over.mean() * basin / 1e6:,.2f} M | {over.mean() * basin / verts:.0%} |")
    relief = np.array([b["z_max"] - b["z_min"] for b in bx])
    print(f"\nRelief inside a box (z max - z min of its grid): median {np.median(relief):.0f} m, "
          f"p10 {np.quantile(relief, 0.1):.0f} m, p90 {np.quantile(relief, 0.9):.0f} m.\n")


def validation_table(path: str, piece: dict, res: dict) -> None:
    d, bx = boxes(path)
    side = (d["side_m"] / 1000.0) ** 2
    rng = np.random.default_rng(1)
    print(f"## Check of the box method: {len(bx)} random {d['side_m'] / 1000:g} km boxes inside "
          "the piece against the piece meshed whole\n")
    print("| tolerance | piece, whole catchment, triangles / km² | boxes inside it, mean (95 %) | "
          "inside the interval |\n|---:|---:|---:|---|")
    for t in bx[0]["tolerances"]:
        x = density(bx, t, side)
        lo, hi = interval(x, rng)
        whole = next(r["tri"] for r in piece["rows"] if r["tol"] == float(t)) / piece["area"]
        print(f"| {t} m | {whole:,.1f} | {x.mean():,.1f} ({lo:,.1f}-{hi:,.1f}) | "
              f"{lo <= whole <= hi} |")
    relief = np.array([b["z_max"] - b["z_min"] for b in bx])
    print(f"\nRelief inside a box: median {np.median(relief):.0f} m.\n")


def scaling_table(path: str) -> None:
    rec = json.loads(Path(path).read_text())
    ts = sorted({s["threads"] for s in rec["samples"]})
    one = med([s["refine_s"] for s in rec["samples"] if s["threads"] == 1])
    print(f"## Refine thread scaling on the piece at {rec['tolerance']:g} m\n")
    print("| threads | refine (median) | speed-up | process wall (median) |\n|---:|---:|---:|---:|")
    for t in ts:
        m = med([s["refine_s"] for s in rec["samples"] if s["threads"] == t])
        p = med([s["proc_s"] for s in rec["samples"] if s["threads"] == t])
        print(f"| {t} | {m:.2f} s | {one / m:.2f} | {p:.2f} s |")
    ac = "AC Power" in rec["pmset_before"] and "AC Power" in rec["pmset_after"]
    print(f"\nPower: AC before and after: {ac}.\n")


def main() -> None:
    res = json.loads(Path(sys.argv[1]).read_text())
    a, basin = area_km2(sys.argv[2]), area_km2(sys.argv[3])
    print(f"Basin (BHO level 2, ottobasin 76): {basin:,.1f} km² geodesic; OAS: {OAS_KM2:,.0f} km².\n")
    piece = piece_table(res, a)
    basin_table(sys.argv[4], basin, piece)
    validation_table(sys.argv[5], piece, res)
    if len(sys.argv) > 6:
        scaling_table(sys.argv[6])


if __name__ == "__main__":
    main()
