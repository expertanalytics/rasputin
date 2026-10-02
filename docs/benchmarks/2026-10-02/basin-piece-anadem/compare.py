"""ANADEM against GLO-30, side by side: triangles per tolerance and the final check's excess.

A measurement script, not production code; nothing imports it. Reads the
JSON of both runs (the GLO-30 one in ``docs/benchmarks/2026-10-01/basin-piece/``)
and prints Markdown. The two box samples use the same seed, so they are the
same 200 boxes; the per-box ratio is paired.

    python docs/benchmarks/2026-10-02/basin-piece-anadem/compare.py GLO_DIR ANADEM_DIR BASIN_OUTLINE

GLO_DIR and ANADEM_DIR are the ``runs/`` folders of the two measurements.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from analyse import area_km2  # noqa: E402


def load(d: Path, src: str) -> tuple[dict, list[dict], list[dict]]:
    piece = json.loads((d / f"velhas76949_{src}" / "results.json").read_text())
    basin = json.loads((d / f"basin_boxes_{src}" / "boxes.json").read_text())["boxes"]
    pboxes = json.loads((d / f"piece_boxes_{src}" / "boxes.json").read_text())["boxes"]
    return piece, basin, pboxes


def main() -> None:
    gd, ad, basin_outline = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    gp, gb, gpb = load(gd, "glo30")
    ap, ab, apb = load(ad, "anadem")
    basin = area_km2(basin_outline)
    gb_by = {b["box"]: b for b in gb}
    assert all(abs(b["lon"] - gb_by[b["box"]]["lon"]) < 1e-12 for b in ab), "box samples differ"
    pair = [(gb_by[b["box"]], b) for b in ab]
    print(f"Paired boxes: {len(pair)} (same seed, same centres).\n")
    print("| tolerance | piece triangles, GLO-30 | piece triangles, ANADEM | ANADEM / GLO-30 | "
          "basin triangles, GLO-30 | basin triangles, ANADEM | ANADEM / GLO-30 | "
          "paired box ratio, median (p10-p90) |\n|---:|---:|---:|---:|---:|---:|---:|---:|")
    for g, a in zip(gp["runs"], ap["runs"]):
        t = f"{g['tolerance']:g}"
        assert g["tolerance"] == a["tolerance"]
        tg, ta = g["timed"][0]["triangles"], a["timed"][0]["triangles"]
        bg = np.mean([x["tolerances"][t]["triangles"] for x, _ in pair]) / 100 * basin
        ba = np.mean([y["tolerances"][t]["triangles"] for _, y in pair]) / 100 * basin
        r = np.array([y["tolerances"][t]["triangles"] / max(x["tolerances"][t]["triangles"], 1)
                      for x, y in pair])  # fmt: skip
        print(f"| {t} m | {tg:,} | {ta:,} | {ta / tg:.3f} | {bg / 1e6:,.1f} M | {ba / 1e6:,.1f} M | "
              f"{ba / bg:.3f} | {np.median(r):.2f} ({np.quantile(r, 0.1):.2f}-{np.quantile(r, 0.9):.2f}) |")
    print("\nThe final check's excess: source-DEM nodes (each DEM's own) off the mesh by more "
          "than the tolerance, interior only. Piece: share of its interior source nodes. "
          "Basin: the boxes' mean per km² times the basin's area, and as a share of the "
          "basin's vertices.\n")
    print("| tolerance | piece, GLO-30 | piece, ANADEM | basin nodes, GLO-30 | basin nodes, ANADEM "
          "| share of vertices, GLO-30 | share of vertices, ANADEM |\n"
          "|---:|---:|---:|---:|---:|---:|---:|")
    for g, a in zip(gp["runs"], ap["runs"]):
        t = f"{g['tolerance']:g}"
        cells = []
        for k in (0, 1):
            fc = (g, a)[k]["final_check"]["source_dem_interior"]
            cells.append(f"{fc['over_share']:.2%} (max {fc['max']:.1f} m)")
        bas = []
        for k in (0, 1):
            over = np.mean([p[k]["tolerances"][t]["final_check"]["source_dem_interior"]["over_tolerance"]
                            for p in pair]) / 100 * basin  # fmt: skip
            verts = np.mean([p[k]["tolerances"][t]["vertices"] for p in pair]) / 100 * basin
            bas.append((over, verts))
        print(f"| {t} m | {cells[0]} | {cells[1]} | {bas[0][0] / 1e6:,.2f} M | {bas[1][0] / 1e6:,.2f} M "
              f"| {bas[0][0] / bas[0][1]:.0%} | {bas[1][0] / bas[1][1]:.0%} |")
    print("\nThe boundary strip on the piece (within 30 m of the outline), worst source node:\n")
    print("| tolerance | GLO-30 | ANADEM |\n|---:|---:|---:|")
    for g, a in zip(gp["runs"], ap["runs"]):
        s = [x["final_check"]["source_dem_strip"] for x in (g, a)]
        print(f"| {g['tolerance']:g} m | {s[0]['max']:.1f} m ({s[0]['over_share']:.1%} over) "
              f"| {s[1]['max']:.1f} m ({s[1]['over_share']:.1%} over) |")
    rg = np.array([x["z_max"] - x["z_min"] for x, _ in pair])
    ra = np.array([y["z_max"] - y["z_min"] for _, y in pair])
    print(f"\nRelief inside a box, median: GLO-30 {np.median(rg):.0f} m, ANADEM {np.median(ra):.0f} m.")


if __name__ == "__main__":
    main()
