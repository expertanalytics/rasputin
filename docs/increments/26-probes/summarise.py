"""Throwaway: print fractions_probe.py's JSON as compact Markdown tables.

python docs/increments/26-probes/summarise.py RUN.json [RUN.json ...]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def fmt(x: list[float] | None) -> str:
    return "-" if x is None else f"{x[1]} / {x[2]}"


for path in sys.argv[1:]:
    r = json.loads(Path(path).read_text())
    q = r["triangle_m2_p50_p90_p99_max"]
    print(
        f"\n### {r['mesh']}: {r['triangles']:,} triangles, "
        f"mean {r['mean_triangle_cells']:.2f} cells ({r['mean_triangle_m2']:,.0f} m2); "
        f"triangle area p50/p90/p99/max "
        f"{q[0]:,.0f} / {q[1]:,.0f} / {q[2]:,.0f} / {q[3]:,.0f} m2; "
        f"exact classes per triangle mean {r['exact_nnz_mean']:.2f}, "
        f"max {r['exact_nnz_max']}; overlap {r['seconds_overlap']} s"
    )
    print(
        f"median triangle, dominant crop {r['median_triangle_m2_dominant_crop']:,.0f} m2 "
        f"({r['triangles_dominant_crop']:,} triangles); dominant other "
        f"{r['median_triangle_m2_dominant_other']:,.0f} m2"
    )
    print(
        "\n| cutoff | floor m2 | rule | classes/tri mean (max) | worst class % | water % "
        "| soybean % | cotton % | coffee % | misplaced 1 km p95/max % | 5 km p95/max % "
        "| 25 km p95/max % | ledger max (mean tri) |"
    )
    print("|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for v in r["variants"]:
        ce = v["crop_errors_pct"]
        g = v["regional_misplaced_pct_p50_p95_max_n"]
        worst = v["worst_class_error"]
        print(
            f"| {v['cutoff']:.0%} | {v['min_m2']:g} | {v['rule']} | "
            f"{v['nnz_mean']:.2f} ({v['nnz_max']}) | {worst[2]} ({worst[0]}) | "
            f"{v['water_error_pct']} | {ce.get('soybean')} | {ce.get('cotton')} | "
            f"{ce.get('coffee')} | {fmt(g['1km'])} | {fmt(g['5km'])} | {fmt(g['25km'])} | "
            f"{v['ledger_max_in_mean_triangles']} |"
        )
