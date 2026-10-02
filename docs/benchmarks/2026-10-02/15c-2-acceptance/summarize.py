"""Tables for README.md from runs/geo/results.json and the four bench.py runs.

    python docs/benchmarks/2026-10-02/15c-2-acceptance/summarize.py > /tmp/tables.md
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
# The hand-resampled ANADEM run of 2026-10-02 (basin-piece-anadem/README.md).
HAND = {50.0: 74_841, 20.0: 264_833, 10.0: 649_845, 5.0: 1_576_681, 2.0: 4_584_784, 1.0: 8_699_303}


def bench_table() -> None:
    runs = ["15c-2-base-e1f6042", "15c-2", "15c-2-r2", "15c-2-base-e1f6042-r2"]
    rec = {n: json.loads((BENCH / n / "run.json").read_text()) for n in runs}
    st = {n: {(s["domain"], s["threads"]): s["median"] for s in r["stats"]} for n, r in rec.items()}
    print("| domain | threads | A base | B 15c-2 | B2 15c-2 | A2 base | B/A | B2/A2 | (B+B2)/(A+A2) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for k in sorted(st[runs[0]]):
        a, b, b2, a2 = (st[n][k] for n in runs)
        print(f"| {k[0]} | {k[1] or 'default'} | {a:.4f} | {b:.4f} | {b2:.4f} | {a2:.4f} | "
              f"{b / a:.3f} | {b2 / a2:.3f} | {(b + b2) / (a + a2):.3f} |")  # fmt: skip
    print()
    print("| run | power | .so sha256 | tile proc_s (default threads) | quarter proc_s | "
          "tile mesh sha256 | quarter mesh sha256 |")  # fmt: skip
    print("|---|---|---|---:|---:|---|---|")
    for n in runs:
        r = rec[n]
        p = lambda d: np.median([s["proc_s"] for s in r["samples"] if s["domain"] == d and s["threads"] == 0])  # noqa: E731
        print(f"| {n} | {r['power']['state']} | {r['build']['so_sha256'][:12]} | {p('tile'):.3f} | "
              f"{p('quarter'):.3f} | {r['quality']['tile']['mesh_sha256'][:12]} | "
              f"{r['quality']['quarter']['mesh_sha256'][:12]} |")  # fmt: skip


def geo_table() -> None:
    res = json.loads((HERE / "runs/geo/results.json").read_text())
    print("| out CRS | tol | triangles | vs hand-resampled | phase 2 inserted (rounds) | proc s | "
          "refine s | other s | peak RSS GiB | worst angle | max deg | CDT viol. | over tol: interior | strip |")  # fmt: skip
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in res["runs"]:
        t = r["timed"]
        med = lambda f: float(np.median([f(x) for x in t]))  # noqa: E731
        tri, tol = t[0]["triangles"], r["tolerance"]
        crs = "suggested tmerc" if r["crs"].startswith("+proj") else r["crs"]
        q, c = r["quality"], r.get("check")
        chk = (f"{c['interior']['over_tolerance']} of {c['interior']['nodes']:,} | "
               f"{c['strip']['over_tolerance']} of {c['strip']['nodes']:,}") if c else "not run | not run"  # fmt: skip
        p2 = r["phase2"]
        print(f"| {crs} | {tol:g} | {tri:,} | {tri / HAND[tol]:.3f} | {p2['inserted']:,} ({p2['rounds']}) | "
              f"{med(lambda x: x['proc_s']):.2f} | {med(lambda x: float(x['rows']['refine'])):.2f} | "
              f"{med(lambda x: float(x['rows']['other'])):.2f} | "
              f"{med(lambda x: x['max_rss_bytes']) / 2**30:.2f} | {q['worst_angle']:.4f} | "
              f"{q['max_degree']} | {q['delaunay_violations']} | {chk} |")  # fmt: skip
    print()
    print(f"Failed attempts (retried): {len(res.get('failures', []))}")
    for f in res.get("failures", []):
        print(f"- {f['utc']} {f['log']}: {f['kind']}")


if __name__ == "__main__":
    bench_table()
    print()
    geo_table()
