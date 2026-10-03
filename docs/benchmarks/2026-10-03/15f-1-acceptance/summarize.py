"""Tables for 15f-1-acceptance.md from the four bench.py run.json files (@perf)."""

import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
RUNS = {"A": "15f-1-base-390b516", "B": "15f-1", "B2": "15f-1-r2", "A2": "15f-1-base-390b516-r2"}
rec = {k: json.loads((HERE / v / "run.json").read_text()) for k, v in RUNS.items()}


def med(k: str, dom: str, t: int) -> float:
    return next(s["median"] for s in rec[k]["stats"] if s["domain"] == dom and s["threads"] == t)


print("| domain | threads | A base | B 15f-1 | B2 15f-1 | A2 base | B/A | B2/A2 | (B+B2)/(A+A2) |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
for dom in ("quarter", "tile"):
    for t in sorted({s["threads"] for s in rec["A"]["stats"]}):
        a, b, b2, a2 = (med(k, dom, t) for k in ("A", "B", "B2", "A2"))
        tl = "default" if t == 0 else str(t)
        print(f"| {dom} | {tl} | {a:.4f} | {b:.4f} | {b2:.4f} | {a2:.4f} | {b / a:.3f} | {b2 / a2:.3f} | {(b + b2) / (a + a2):.3f} |")
print()
print("| run | commit | power | .so sha256 | tile proc_s | quarter proc_s | tile mesh sha256 | quarter mesh sha256 |")
print("|---|---|---|---|---:|---:|---|---|")
for k, lab in RUNS.items():
    r = rec[k]

    def proc(dom: str) -> float:
        return statistics.median(s["proc_s"] for s in r["samples"] if s["domain"] == dom and s["threads"] == 0)

    q = r["quality"]
    print(f"| {lab} | {r['tree']['commit'][:7]} | {r['power']['state']} | {r['build']['so_sha256'][:12]} | "
          f"{proc('tile'):.3f} | {proc('quarter'):.3f} | {q['tile']['mesh_sha256'][:12]} | {q['quarter']['mesh_sha256'][:12]} |")
print()
print("| run | domain | worst angle deg | share < 1 deg | max degree | within tol | Delaunay viol. | max_error | rounds | inserted | flips |")
print("|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|")
for k, lab in RUNS.items():
    r = rec[k]
    for dom in ("tile", "quarter"):
        q = r["quality"][dom]
        c = {(s["max_error"], s["rounds"], s["inserted"], s["flips"]) for s in r["samples"] if s["domain"] == dom}
        assert len(c) == 1, (lab, dom, c)  # every sample of a domain, every thread count, agrees
        me, ro, ins, fl = c.pop()
        print(f"| {lab} | {dom} | {q['worst_angle']:.4f} | {q['share_under_1']:.2e} | {q['max_degree']} | "
              f"{q['within_tolerance']} | {q['delaunay_violations']} | {me:.7f} | {ro} | {ins} | {fl} |")
