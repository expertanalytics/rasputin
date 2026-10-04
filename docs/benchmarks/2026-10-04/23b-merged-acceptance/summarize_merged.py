"""23b at the merged head 4541e38 (F) against master d20126b (B), per side (@perf).
python summarize_merged.py [DIR] > tables-merged.md. Derived from 23b-acceptance/summarize_fix2.py.

Sides: B (master d20126b), F (23b merged head 4541e38). Order: B F F B, then F B B F.
"pooled" is the median over all samples of a side's runs (5 per run).
"""
import json
import statistics as st
import sys
from pathlib import Path

HERE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent
SIDES = {
    "B": ["23bm-base-r1", "23bm-base-r2", "23bm-base-r3", "23bm-base-r4"],
    "F": ["23bm-r1", "23bm-r2", "23bm-r3", "23bm-r4"],
}
EXTRA = []
rec = {k: json.loads((HERE / k / "run.json").read_text()) for v in SIDES.values() for k in v + EXTRA if (HERE / k).exists()}
THREADS = sorted({s["threads"] for s in rec["23bm-r1"]["samples"]})


def smp(k, d, t, key="refine_s"):
    return [s[key] for s in rec[k]["samples"] if s["domain"] == d and s["threads"] == t]


def pool(side, d, t, key="refine_s"):
    return st.median([x for k in SIDES[side] for x in smp(k, d, t, key)])


def rng(side, d, t):
    m = [st.median(smp(k, d, t)) for k in SIDES[side]]
    return min(m), max(m)


pc = lambda x: f"{(x - 1) * 100:+.1f}"
print("### refine_s per side (pooled median, s; run-median range; change against B in %)\n")
print("| domain | threads | B pooled | B runs | F pooled | F runs | F/B |")
print("|---|---:|---:|---|---:|---|---:|")
for d in ("tile", "quarter"):
    for t in THREADS:
        b, n = (pool(s, d, t) for s in "BF")
        r = {s: rng(s, d, t) for s in "BF"}
        f = lambda s: f"{r[s][0]:.4f}-{r[s][1]:.4f}"
        print(f"| {d} | {'default' if t == 0 else t} | {b:.4f} | {f('B')} | {n:.4f} | {f('F')} | {pc(n / b)} |")
print()
print("### Every run\n")
print("| run | commit | dirty | `_core` sha256 | power before / after | start (UTC) | tile mesh | quarter mesh |")
print("|---|---|---|---|---|---|---|---|")
for k, r in rec.items():
    q = r["quality"]
    pw = r["power"]
    print(f"| {k} | {r['tree']['commit'][:7]} | {r['tree']['dirty']} | {r['build']['so_sha256'][:12]} | "
          f"{pw['state']} {pw['percent']} % | {r.get('started', r.get('date', ''))} | "
          f"`{q['tile']['mesh_sha256'][:16]}` | `{q['quarter']['mesh_sha256'][:16]}` |")
tup = {}
for k, r in rec.items():
    for s in r["samples"]:
        tup.setdefault((s["domain"], s["threads"]), set()).add((s["max_error"], s["rounds"], s["inserted"], s["flips"]))
print(f"\nCounter tuples (max_error, rounds, inserted, flips) per (domain, threads) over all {len(rec)} runs: "
      f"{len(tup)} cells, {sum(len(v) != 1 for v in tup.values())} with more than one value.")
for d in ("tile", "quarter"):
    hs = {r["quality"][d]["mesh_sha256"] for r in rec.values()}
    qs = {(r["quality"][d]["worst_angle"], r["quality"][d]["max_degree"], r["quality"][d]["within_tolerance"],
           r["quality"][d]["delaunay_violations"]) for r in rec.values()}
    print(f"- {d}: {len(hs)} distinct mesh sha256 ({', '.join(hs)}); quality tuples {qs}")
