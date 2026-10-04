"""Tables for 20-nodata-acceptance.md, from the eight bench.py runs (@perf).
python summarize_bench.py [DIR] > tables-bench.md. B = master 2060f14, N = the fix.
Order B N N B, N B B N. "pooled": median over all samples of a side (20 per cell)."""
import json
import statistics as st
import sys
from pathlib import Path

HERE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent
SIDES = {"B": [f"20nd-base-r{i}" for i in range(1, 5)], "N": [f"20nd-r{i}" for i in range(1, 5)]}
rec = {k: json.loads((HERE / k / "run.json").read_text()) for v in SIDES.values() for k in v}
THREADS = sorted({s["threads"] for s in rec["20nd-r1"]["samples"]})
pc = lambda x: f"{(x - 1) * 100:+.1f}"


def smp(k, d, t, key):
    return [s[key] for s in rec[k]["samples"] if s["domain"] == d and s["threads"] == t]


def pool(side, d, t, key):
    return st.median([x for k in SIDES[side] for x in smp(k, d, t, key)])


def rng(side, d, t, key):
    m = [st.median(smp(k, d, t, key)) for k in SIDES[side]]
    return f"{min(m):.4f}-{max(m):.4f}"


for key in ("refine_s", "app_s", "proc_s"):
    print(f"### {key} (pooled median, s; run-median range; change)\n")
    print("| domain | threads | B | B runs | N | N runs | N/B |")
    print("|---|---:|---:|---|---:|---|---:|")
    for d in ("tile", "quarter"):
        for t in THREADS:
            b, n = pool("B", d, t, key), pool("N", d, t, key)
            print(f"| {d} | {'default' if t == 0 else t} | {b:.4f} | {rng('B', d, t, key)} | {n:.4f} | "
                  f"{rng('N', d, t, key)} | {pc(n / b)} |")
    print()
print("### Every run\n")
print("| run | commit | dirty | `_core` | power | tile mesh | quarter mesh |")
print("|---|---|---|---|---|---|---|")
for k, r in rec.items():
    q = r["quality"]
    print(f"| {k} | {r['tree']['commit'][:7]} | {r['tree']['dirty']} | {r['build']['so_sha256'][:12]} | "
          f"{r['power']['state']} {r['power']['percent']} % | `{q['tile']['mesh_sha256'][:16]}` | "
          f"`{q['quarter']['mesh_sha256'][:16]}` |")
print("\n### Quality per side (identical across a side's runs?)\n")
for side, ks in SIDES.items():
    for d in ("tile", "quarter"):
        qs = {json.dumps({k: v for k, v in rec[k]["quality"][d].items()}, sort_keys=True) for k in ks}
        print(f"- {side} {d}: {len(qs)} distinct; {next(iter(qs))}")
        tup = {(s["max_error"], s["rounds"], s["inserted"], s["flips"]) for k in ks for s in rec[k]["samples"]
               if s["domain"] == d and s["threads"] == 0}
        print(f"  counters at default threads: {tup}")
