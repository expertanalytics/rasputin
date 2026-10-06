"""Pool the base and head runs of one prefix: medians per (domain, threads), head vs base.

    python summarize.py <runs-dir> <prefix>     e.g. . pra    or . vel

Reads <runs-dir>/<prefix>-base-r*/run.json and <prefix>-head-r*/run.json.
Prints a Markdown table of pooled refine_s, app_s and proc_s medians and the
head/base change, the per-run range of the refine change, and the mesh hashes.
"""

import json
import statistics
import sys
from pathlib import Path

root, prefix = Path(sys.argv[1]), sys.argv[2]


def load(side: str) -> list[dict]:
    return [json.loads(p.read_text()) for p in sorted(root.glob(f"{prefix}-{side}-r[1-9]*/run.json"))]


runs = {s: load(s) for s in ("base", "head")}
for s, rs in runs.items():
    print(f"{s}: {len(rs)} runs, commits {sorted({r['tree']['commit'][:7] for r in rs})}, "
          f"power {sorted({r['power']['state'] for r in rs})}, hardening {sorted({r['hardening'] for r in rs})}")
print()
print("| domain | mesh sha256 (base) | mesh sha256 (head) | equal |")
print("|---|---|---|---|")
doms = list(runs["base"][0]["quality"])
for d in doms:
    b = {r["quality"][d]["mesh_sha256"] for r in runs["base"]}
    h = {r["quality"][d]["mesh_sha256"] for r in runs["head"]}
    print(f"| {d} | {', '.join(sorted(b))} | {', '.join(sorted(h))} | {b == h and len(b) == 1} |")
print()


def pooled(side: str, d: str, t: int, key: str) -> float:
    return statistics.median(s[key] for r in runs[side] for s in r["samples"] if s["domain"] == d and s["threads"] == t)


def per_run(side: str, d: str, t: int, key: str) -> list[float]:
    return [statistics.median(s[key] for s in r["samples"] if s["domain"] == d and s["threads"] == t) for r in runs[side]]


threads = sorted({s["threads"] for s in runs["base"][0]["samples"]})
print("| domain | threads | refine_s base | refine_s head | change | app_s change | proc_s change |")
print("|---|---:|---:|---:|---:|---:|---:|")
allch = []
for d in doms:
    for t in threads:
        b, h = pooled("base", d, t, "refine_s"), pooled("head", d, t, "refine_s")
        ch = 100 * (h / b - 1)
        allch.append(ch)
        a = 100 * (pooled("head", d, t, "app_s") / pooled("base", d, t, "app_s") - 1)
        p = 100 * (pooled("head", d, t, "proc_s") / pooled("base", d, t, "proc_s") - 1)
        print(f"| {d} | {t or 'default'} | {b:.4f} | {h:.4f} | {ch:+.1f} % | {a:+.1f} % | {p:+.1f} % |")
print()
print(f"refine_s pooled change over every cell: {min(allch):+.1f} .. {max(allch):+.1f} %")
bb = [x for d in doms for t in threads for x in per_run("base", d, t, "refine_s")]
spread = []
for d in doms:
    for t in threads:
        v = per_run("base", d, t, "refine_s")
        spread.append(100 * (max(v) / min(v) - 1))
print(f"base-vs-base noise (per cell, max/min of the base runs' medians): median {statistics.median(spread):.1f} %, max {max(spread):.1f} %")
