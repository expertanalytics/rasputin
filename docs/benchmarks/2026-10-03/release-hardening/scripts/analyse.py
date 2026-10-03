"""Pairwise overhead of a hardened Release over plain Release ("fast" below is
the hardened side: libc++ FAST mode, or libstdc++ _GLIBCXX_ASSERTIONS).

Usage: analyse.py ROOT PAIR... where PAIR is plainlabel:hardenedlabel.
Reads each bench.py run.json under ROOT/*/<label>/ (here: ROOT=. and raw/), prints:
 - per run: power, .so sha256, mesh sha256 per domain, quality;
 - per (domain, threads, measure): each pair's median ratio fast/plain - 1,
   and the pooled-sample median ratio over all pairs.
"""

import json
import statistics as st
import sys
from pathlib import Path

root = Path(sys.argv[1])
pairs = [p.split(":") for p in sys.argv[2:]]
runs = {p.parent.name: json.loads(p.read_text()) for p in root.glob("*/*/run.json")}

print("## runs\n")
print("| run | started | power | _core sha256 | tile mesh sha256 | quarter mesh sha256 |")
print("|---|---|---|---|---|---|")
for label in [x for pair in pairs for x in pair]:
    r = runs[label]
    q = r["quality"]
    print(f"| {label} | {r['started'][11:19]} | {r['power']['state']} {r['power'].get('percent')}% | "
          f"`{r['build']['so_sha256'][:12]}` | `{q['tile']['mesh_sha256'][:12]}` | "
          f"`{q['quarter']['mesh_sha256'][:12]}` |")

print("\n## quality (identical across runs?)\n")
keys = ["worst_angle", "max_degree", "within_tolerance", "delaunay_violations", "mesh_sha256"]
for dom in ("tile", "quarter"):
    seen = {json.dumps({k: runs[l]["quality"][dom][k] for k in keys}) for p in pairs for l in p}
    print(f"- {dom}: {len(seen)} distinct; {sorted(seen)[0]}")
counters = {(s["domain"], s["threads"], s["rounds"], s["inserted"], s["flips"], s["max_error"])
            for p in pairs for l in p for s in runs[l]["samples"]}
print(f"- distinct (domain, threads, rounds, inserted, flips, max_error) tuples: {len(counters)}")


def samples(label, dom, t, m):
    return [s[m] for s in runs[label]["samples"] if s["domain"] == dom and s["threads"] == t]


threads = sorted({s["threads"] for s in runs[pairs[0][0]]["samples"]})


def pooled(dom, t, m):
    pp = st.median([x for p, _ in pairs for x in samples(p, dom, t, m)])
    ff = st.median([x for _, f in pairs for x in samples(f, dom, t, m)])
    return 100 * (ff / pp - 1)


def ceiling(side, dom):
    med = {t: st.median([x for pair in pairs for x in samples(pair[side], dom, t, "refine_s")])
           for t in threads if t}
    best = min(med, key=med.get)
    return f"{med[1] / med[best]:.2f}x at {best}"


print("\n## summary: pooled overhead %, t=1 | t=0 (CLI default) | min..max over t=2..20\n")
for m in ("refine_s", "app_s", "proc_s"):
    for dom in ("tile", "quarter"):
        rest = [pooled(dom, t, m) for t in threads if t >= 2]
        print(f"- {m} {dom}: {pooled(dom, 1, m):+.1f} | {pooled(dom, 0, m):+.1f} | "
              f"{min(rest):+.1f}..{max(rest):+.1f}")
for dom in ("tile", "quarter"):
    print(f"- refine ceiling {dom} (pooled medians): plain {ceiling(0, dom)}, "
          f"hardened {ceiling(1, dom)}")

for m in ("refine_s", "app_s", "proc_s"):
    print(f"\n## {m}: overhead fast/plain - 1, % (per pair; pooled medians)\n")
    print("| domain | threads | plain median s | fast median s | " +
          " | ".join(f"pair {i + 1}" for i in range(len(pairs))) + " | pooled | spread (min..max pair) |")
    print("|---|---:|---:|---:|" + "---:|" * (len(pairs) + 2))
    for dom in ("tile", "quarter"):
        for t in threads:
            per = [st.median(samples(f, dom, t, m)) / st.median(samples(p, dom, t, m)) - 1
                   for p, f in pairs]
            pp = st.median([x for p, _ in pairs for x in samples(p, dom, t, m)])
            ff = st.median([x for _, f in pairs for x in samples(f, dom, t, m)])
            print(f"| {dom} | {t} | {pp:.4f} | {ff:.4f} | " +
                  " | ".join(f"{100 * x:+.1f}" for x in per) +
                  f" | {100 * (ff / pp - 1):+.1f} | {100 * min(per):+.1f}..{100 * max(per):+.1f} |")
