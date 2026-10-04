"""Tables for 23b-acceptance.md from the four bench.py run.json files (@perf).

Run from anywhere: python summarize.py [DIR] > tables.md. Run order B N N2 B2:
the pairs are (B, N) and (B2, N2); "pooled" is the median over both runs'
samples of a side (10 per thread count) divided likewise. "same-binary" is
B2 against B and N2 against N: the session's drift, i.e. the noise.
"""

import json
import statistics
import sys
from pathlib import Path

HERE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent
RUNS = {"B": "23b-base-b4bcdc3", "N": "23b", "N2": "23b-r2", "B2": "23b-base-b4bcdc3-r2"}
rec = {k: json.loads((HERE / v / "run.json").read_text()) for k, v in RUNS.items()}
THREADS = sorted({s["threads"] for s in rec["B"]["samples"]})
DOMS = sorted({s["domain"] for s in rec["B"]["samples"]}, key=lambda d: d != "tile")


def samples(k, dom, t, key="refine_s"):
    return [s[key] for s in rec[k]["samples"] if s["domain"] == dom and s["threads"] == t]


def med(k, dom, t, key="refine_s"):
    return statistics.median(samples(k, dom, t, key))


def pooled(ks, dom, t, key="refine_s"):
    return statistics.median(samples(ks[0], dom, t, key) + samples(ks[1], dom, t, key))


def pct(x):
    return f"{(x - 1) * 100:+.1f}"


def table(key="refine_s"):
    print(f"### 23b against the merge base ({key}, median of 5 per run, seconds; change in %)\n")
    print("| domain | threads | B | N | N2 | B2 | N/B | N2/B2 | pooled | B2/B | N2/N |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    out = {}
    for dom in DOMS:
        for t in THREADS:
            b, n, n2, b2 = (med(k, dom, t, key) for k in ("B", "N", "N2", "B2"))
            p = pooled(("N", "N2"), dom, t, key) / pooled(("B", "B2"), dom, t, key)
            out[(dom, t)] = (n / b, n2 / b2, p, b2 / b, n2 / n)
            tl = "default" if t == 0 else str(t)
            print(f"| {dom} | {tl} | {b:.4f} | {n:.4f} | {n2:.4f} | {b2:.4f} | {pct(n / b)} | "
                  f"{pct(n2 / b2)} | {pct(p)} | {pct(b2 / b)} | {pct(n2 / n)} |")
    print()
    return out


res = table()
table("app_s")
table("proc_s")
print("### Summary\n")
for dom in DOMS:
    cells = {t: v for (d, t), v in res.items() if d == dom}
    pair = [x for v in cells.values() for x in v[:2]]
    pool = [v[2] for v in cells.values()]
    same = [x for v in cells.values() for x in v[3:]]
    print(f"- {dom}: per-pair {pct(min(pair))}..{pct(max(pair))} %; pooled {pct(min(pool))}..{pct(max(pool))} %; "
          f"same-binary drift {pct(min(same))}..{pct(max(same))} %")
for side, ks in (("base", ("B", "B2")), ("23b", ("N", "N2"))):
    for dom in DOMS:
        m = {t: pooled(ks, dom, t) for t in THREADS if t >= 1}
        best = min(m, key=m.__getitem__)
        print(f"- ceiling {side} {dom}: 1 thread over best {m[1] / m[best]:.2f}x (best at {best}), "
              f"over {max(m)} {m[1] / m[max(m)]:.2f}x")
print()
print("### Builds and geometry\n")
print("| run | commit | dirty | hardening | power before/after | `_core` sha256 | "
      + " | ".join(f"{d} mesh sha256" for d in DOMS) + " |")
print("|---|---|---|---|---|---|" + "---|" * len(DOMS))
for k, lab in RUNS.items():
    r = rec[k]
    q = r["quality"]
    print(f"| {lab} | {r['tree']['commit'][:7]} | {r['tree']['dirty']} | {r.get('hardening')} | "
          f"{r['power']['state']} {r['power']['percent']} % | {r['build']['so_sha256'][:12]} | "
          + " | ".join(f"`{q[d]['mesh_sha256'][:16]}`" for d in DOMS) + " |")
print()
print("### Quality\n")
print("| run | domain | worst angle deg | max degree | within tol | Delaunay viol. | max_error | rounds | inserted | flips |")
print("|---|---|---:|---:|---|---:|---:|---:|---:|---:|")
tuples = {}
for k, lab in RUNS.items():
    r = rec[k]
    for dom in DOMS:
        q = r["quality"][dom]
        for s in r["samples"]:
            if s["domain"] == dom:
                tuples.setdefault((dom, s["threads"]), set()).add(
                    (s["max_error"], s["rounds"], s["inserted"], s["flips"]))
        c0 = {(s["max_error"], s["rounds"], s["inserted"], s["flips"])
              for s in r["samples"] if s["domain"] == dom and s["threads"] == 0}
        assert len(c0) == 1, (lab, dom, c0)
        me, ro, ins, fl = c0.pop()
        print(f"| {lab} | {dom} | {q['worst_angle']:.4f} | {q['max_degree']} | {q['within_tolerance']} | "
              f"{q['delaunay_violations']} of {q['delaunay_checked']:,} | {me:.7f} | {ro} | {ins} | {fl} |")
bad = {key: v for key, v in tuples.items() if len(v) != 1}
print(f"\nCounter tuples (max_error, rounds, inserted, flips) per (domain, threads) over all four runs: "
      f"{len(tuples)} cells, {len(bad)} with more than one value.")
