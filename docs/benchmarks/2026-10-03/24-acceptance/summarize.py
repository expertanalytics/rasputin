"""Tables for 24-acceptance.md from the six bench.py run.json files (@perf).

Run from anywhere: python summarize.py [DIR] > tables.md. Pairs, in run order
M OFF ON ON2 OFF2 M2: "nothing else moved" is (M, OFF) and (M2, OFF2); the
hardening cost is (OFF, ON) and (OFF2, ON2). "pooled" is the median over both
runs' samples of a side (10 per thread count) divided likewise.
"""

import json
import statistics
import sys
from pathlib import Path

# The directory holding the six run directories: argv[1], default the date folder.
HERE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent
RUNS = {
    "M": "24-base-6cdc8cc", "OFF": "24-off", "ON": "24-on",
    "ON2": "24-on-r2", "OFF2": "24-off-r2", "M2": "24-base-6cdc8cc-r2",
}  # fmt: skip
rec = {k: json.loads((HERE / v / "run.json").read_text()) for k, v in RUNS.items()}
THREADS = sorted({s["threads"] for s in rec["M"]["stats"]})


def samples(k: str, dom: str, t: int, key: str = "refine_s") -> list[float]:
    return [s[key] for s in rec[k]["samples"] if s["domain"] == dom and s["threads"] == t]


def med(k: str, dom: str, t: int, key: str = "refine_s") -> float:
    return statistics.median(samples(k, dom, t, key))


def pooled(ks: tuple[str, str], dom: str, t: int, key: str = "refine_s") -> float:
    return statistics.median(samples(ks[0], dom, t, key) + samples(ks[1], dom, t, key))


def pct(x: float) -> str:
    return f"{(x - 1) * 100:+.1f}"


def table(title: str, base: tuple[str, str], new: tuple[str, str], key: str = "refine_s") -> dict:
    """Per-pair and pooled change of `new` over `base`, percent; returns pooled per (dom, t)."""
    b1, b2 = base
    n1, n2 = new
    print(f"### {title} ({key}, median of 5 per run, seconds; change in %)\n")
    print(f"| domain | threads | {b1} | {n1} | {n2} | {b2} | {n1}/{b1} | {n2}/{b2} | pooled |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    out = {}
    for dom in ("tile", "quarter"):
        for t in THREADS:
            a, b, c, d = (med(k, dom, t, key) for k in (b1, n1, n2, b2))
            p = pooled(new, dom, t, key) / pooled(base, dom, t, key)
            out[(dom, t)] = (b / a, c / d, p)
            tl = "default" if t == 0 else str(t)
            print(f"| {dom} | {tl} | {a:.4f} | {b:.4f} | {c:.4f} | {d:.4f} | {pct(b / a)} | {pct(c / d)} | {pct(p)} |")
    print()
    return out


moved = table("Nothing else moved: OFF against the previous merge", ("M", "M2"), ("OFF", "OFF2"))
cost = table("Hardening cost: ON against OFF", ("OFF", "OFF2"), ("ON", "ON2"))
table("Hardening cost, in process (read, mesh, write)", ("OFF", "OFF2"), ("ON", "ON2"), "app_s")
table("Hardening cost, whole child process", ("OFF", "OFF2"), ("ON", "ON2"), "proc_s")

print("### Summary\n")
for name, res in (("OFF vs M", moved), ("ON vs OFF", cost)):
    for dom in ("tile", "quarter"):
        cells = {t: v for (d, t), v in res.items() if d == dom}
        per_pair = [x for v in cells.values() for x in v[:2]]
        multi = [cells[t][2] for t in THREADS if t >= 2]
        print(f"- {name}, {dom}: per-pair {pct(min(per_pair))}..{pct(max(per_pair))} %; pooled "
              f"default {pct(cells[0][2])}, 1 thread {pct(cells[1][2])}, 2-20 threads "
              f"{pct(min(multi))}..{pct(max(multi))} %")
for side, ks in (("OFF", ("OFF", "OFF2")), ("ON", ("ON", "ON2")), ("M", ("M", "M2"))):
    for dom in ("tile", "quarter"):
        m = {t: pooled(ks, dom, t) for t in THREADS if t >= 1}
        best = min(m, key=m.__getitem__)
        print(f"- ceiling {side} {dom}: 1 thread over best {m[1] / m[best]:.2f}x (best at {best}), "
              f"over 20 {m[1] / m[20]:.2f}x")
print()

print("### Builds and geometry\n")
print("| run | commit | dirty | hardening | power | `_core` sha256 | tile mesh sha256 | quarter mesh sha256 |")
print("|---|---|---|---|---|---|---|---|")
for k, lab in RUNS.items():
    r = rec[k]
    q = r["quality"]
    print(f"| {lab} | {r['tree']['commit'][:7]} | {r['tree']['dirty']} | {r.get('hardening')} | "
          f"{r['power']['state']} {r['power']['percent']} % | {r['build']['so_sha256'][:12]} | "
          f"{q['tile']['mesh_sha256'][:12]} | {q['quarter']['mesh_sha256'][:12]} |")
print()

print("### Quality\n")
print("| run | domain | worst angle deg | max degree | within tol | Delaunay viol. | max_error | rounds | inserted | flips |")
print("|---|---|---:|---:|---|---:|---:|---:|---:|---:|")
tuples: dict[tuple[str, int], set] = {}
for k, lab in RUNS.items():
    r = rec[k]
    for dom in ("tile", "quarter"):
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
              f"{q['delaunay_violations']} of {q['delaunay_checked']:,} | "
              f"{me:.7f} | {ro} | {ins} | {fl} |")
bad = {key: v for key, v in tuples.items() if len(v) != 1}
print(f"\nCounter tuples (max_error, rounds, inserted, flips) per (domain, threads) over all six runs: "
      f"{len(tuples)} cells, {len(bad)} with more than one value.")
