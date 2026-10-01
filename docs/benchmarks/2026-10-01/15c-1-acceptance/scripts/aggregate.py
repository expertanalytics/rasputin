"""Aggregate the 15c-1 acceptance pairs: per (domain, threads), base vs branch.

Run from the 15c-1 worktree root with any Python 3.12+:
    python3 docs/benchmarks/2026-10-01/15c-1-acceptance/scripts/aggregate.py
Reads every run.json in the acceptance directory; prints markdown tables.
"""

import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]


def load() -> dict[str, dict]:
    return {p.parent.name: json.loads(p.read_text()) for p in sorted(HERE.glob("*/run.json"))}


def medians(rec: dict) -> dict[tuple[str, int], float]:
    return {(s["domain"], s["threads"]): s["median"] for s in rec["stats"]}


def main() -> None:
    runs = load()
    print("## Runs\n")
    print("| run | commit | dirty | power | pmset % | so sha256 | started |")
    print("|---|---|---|---|---|---|---|")
    for name, r in runs.items():
        print(
            f"| {name} | {r['tree']['commit'][:7]} | {r['tree']['dirty']} | {r['power']['state']} "
            f"| {r['power']['percent']} | {(r['build']['so_sha256'] or '')[:12]} | {r['started'][11:19]} |"
        )
    print("\n## Quality\n")
    print("| run | domain | worst angle | max degree | within tol | Delaunay checked / ambiguous / violations | mesh sha256 |")
    print("|---|---|---:|---:|---|---|---|")
    for name, r in runs.items():
        for d, q in sorted(r["quality"].items()):
            print(
                f"| {name} | {d} | {q['worst_angle']:.4f} | {q['max_degree']} | {q['within_tolerance']} "
                f"| {q['delaunay_checked']} / {q['delaunay_ambiguous']} / {q['delaunay_violations']} "
                f"| {q['mesh_sha256'][:16]} |"
            )
    base = [medians(r) for n, r in runs.items() if n.startswith("base-")]
    new = [medians(r) for n, r in runs.items() if n.startswith("15c-1-")]
    keys = sorted(base[0])
    print(f"\n## Refine time, median of the per-run medians ({len(base)} base runs, {len(new)} branch runs)\n")
    print("threads 0 is the CLI default (hardware_concurrency).\n")
    print("| domain | threads | base a130f7c (s) | base min-max | 15c-1 (s) | 15c-1 min-max | change |")
    print("|---|---:|---:|---|---:|---|---:|")
    changes: list[float] = []
    for k in keys:
        b = [m[k] for m in base]
        n = [m[k] for m in new]
        bm, nm = statistics.median(b), statistics.median(n)
        ch = (nm / bm - 1) * 100
        changes.append(ch)
        print(
            f"| {k[0]} | {k[1]} | {bm:.4f} | {min(b):.4f}-{max(b):.4f} | {nm:.4f} "
            f"| {min(n):.4f}-{max(n):.4f} | {ch:+.1f} % |"
        )
    print(
        f"\nOver {len(changes)} cells: median change {statistics.median(changes):+.2f} %, "
        f"range {min(changes):+.1f} % to {max(changes):+.1f} %, "
        f"{sum(c > 5 for c in changes)} cells above +5 %, {sum(c < -5 for c in changes)} below -5 %."
    )
    # Per pair and per domain: the median over thread counts of the change
    # branch/base, so an order effect (pairs 1-3 base first, 4-5 branch first)
    # shows. Pair n is base-a130f7c-r<n> with 15c-1-r<n>.
    print("\n## Per pair: median over the 21 thread cells of (branch / base - 1)\n")
    print("| pair | order | quarter | tile | verdict stored in the branch run.json |")
    print("|---|---|---:|---:|---|")
    for n in sorted({k.rsplit("-r", 1)[1] for k in runs}):
        b, w = medians(runs[f"base-a130f7c-r{n}"]), medians(runs[f"15c-1-r{n}"])
        order = "base first" if runs[f"base-a130f7c-r{n}"]["started"] < runs[f"15c-1-r{n}"]["started"] else "branch first"
        per = {
            d: statistics.median((w[k] / b[k] - 1) * 100 for k in keys if k[0] == d)
            for d in ("quarter", "tile")
        }
        verdict = "; ".join(runs[f"15c-1-r{n}"]["verdict"])
        print(f"| {n} | {order} | {per['quarter']:+.1f} % | {per['tile']:+.1f} % | {verdict} |")
    print(
        "\nPairs 4-5 ran branch first, so each branch run stored NO BASELINE; "
        "`bench.py compare` against its pair's base is in logs/pairs_reversed.log."
    )
    # Ceiling: speed-up over 1 thread at 20, and the best, on the pooled medians.
    print("\n## Ceiling (pooled medians; 2026-09-26 reference 2.2x)\n")
    print("| domain | build | 1 -> 20 threads | best | at threads |")
    print("|---|---|---:|---:|---:|")
    for d in ("quarter", "tile"):
        for name, group in (("base a130f7c", base), ("15c-1", new)):
            med = {k[1]: statistics.median(m[k] for m in group) for k in keys if k[0] == d and k[1] > 0}
            best = max(med, key=lambda t: med[1] / med[t])
            print(f"| {d} | {name} | {med[1] / med[20]:.2f}x | {med[1] / med[best]:.2f}x | {best} |")
    # Noise floor: the same build against itself, every pair of base runs.
    self_ch = [
        (b2[k] / b1[k] - 1) * 100 for i, b1 in enumerate(base) for b2 in base[i + 1 :] for k in keys
    ] + [(n2[k] / n1[k] - 1) * 100 for i, n1 in enumerate(new) for n2 in new[i + 1 :] for k in keys]
    print(
        f"Same build against itself (base-base and branch-branch run pairs, {len(self_ch)} cells): "
        f"range {min(self_ch):+.1f} % to {max(self_ch):+.1f} %, "
        f"{sum(abs(c) > 5 for c in self_ch)} cells beyond +-5 %."
    )


if __name__ == "__main__":
    main()
