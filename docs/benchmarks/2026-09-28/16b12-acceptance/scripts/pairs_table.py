"""The four bench.py pairs: per-cell medians, pair changes, the pooled change,
same-build drift, the ceiling and quality. Run from this directory."""
import json, statistics as S
runs = ["base-14f5fe3", "16b12", "base-14f5fe3-r2", "16b12-r2", "base-14f5fe3-r3", "16b12-r3", "base-14f5fe3-r4", "16b12-r4"]
R = {r: json.load(open(f"{r}/run.json")) for r in runs}
for r in runs:
    d = R[r]
    print(f"{r}: {d['tree']['commit'][:7]} so {d['build']['so_sha256'][:12]} {d['power']['state']} {d['power']['raw'].split(chr(9))[1].split(';')[0]}->{d['power']['raw'].split(chr(9))[-1].split(';')[0]} verdict {'; '.join(d['verdict'])[:120]}")
    for k, v in d["quality"].items():
        print(f"   {k}: worst {v['worst_angle']:.4f} maxdeg {v['max_degree']} tol {v['within_tolerance']} delaunay {v['delaunay_violations']}/{v['delaunay_checked']} sha {v['mesh_sha256'][:16]}")
med = {r: {(s["domain"], s["threads"]): s["median"] for s in R[r]["stats"]} for r in runs}
keys = sorted(med[runs[0]]); B, N = runs[0::2], runs[1::2]
pooled = []
print("| domain | threads | base r1..r4 (ms) | 16b-1/2 r1..r4 (ms) | pair changes | pooled |")
print("|---|---:|---|---|---|---:|")
for k in keys:
    b = [med[r][k] for r in B]; n = [med[r][k] for r in N]
    ch = [(y / x - 1) * 100 for x, y in zip(b, n)]
    p = (S.median(n) / S.median(b) - 1) * 100; pooled.append(p)
    if k[1] in (0, 1, 2, 4, 8, 10, 16, 20):
        print(f"| {k[0]} | {k[1] or 'default (10)'} | {' / '.join(f'{x*1000:.1f}' for x in b)} | {' / '.join(f'{x*1000:.1f}' for x in n)} | {', '.join(f'{c:+.1f}' for c in ch)} | {p:+.1f} % |")
print(f"pooled over {len(pooled)} cells: median {S.median(pooled):+.2f} %, min {min(pooled):+.1f}, max {max(pooled):+.1f}, cells > +5 %: {sum(p > 5 for p in pooled)}")
for grp in (B, N):
    worst = max(((med[b][k] / med[a][k] - 1) * 100, a, b, k) for a in grp for b in grp if a != b for k in keys)
    print(f"same-build drift within {grp[0]}*: max {worst[0]:+.1f} % ({worst[1]} -> {worst[2]}, {worst[3]})")
for r in runs:
    out = []
    for dom in ("tile", "quarter"):
        one = med[r][(dom, 1)]; best = min((med[r][(dom, t)], t) for t in range(1, 21))
        out.append(f"{dom} {one/best[0]:.2f}x @{best[1]} ({one/med[r][(dom, 20)]:.2f}x @20)")
    print(f"ceiling {r}: {'; '.join(out)}")
