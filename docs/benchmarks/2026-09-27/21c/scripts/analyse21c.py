"""analyse21c.py DATA_DIR -- the README's tables, from data/sims/*.json."""
import json
import statistics as st
import sys
from pathlib import Path

d = Path(sys.argv[1])
runs = {p.stem: json.loads(p.read_text()) for p in sorted(d.glob("*.json"))}
order = ["today_foot", "C", "Cthin", "Cmis", "Cthin2", "A1", "A1_seed1", "A1_seed2", "A1_seed3",
         "Hser", "Hser_seed1", "Hser_seed2", "Hser_seed3"]


def rs(r):
    return r["sim"]["rounds"]


print("## Summary\n")
print("Triangles: the refine output's count (`RefineOutcome.triangles`); the written tile VTK has "
      "8,078-8,084 fewer in every mode (column `vtk_triangles` in the data).\n")
print("| domain | mode | triangles | vs today | rounds | inserted | flips | worst angle | "
      "share < 1 deg | max degree | max error | full rescan: max error / over tol | "
      "Delaunay violations / edges | consistency |")
print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|")
for dom in ["quarter", "tile"]:
    base = runs[f"{dom}_today_foot"]["sim"]["triangles"]
    for m in order:
        r = runs[f"{dom}_{m}"]
        s, q, c = r["sim"], r["quality"], r["child"]
        viol = sum(x["cavity_violations"] for x in rs(r))
        print(f"| {dom} | {m.replace('_foot', '')} | {s['triangles']:,} | "
              f"{100 * (s['triangles'] / base - 1):+.2f} % | {c['rounds']} | {c['inserted']:,} | "
              f"{c['flips']:,} | {q['worst_angle']:.4f} | {q['share_under_1']:.2e} | "
              f"{q['max_degree']} | {c['max_error']:.4f} | {s['full_rescan_max_error']:.4f} / "
              f"{s['full_rescan_needs_split']} | {q['delaunay_violations']} / "
              f"{q['delaunay_checked']:,} | {viol} |")

print("\n## Batch semantics per mode (sums over rounds)\n")
print("| domain | mode | marked | inserted | skipped (slot written) | skipped (edge) | "
      "deferred (split-set conflict) | deferred (thinning) |")
print("|---|---|---:|---:|---:|---:|---:|---:|")
for dom in ["quarter", "tile"]:
    for m in order:
        R = rs(runs[f"{dom}_{m}"])
        f = lambda k: sum(x[k] for x in R)  # noqa: E731
        print(f"| {dom} | {m.replace('_foot', '')} | {f('marked'):,} | {f('inserted'):,} | "
              f"{f('skipped_touched'):,} | {f('skipped_edge'):,} | {f('deferred_conflict'):,} | "
              f"{f('deferred_thin'):,} |")

print("\n## Option C: flip rounds\n")
print("| domain | mode | refine rounds | flip rounds, total | flip rounds per refine round: "
      "median / max | flips | first flip round: candidate edges / flips, largest round |")
print("|---|---|---:|---:|---|---:|---|")
for dom in ["quarter", "tile"]:
    for m in ["C", "Cthin", "Cmis", "Cthin2"]:
        R = [x for x in rs(runs[f"{dom}_{m}"]) if x["inserted"]]
        big = max(R, key=lambda x: x["inserted"])
        fr = [x["flip_rounds"] for x in R]
        print(f"| {dom} | {m} | {len(rs(runs[f'{dom}_{m}']))} | {sum(fr)} | {st.median(fr)} / "
              f"{max(fr)} | {sum(x['flips'] for x in R):,} | round {big['round']} "
              f"({big['inserted']:,} inserted): {big['flip_candidates_first']:,} / "
              f"{big['flips_first']:,} |")

print("\n## Option A1: sub-rounds\n")
print("| domain | seed | rounds | sub-rounds, total | per round: median / max | "
      "big rounds (>= 5,000 marks): first sub-round winners as a share of marks | "
      "big rounds: sub-rounds until 90 % of that round's inserts | footprint (cavity + ring) "
      "mean, p99, max |")
print("|---|---|---:|---:|---|---|---|---|")
for dom in ["quarter", "tile"]:
    for m in ["A1", "A1_seed1", "A1_seed2", "A1_seed3"]:
        r = runs[f"{dom}_{m}"]
        R = [x for x in rs(r) if x["subrounds"]]
        sub = [x["subrounds"] for x in R]
        big = [x for x in R if x["marked"] >= 5000]
        w1 = [x["winners_first"] / x["marked"] for x in big]
        to90 = []
        for x in big:
            acc, tot = 0, x["inserted"]
            for k, w in enumerate(x["subround_winners"], 1):
                acc += w
                if acc >= 0.9 * tot:
                    to90.append(k)
                    break
        s = r["sim"]
        fp_mean = sum(x["fp_sum"] for x in R) / max(1, sum(x["inserted"] for x in R))
        print(f"| {dom} | {m.split('seed')[-1] if 'seed' in m else 0} | {len(rs(r))} | {sum(sub)} | "
              f"{st.median(sub)} / {max(sub)} | {100 * min(w1):.0f}-{100 * max(w1):.0f} % of marks "
              f"({len(big)} rounds) | {min(to90)}-{max(to90)} | {fp_mean:.1f}, {s['fp_p99']}, "
              f"{s['fp_max']} |")

print("\n## Today's order: dependence depth and footprint extent, per round\n")
bs = [16, 32, 64, 128]
for dom in ["quarter", "tile"]:
    r = runs[f"{dom}_today_foot"]
    s = r["sim"]
    print(f"\n### {dom}\n")
    print(f"Over all insertions: footprint slots p50 {s['fp_p50']}, p90 {s['fp_p90']}, "
          f"p99 {s['fp_p99']}, max {s['fp_max']}; bounding-box larger side (nodes) p50 "
          f"{s['side_p50']}, p90 {s['side_p90']}, p99 {s['side_p99']}, max {s['side_max']}.\n")
    R = [x for x in rs(r) if x["inserted"]]
    tot = sum(x["inserted"] for x in R)
    print("| B | leaves the B/2 halo | crosses a block boundary (no halo) |")
    print("|---:|---:|---:|")
    for b in [16, 32, 64, 128, 256, 512, 1024]:
        lh = sum(x.get(f"leaves_halo_{b}", 0) for x in R)
        cr = sum(x.get(f"crosses_{b}", 0) for x in R)
        print(f"| {b} | {'%d (%.3f %%)' % (lh, 100 * lh / tot) if b <= 128 else '-'} | "
              f"{cr:,} ({100 * cr / tot:.1f} %) |")
    print()
    print("| round | inserted | depth | depth / inserted | fp mean | side p50 / p99 / max | "
          + " | ".join(f"halo B={b}" for b in bs) + " | ins./block B=32 max / mean | "
          "B=64 max / mean |")
    print("|" + "---:|" * (8 + len(bs)))
    for x in R:
        n = x["inserted"]
        print(f"| {x['round']} | {n:,} | {x['depth']:,} | {x['depth'] / n:.2f} | "
              f"{x['fp_sum'] / n:.1f} | {x['side_p50']} / {x['side_p99']} / {x['side_max']} | "
              + " | ".join(f"{100 * x.get(f'leaves_halo_{b}', 0) / n:.2f} %" for b in bs)
              + f" | {x['block_max_32']} / {x['block_mean_32']:.1f} | "
              f"{x['block_max_64']} / {x['block_mean_64']:.1f} |")
