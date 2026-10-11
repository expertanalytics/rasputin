"""Quick check: the head record against the base record, both saved by
bench_quick.py --save-baseline and copied to quick/, per case: total, the features
phases, decode and refine; then the section 7 ratios and each case's phases at 10 %
or more of its run (the hotspot line is 40 %).
Usage: python quick_delta.py [base.json head.json]; default: the back-to-back pair
base-3e6f87eb.json and head-357808b1.json. Run from anywhere."""
import json
import sys
from pathlib import Path

here = Path(__file__).resolve().parent.parent
names = sys.argv[1:3] or ["base-3e6f87eb.json", "head-357808b1.json"]
b = json.loads((here / "quick" / names[0]).read_text())
h = json.loads((here / "quick" / names[1]).read_text())
PHASES = ["features read", "features clip", "decode", "refine"]


def cell(c: dict, hc: dict, k: str) -> str:
    if k not in c["measures"]:
        return "-"
    pb, ph = c["measures"][k]["median"], hc["measures"][k]["median"]
    return f"{pb:.3f} / {ph:.3f} ({100 * (ph / pb - 1):+.1f} %)"


print("| case | base total s | head total s | change | mesh same | " + " | ".join(PHASES) + " |")
print("|---|---:|---:|---:|---|" + "---|" * len(PHASES))
for n, c in b["cases"].items():
    hc = h["cases"][n]
    tb, th = c["measures"]["total"]["median"], hc["measures"]["total"]["median"]
    print(f"| {n} | {tb:.3f} | {th:.3f} | {100 * (th / tb - 1):+.1f} % | "
          f"{c['mesh_sha256'] == hc['mesh_sha256']} | " + " | ".join(cell(c, hc, k) for k in PHASES) + " |")  # fmt: skip

print("\nSection 7 targets (ratio head median / base median):")
for n, k, limit in [("lagan", "features clip", 0.65), ("lagan", "total", 0.85),
                    ("numedalslagen", "features clip", 1.0)]:  # fmt: skip
    pb = b["cases"][n]["measures"][k]["median"]
    ph = h["cases"][n]["measures"][k]["median"]
    print(f"  {n} {k}: {pb:.3f} -> {ph:.3f} s, ratio {ph / pb:.3f}, limit {limit}: "
          f"{'met' if ph / pb <= limit else 'NOT MET'}")

for n in ("lagan", "numedalslagen"):
    for tag, rec in (("base", b), ("head", h)):
        m = rec["cases"][n]["measures"]
        tot = m["total"]["median"]
        print(f"\n{n} {tag} ({rec['commit'][:8]}): total {tot:.3f} s; phases at 10 % or more:")
        for k, s in sorted(m.items(), key=lambda kv: -kv[1]["median"]):
            if k != "total" and s["median"] / tot >= 0.10:
                print(f"  {k}: {s['median']:.3f} s ({100 * s['median'] / tot:.1f} %)")
