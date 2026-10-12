"""Quick check: a head record against a base record, both saved by
bench_quick.py --save-baseline and copied to quick/, per case: total, mesh
SHA-256 equality, max_error, and the features clip, decode, refine and write
phases. Adapted from 2026-10-11/cleanup-speed/scripts/quick_delta.py.
Usage: python quick_delta.py [base.json head.json]; default: the back-to-back
pair base-81a551b3.json and head-e2470ae7.json. Run from anywhere."""
import json
import sys
from pathlib import Path

here = Path(__file__).resolve().parent.parent
names = sys.argv[1:3] or ["base-81a551b3.json", "head-e2470ae7.json"]
b = json.loads((here / "quick" / names[0]).read_text())
h = json.loads((here / "quick" / names[1]).read_text())
PHASES = ["features clip", "decode", "refine", "write: encode", "write: disk"]


def cell(c: dict, hc: dict, k: str) -> str:
    if k not in c["measures"]:
        return "-"
    pb, ph = c["measures"][k]["median"], hc["measures"][k]["median"]
    return f"{pb:.3f} / {ph:.3f} ({100 * (ph / pb - 1):+.1f} %)"


print(f"base {b['commit'][:8]} (numpy {b['machine']['numpy']}, {b['power']}) -> "
      f"head {h['commit'][:8]} (numpy {h['machine']['numpy']}, {h['power']})\n")  # fmt: skip
print("| case | base total s | head total s | change | mesh same | max_error base / head | "
      + " | ".join(PHASES) + " |")  # fmt: skip
print("|---|---:|---:|---:|---|---|" + "---|" * len(PHASES))
for n, c in b["cases"].items():
    hc = h["cases"][n]
    tb, th = c["measures"]["total"]["median"], hc["measures"]["total"]["median"]
    print(f"| {n} | {tb:.3f} | {th:.3f} | {100 * (th / tb - 1):+.1f} % | "
          f"{c['mesh_sha256'] == hc['mesh_sha256']} | {c['max_error']:.3f} / {hc['max_error']:.3f} | "
          + " | ".join(cell(c, hc, k) for k in PHASES) + " |")  # fmt: skip

print("\nSpreads (max - min) / median of the total, base / head:")
for n, c in b["cases"].items():
    sb, sh = c["measures"]["total"], h["cases"][n]["measures"]["total"]
    print(f"  {n}: {100 * (sb['max'] - sb['min']) / sb['median']:.1f} % / "
          f"{100 * (sh['max'] - sh['min']) / sh['median']:.1f} %")  # fmt: skip
