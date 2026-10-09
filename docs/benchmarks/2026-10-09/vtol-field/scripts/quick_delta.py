"""Quick check: the head record (docs/benchmarks/quick/baseline-ac.json) against the base's."""
import json
from pathlib import Path

here = Path(__file__).resolve().parent.parent
b = json.loads((here / "quick/base-4cd7e050.json").read_text())
h = json.loads((here.parent.parent / "quick/baseline-ac.json").read_text())
print("| case | base total s | head total s | change | refine base s | refine head s | mesh same |")
print("|---|---:|---:|---:|---:|---:|---|")
for n, c in b["cases"].items():
    hc = h["cases"][n]
    tb, th = c["measures"]["total"]["median"], hc["measures"]["total"]["median"]
    rb, rh = c["measures"]["refine"]["median"], hc["measures"]["refine"]["median"]
    print(f"| {n} | {tb:.3f} | {th:.3f} | {100 * (th / tb - 1):+.1f} % | {rb:.3f} | {rh:.3f} | "
          f"{c['mesh_sha256'] == hc['mesh_sha256']} |")
