"""Quick check: the head record (docs/benchmarks/quick/baseline-ac.json, saved at head)
against the base's (quick/base-67af3081.json, the baseline master carried), per case:
total, refine, and the two point-scan paths (edge strip and final check / check points).
Run from the worktree root."""
import json
from pathlib import Path

here = Path(__file__).resolve().parent.parent
b = json.loads((here / "quick/base-67af3081.json").read_text())
h = json.loads((here.parent.parent / "quick/baseline-ac.json").read_text())
PHASES = ["refine", "edge strip: scan (parallel)", "edge strip: split + flip (serial)",
          "check points: store", "final check: scan (parallel)", "final check: split + flip (serial)"]  # fmt: skip


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
for n in h["cases"].keys() - b["cases"].keys():
    hc = h["cases"][n]
    m = hc["measures"]
    print(f"\nnew case {n}: total {m['total']['median']:.3f} s (min {m['total']['min']:.3f}, max {m['total']['max']:.3f}), "
          f"refine {m['refine']['median']:.3f} s, slope {m['slope']['median']:.4f} s, "
          f"max_error {hc['max_error']:.4f} of tolerance {hc['tolerance']}")  # fmt: skip
    for k, s in sorted(m.items(), key=lambda kv: -kv[1]["median"])[:8]:
        print(f"  {k}: {s['median']:.3f} s ({100 * s['median'] / m['total']['median']:.0f} %)")
