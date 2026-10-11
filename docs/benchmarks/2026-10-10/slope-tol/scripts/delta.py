"""Per (domain, threads): median refine_s of bench-base and bench-head, and the change."""
import json
import sys
from pathlib import Path

here = Path(__file__).resolve().parent.parent
b, h = (json.loads((here / d / "run.json").read_text()) for d in ("bench-base", "bench-head"))


def med(r):
    return {(s["domain"], s["threads"]): s["median"] for s in r["stats"]}


mb, mh = med(b), med(h)
print("| domain | threads | base median s | head median s | change |")
print("|---|---:|---:|---:|---:|")
ch = []
for k in sorted(mb):
    d = 100 * (mh[k] / mb[k] - 1)
    ch.append(d)
    print(f"| {k[0]} | {k[1]} | {mb[k]:.4f} | {mh[k]:.4f} | {d:+.1f} % |")
print(f"\ncells {len(ch)}; change min {min(ch):+.1f} %, max {max(ch):+.1f} %, "
      f"median {sorted(ch)[len(ch)//2]:+.1f} %", file=sys.stderr)
for d in ("tile", "quarter"):
    print(d, b["quality"][d]["mesh_sha256"] == h["quality"][d]["mesh_sha256"], file=sys.stderr)
