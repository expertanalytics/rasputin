"""Tables for README.md (@perf). python summarize.py [MESHES] > tables.md
Part 1 from the bench runs in 2026-10-19/ and the layout sweep; part 2 from
logs/ (velhas.sh). MESHES: velhas.sh's scratch directory, for the uncut runs'
point counts (left out when absent)."""
import json
import re
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MESHES = Path(sys.argv[1]) if len(sys.argv) > 1 else None
MB = 1e6
pc = lambda n, b: f"{(n / b - 1) * 100:+.1f} %"

# --- part 1: bench.py ---
RUNS = HERE / "2026-10-19"
SIDES = {"M": ["master-r1", "master-r2"], "B": [f"base-r{i}" for i in range(1, 5)],
         "N": [f"head-r{i}" for i in range(1, 5)]}  # fmt: skip
rec = {k: json.loads((RUNS / k / "run.json").read_text()) for v in SIDES.values() for k in v}


def pool(side, d, t, key):
    return st.median([s[key] for k in SIDES[side] for s in rec[k]["samples"]
                      if s["domain"] == d and s["threads"] == t])  # fmt: skip


print("## Bench (bench.py), pooled medians: M 2 runs, B 4, N 4, 5 samples a cell each\n")
print("| domain | threads | measure | M | B | N | N/B | N/M |")
print("|---|---:|---|---:|---:|---:|---:|---:|")
for d in ("tile", "quarter"):
    for t in (0, 1, 8, 20):
        for key in ("refine_s", "proc_s"):
            m, b, n = (pool(s, d, t, key) for s in "MBN")
            print(f"| {d} | {'default' if t == 0 else t} | {key} | {m:.4f} | {b:.4f} | {n:.4f} | {pc(n, b)} | {pc(n, m)} |")
print("\n| run | commit | `_core` sha256 | power | tile mesh | quarter mesh |")
print("|---|---|---|---|---|---|")
for k, r in rec.items():
    q = r["quality"]
    print(f"| {k} | {r['tree']['commit'][:8]} | {r['build']['so_sha256'][:8]} | {r['power']['state']} | "
          f"{q['tile']['mesh_sha256'][:8]} | {q['quarter']['mesh_sha256'][:8]} |")  # fmt: skip
for f in sorted(HERE.glob("layout_sweep_t*.tsv")):
    rows = [ln.split("\t") for ln in f.read_text().splitlines()]
    print(f"\n### {f.name}: 32 layouts (PAD 0 to 248 bytes), one sample each\n")
    print("| domain | tree | refine median | refine mean | proc median | proc mean |")
    print("|---|---|---:|---:|---:|---:|")
    for d in ("tile", "quarter"):
        for name in ("master", "base", "head"):
            r = [float(x[4]) for x in rows if x[0] == name and x[1] == d]
            p = [float(x[6]) for x in rows if x[0] == name and x[1] == d]
            print(f"| {d} | {name} | {st.median(r):.4f} | {st.mean(r):.4f} | {st.median(p):.4f} | {st.mean(p):.4f} |")

# --- part 2: Velhas ---
LOGS = HERE / "logs"


def timel(err):
    t = err.read_text()
    g = lambda pat: float(re.search(pat, t).group(1))
    return g(r"([\d.]+) real"), g(r"(\d+)\s+maximum resident set size"), g(r"(\d+)\s+peak memory footprint")


def points(name):
    if MESHES is None:
        return None
    f = MESHES / name / "velhas.vtk"
    if not f.exists():
        return None
    with f.open("rb") as fh:
        head = fh.read(4096)
    m = re.search(rb"POINTS (\d+)", head)
    return int(m.group(1)) if m else None


print("\n## Velhas: every run\n")
print("| run | wall s | max RSS MB | peak footprint MB | pieces | points | start inserts (sum; median; range) "
      "| refine s (sum) | decode s | largest piece peak MB | setup peak MB | largest cell dx·dy·b(T) MB |")
print("|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|")
for err in sorted(LOGS.glob("*.err"), key=lambda p: (p.stem.split("_")[0], -float(p.stem.split("_")[1][1:]),
                                                       int(p.stem.split("_")[2][1:]))):  # fmt: skip
    n = err.stem
    wall, rss, foot = timel(err)
    lines = [json.loads(x) for x in (LOGS / f"{n}.jsonl").read_text().splitlines()]
    calls, tail = lines[:-1], lines[-1]
    ins = [c["quality_inserted"] for c in calls]
    peaks = [max(c["pre_peak"], c["in_peak"]) for c in calls[1:]]  # the first "pre" holds the setup
    first = calls[0]
    setup = first["pre_peak"]
    p1 = first["in_peak"]
    piece_peak = max([p1, *peaks]) if len(calls) > 1 else None
    decode = sum(ph.get("decode", 0.0) for ph in tail["phases"])
    idx = LOGS / f"{n}.index.json"
    cell, count, pts = "", len(calls), points(n)
    if idx.exists():
        j = json.loads(idx.read_text())
        p = j["partition"]
        cell = f"{p['dx'] * p['dy'] * p['bytes_per_node'] / MB:.0f}"
        count = len(j["pieces"])
        pts = sum(e["counts"]["vertices"] for e in j["pieces"])  # shared seam vertices counted twice
    print(f"| {n} | {wall:.2f} | {rss / MB:.0f} | {foot / MB:.0f} | {count} | {pts if pts is not None else ''} | "
          f"{sum(ins)}; {st.median(ins):.0f}; {min(ins)}-{max(ins)} | {sum(c['seconds'] for c in calls):.2f} | "
          f"{decode:.2f} | {piece_peak / MB if piece_peak else 0:.0f} | {setup / MB:.0f} | {cell} |")  # fmt: skip
