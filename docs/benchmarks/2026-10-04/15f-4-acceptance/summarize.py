"""Tables for 15f-4-acceptance.md (@perf). python summarize.py SCRATCH > tables.md
Sides: base c193cb1 (B), 15f-3 107cdb4 (P), 15f-4 5eeb87a (N); batch 1 B P N, batch 2 N P B.
bench: pooled medians over a side's two runs (10 samples a cell). Bygdin 1 m: the six
timed --binary runs a side (time -l real, peak RSS) and their --stats phases. Velhas 5 m:
velhas.py's four timed runs a side. SCRATCH holds the Bygdin ASCII meshes for the hash."""
import hashlib
import json
import re
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = HERE.parent
S = Path(sys.argv[1])
SIDES = {"base": "B c193cb1", "15f3": "P 15f-3", "15f4": "N 15f-4"}
pc = lambda a, b: f"{(a / b - 1) * 100:+.1f} %"  # noqa: E731

rec = {s: [json.loads((D / f"15f-4-{s}-b{b}" / "run.json").read_text()) for b in (1, 2)] for s in SIDES}
THREADS = sorted({x["threads"] for x in rec["base"][0]["samples"]})


def pool(s, d, t, key):
    return st.median([x[key] for r in rec[s] for x in r["samples"] if x["domain"] == d and x["threads"] == t])


print("### 1 m benchmark (bench.py), pooled medians, s\n")
for key in ("refine_s", "app_s", "proc_s"):
    print(f"#### {key}\n\n| domain | threads | B | P | N | N/P | N/B | P/B |\n|---|---:|---:|---:|---:|---:|---:|---:|")
    for d in ("tile", "quarter"):
        for t in THREADS:
            b, p, n = (pool(s, d, t, key) for s in SIDES)
            print(f"| {d} | {'default' if t == 0 else t} | {b:.4f} | {p:.4f} | {n:.4f} | {pc(n, p)} | {pc(n, b)} | {pc(p, b)} |")
    print()
print("| run | commit | `_core` | power | tile mesh | quarter mesh | quarter worst angle | Delaunay viol. tile/quarter |")
print("|---|---|---|---|---|---|---:|---|")
for s in SIDES:
    for r in rec[s]:
        q = r["quality"]
        print(f"| {r['label']} | {r['tree']['commit'][:7]} | {r['build']['so_sha256'][:12]} | {r['power']['state']} {r['power']['percent']} % | "
              f"`{q['tile']['mesh_sha256'][:16]}` | `{q['quarter']['mesh_sha256'][:16]}` | {q['quarter']['worst_angle']:.4f} | "
              f"{q['tile']['delaunay_violations']}/{q['quarter']['delaunay_violations']} |")
print()

L = HERE / "bygdin"


def phases(p):
    t = p.read_text()
    out = {m[1]: float(m[2]) for m in re.finditer(r"^\| ([^|*]+?) \| ([\d.]+) \| [\d.]+ % \|$", t, re.M)}
    tot = re.search(r"^\| \*\*total\*\* \| \*\*([\d.]+)\*\*", t, re.M)
    out["total"] = float(tot[1]) if tot else float("nan")
    return out


print("### Bygdin 1 m with CORINE (6 timed runs a side)\n")
print("| item | B | P | N | N/P | N/B |\n|---|---:|---:|---:|---:|---:|")
col = {}
for s in SIDES:
    runs = sorted(L.glob(f"{s}_b*_run*.err"))
    real = [float(re.search(r"([\d.]+) real", r.read_text())[1]) for r in runs]
    rss = [int(re.search(r"(\d+)\s+maximum resident", r.read_text())[1]) / 2**30 for r in runs]
    ph = [phases(Path(str(r)[:-4] + ".stats.md")) for r in runs]
    col[s] = {"process wall (s)": st.median(real), "peak RSS (GiB), max": max(rss),
              **{f"phase {k} (s)": st.median(x.get(k, 0.0) for x in ph)
                 for k in ("refine", "edge strip: generate", "edge strip: scan (parallel)",
                           "edge strip: split + flip (serial)", "other", "total")}, "n": len(runs),
              "range": f"{min(real):.2f}-{max(real):.2f}"}
for k in [k for k in col["base"] if k not in ("n", "range")]:
    b, p, n = (col[s][k] for s in SIDES)
    print(f"| {k} | {b:.3f} | {p:.3f} | {n:.3f} | {pc(n, p) if p else '-'} | {pc(n, b) if b else '-'} |")
print(f"| process wall range | {col['base']['range']} | {col['15f3']['range']} | {col['15f4']['range']} | | |")
for s in SIDES:
    data = (S / "bygdin" / f"{s}_ascii.vtk").read_bytes()
    print(f"\n- Bygdin 1 m ASCII mesh sha256 (from POINTS) {s}: `{hashlib.sha256(data[data.index(b'\nPOINTS ') + 1:]).hexdigest()}`", end="")
print("\n")

print("### Velhas 5 m, EPSG:31983 (4 timed runs a side)\n")
V = HERE / "velhas"
print("| item | B | P | N | N/P | N/B |\n|---|---:|---:|---:|---:|---:|")
vc = {}
for s in SIDES:
    timed, ph, sha = [], [], set()
    for b in (1, 2):
        r = json.loads((V / f"{s}-b{b}" / "results.json").read_text())["runs"][0]
        timed += [t["real_s"] for t in r["timed"]]
        ph += [phases(V / f"{s}-b{b}" / "logs" / f"t5_run{j}.stats.md") for j in (1, 2)]
        sha.add(r["mesh_sha256"])
    vc[s] = {"process wall (s)": st.median(timed),
             **{f"phase {k} (s)": st.median(x.get(k, 0.0) for x in ph)
                for k in ("refine", "final check: scan (parallel)", "final check: split + flip (serial)", "other", "total")},
             "sha": sha}
for k in [k for k in vc["base"] if k != "sha"]:
    b, p, n = (vc[s][k] for s in SIDES)
    print(f"| {k} | {b:.3f} | {p:.3f} | {n:.3f} | {pc(n, p) if p else '-'} | {pc(n, b) if b else '-'} |")
for s in SIDES:
    print(f"\n- Velhas 5 m mesh sha256 {s}: {', '.join('`' + x + '`' for x in sorted(vc[s]['sha']))}", end="")
print()
