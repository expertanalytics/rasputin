"""Q7: refine time per output triangle, the slope's ramp against uniform 2 m on tile 6901_3.

Run from the worktree root with the worktree's venv after `bench.py run` has built
build-bench/pkg:  python <this> RASPUTIN_DATA OUT_DIR RAW_JSONL
One warm-up of each case, then 3 timed runs alternating uniform and slope, default threads.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

data, out, raw = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
pkg = Path.cwd() / "build-bench/pkg"
base = ["mesh", "--dem", f"{data}/DTM10_UTM33_20260925/6901_3_10m_z33.tif", "--binary"]
cases = {
    "uniform-2m": ["--tolerance", "2"],
    "slope": ["--tolerance", "10", "--tolerance-slope", "2", "25", "35"],
}


def run(name: str) -> dict:
    stats = out / f"{name}.md"
    argv = [sys.executable, "tools/bench.py", "_child", "--pkg", str(pkg), "--threads", "0", "--",
            *base, *cases[name], "--out", str(out / f"{name}.vtk"), "--stats", str(stats)]  # fmt: skip
    p = subprocess.run(argv, capture_output=True, text=True, check=True)
    line = next(x for x in p.stderr.splitlines() if x.startswith("BENCH "))
    rec = json.loads(line[6:])
    tri = int(re.search(r"\| output triangles \| (\d+) \|", stats.read_text()).group(1))
    pm = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout.splitlines()[0]
    return {"case": name, "triangles": tri, "power": pm, **rec}


for name in cases:
    run(name)  # warm-up
with raw.open("w") as f:
    for _ in range(3):
        for name in cases:
            f.write(json.dumps(run(name)) + "\n")
rows = [json.loads(x) for x in raw.read_text().splitlines()]
med = {}
for name in cases:
    rs = sorted(r["phases"]["refine"] for r in rows if r["case"] == name)
    tot = sorted(r["total_s"] for r in rows if r["case"] == name)
    tri = {r["triangles"] for r in rows if r["case"] == name}
    med[name] = (rs[1], tri.pop() if len(tri) == 1 else tri)
    print(name, "refine median s", round(rs[1], 4), "min", round(rs[0], 4), "max", round(rs[-1], 4),
          "total median s", round(tot[1], 3), "triangles", med[name][1],
          "us per triangle", round(1e6 * rs[1] / med[name][1], 3))  # fmt: skip
    if name == "slope":
        sl = sorted(r["phases"].get("slope", float("nan")) for r in rows if r["case"] == name)
        print("slope phase median s", round(sl[1], 4), "min", round(sl[0], 4), "max", round(sl[-1], 4))
u, s = med["uniform-2m"], med["slope"]
print("ratio slope / uniform-2m, refine per output triangle", round((s[0] / s[1]) / (u[0] / u[1]), 3))
