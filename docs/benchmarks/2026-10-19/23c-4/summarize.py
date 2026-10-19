"""Refine medians (s) from each bench run's README: t=0 (default), 1 and 12, and the ceiling t1/t0."""
import re
runs = {"23c-3 head-r4 (accepted)": "23c-3/2026-10-19/head-r4", "23c-4 head-r1": "23c-4/2026-10-19/head-r1",
        "23c-3 base-r1": "23c-4/2026-10-19/base-r1", "23c-4 head-r2": "23c-4/2026-10-19/head-r2"}
print("| run | tile t0 | tile t1 | tile t12 | tile ceiling | quarter t0 | quarter t1 | quarter t12 | quarter ceiling |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
for k, d in runs.items():
    t = open(d + "/README.md").read()
    cells = []
    for dom in ("tile", "quarter"):
        sec = t.split(f"## `{dom}`")[1].split("\n## ")[0]
        rows = {int(a): float(b) for a, b in re.findall(r"^\| (\d+) \| \d+ \| ([\d.]+) \|", sec, re.M)}
        cells += [f"{rows[0]:.4f}", f"{rows[1]:.4f}", f"{rows[12]:.4f}", f"{rows[1] / rows[0]:.2f}x"]
    print(f"| {k} | " + " | ".join(cells) + " |")
