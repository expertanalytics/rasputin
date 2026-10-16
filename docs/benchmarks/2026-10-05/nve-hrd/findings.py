"""Step 4's per-station tables (@perf), from results.csv, checks.csv and,
when present, window_check_12km.csv and bypass.csv. Markdown on stdout.

Usage: findings.py HERE"""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def load(path: Path) -> dict[str, dict[str, str]]:
    return {r["station"]: r for r in csv.DictReader(open(path, encoding="utf-8"))} if path.exists() else {}


def pct(v: str) -> str:
    return "" if v == "" else f"{100 * float(v):.1f}"


def main(here: Path) -> None:
    rows = load(here / "results.csv")
    chk = load(here / "checks.csv")
    byp = load(here / "bypass.csv")
    print("### Every `uncertain` station, with its causes\n")
    print("Causes from the sensitivity (`causes`): `swing` the area changes over 5 % within "
          "the position uncertainty `U` up or down the river; `downstream_unread` the river "
          "was not read to `U` below the gauge; `chain_not_draining` a node of the burnt path "
          "does not drain into the next, or the counts do not rise; `chain_end_open` the "
          "path's lowered end found no way down within the cap. Largest step: the largest rise "
          "in area between two path nodes within `U`, and where (m, + downstream). "
          "Bypass: in the burnt DEM, the largest count within 30 m of the placed node, against "
          "the placed node's own (`bypass.py`).\n")
    print("| station | name | NVE km² | ours km² | NVE's in ours % | causes | largest step km² at m | U m | lake line / lake above P | bypass km² |")
    print("|---|---|---:|---:|---:|---|---|---:|---|---:|")
    for s, r in rows.items():
        if r["class"] != "uncertain":
            continue
        lake = "lake line" if r["lake"] == "True" else ""
        if chk.get(s, {}).get("lake_above_p") == "True":
            lake = (lake + ", " if lake else "") + "lake above P"
        step = f"{float(r['largest_step']):.2f} at {float(r['largest_step_at_m']):.0f}"
        b = byp.get(s)
        bp = "" if b is None else f"{float(b['near_max_km2']):.2f}"
        print(f"| {s} | {r['name']} | {float(r['reference_area_km2']):.1f} | {float(r['fine_area_km2']):.3f} | "
              f"{pct(r['nve_in_ours'])} | {r['causes'].replace(';', ', ')} | {step} | "
              f"{float(r['uncertainty_m']):.0f} | {lake} | {bp} |")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
