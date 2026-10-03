"""Per-phase peaks and per-unit bytes from `run_phases.py`'s records. A
measurement script, not production code.

    python analyse.py RUNS_DIR TAG... > tables.md

For each TAG reads ``TAG.markers.jsonl``, ``TAG.mem.csv``, ``TAG.json`` and
``TAG.counts.json`` (the counts the CLI printed: final triangles, phase 2's
insertions, phase 1's vertices) and prints, per phase, the wall time, the
footprint when it starts, its peak, and the rise. A phase's peak is the largest
of its samples, of the footprint read at every marker inside it, and of the
lifetime peak read at its last marker when that rose inside the phase (exact).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

PHASES = [  # name, start marker, end marker ((name, ev) pairs)
    ("P1 decode + assemble", ("clock: decode", "B"), ("assemble", "E")),
    ("P2 resample (blocks)", ("resample", "B"), ("resample: DemTile copy", "B")),
    ("P2 end (DemTile copy)", ("resample: DemTile copy", "B"), ("resample", "E")),
    ("start mesh", ("resample", "E"), ("refine phase 1", "B")),
    ("P3 refine phase 1", ("refine phase 1", "B"), ("refine phase 1", "E")),
    ("P4 store fill + freeze", ("final_check.run", "B"), ("refine_points", "B")),
    ("P5 phase 2 (refine_points)", ("refine_points", "B"), ("refine_points", "E")),
    ("P6 trim + write", ("refine_points", "E"), ("process", "E")),
]


def load(runs: Path, tag: str) -> tuple[list[dict], list[tuple[float, int]], dict, dict]:
    markers = [json.loads(x) for x in (runs / f"{tag}.markers.jsonl").read_text().splitlines()]
    with (runs / f"{tag}.mem.csv").open() as f:
        samples = [(float(r["t"]), int(r["footprint_bytes"])) for r in csv.DictReader(f)]
    rec = json.loads((runs / f"{tag}.json").read_text())
    counts = json.loads((runs / f"{tag}.counts.json").read_text())
    return markers, samples, rec, counts


def phases(markers: list[dict], samples: list[tuple[float, int]]) -> list[dict]:
    def find(name: str, ev: str) -> int:
        return next(i for i, m in enumerate(markers) if m["name"] == name and m["ev"] == ev)

    out = []
    for label, (bn, be), (en, ee) in PHASES:
        i, j = find(bn, be), find(en, ee)
        b, e = markers[i], markers[j]
        inner = markers[i : j + 1]
        peak = max([m["fp"] for m in inner] + [fp for t, fp in samples if b["t"] <= t <= e["t"]])
        exact = e["life_max"] > b["life_max"]
        if exact:
            peak = max(peak, e["life_max"])
        out.append({"phase": label, "seconds": e["t"] - b["t"], "start": b["fp"], "end": e["fp"],
                    "peak": peak, "exact": exact})  # fmt: skip
    return out


def main() -> None:
    runs = Path(sys.argv[1])
    for tag in sys.argv[2:]:
        markers, samples, rec, counts = load(runs, tag)
        rows = phases(markers, samples)
        print(f"\n### `{tag}`\n")
        power = rec["pmset_before"].splitlines()[0].split("'")[1]
        print(f"`time -l` peak memory footprint {rec['time_peak_footprint_bytes'] / 1e9:.2f} GB, "
              f"wall {rec['wall_s']:.1f} s, {rec['n_samples']} samples, power: {power}\n")  # fmt: skip
        print("| phase | s | footprint at start GB | peak GB | rise over start GB |")
        print("|---|---:|---:|---:|---:|")
        for r in rows:
            star = "" if r["exact"] else " (sampled)"
            print(f"| {r['phase']} | {r['seconds']:.1f} | {r['start'] / 1e9:.2f} | "
                  f"{r['peak'] / 1e9:.2f}{star} | {(r['peak'] - r['start']) / 1e9:.2f} |")  # fmt: skip
        by = {r["phase"]: r for r in rows}
        m = next(x for x in markers if x["name"] == "assemble" and x["ev"] == "E")
        rs = next(x for x in markers if x["name"] == "resample" and x["ev"] == "B")
        store = next(x for x in markers if x["name"] == "refine_points" and x["ev"] == "B")["store_size"]
        mosaic_nodes, grid_nodes = m["shape"][0] * m["shape"][1], rs["rows"] * rs["cols"]
        canvas = grid_nodes * 4
        threads, block_rows = 10, rs.get("block_rows", 256)
        p1, p2, p3 = by["P1 decode + assemble"], by["P2 resample (blocks)"], by["P3 refine phase 1"]
        p4, p5 = by["P4 store fill + freeze"], by["P5 phase 2 (refine_points)"]
        p2e = by["P2 end (DemTile copy)"]
        tri1, tri = counts["phase1_triangles"], counts["final_triangles"]
        unit = {
            "source mosaic nodes": mosaic_nodes, "target grid nodes": grid_nodes,
            "check points": store, "phase 1 triangles (2 V - hull)": tri1, "final triangles": tri,
            "P1 rise per mosaic node, B": (p1["end"] - p1["start"]) / mosaic_nodes,
            "P2 transient per block node, B (peak - start - canvas) / (threads x block rows x cols)":
                (p2["peak"] - p2["start"] - canvas) / (threads * block_rows * rs["cols"]),
            "P2 end rise per grid node, B (copy)": (p2e["peak"] - p2["start"]) / grid_nodes,
            "P4 rise per check point, B": (p4["peak"] - p4["start"]) / store,
            "P4 held after freeze per check point, B": (p4["end"] - p4["start"]) / store,
            "P3 rise per phase-1 triangle, B": (p3["peak"] - p3["start"]) / tri1,
            "P5 rise per final triangle, B": (p5["peak"] - p5["start"]) / tri,
        }  # fmt: skip
        print("\n| measure | value |\n|---|---:|")
        for k, v in unit.items():
            print(f"| {k} | {v:,.1f} |" if isinstance(v, float) else f"| {k} | {v:,} |")
        (runs / f"{tag}.phases.json").write_text(json.dumps({"phases": rows, "units": unit}, indent=1))


if __name__ == "__main__":
    main()
