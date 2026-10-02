"""The per-unit table of the level-3 runs at one tolerance. A measurement script.

    python tabulate.py RUNS_DIR UNITS_JSON TOLERANCE MESH_DIR...

Reads, per unit: area and box (target grid) nodes from UNITS_JSON; output
vertices, triangles, worst angle and max degree from ``<p>.stats.md``; the
final check's points, insertions and rounds from the mesh header's elevation
sentence (URL-encoded in the text field); wall time and peak footprint from
``/usr/bin/time -l`` in ``<p>.log``. The mesh
``sub_basin_<p>_anadem_tol<TOLERANCE>m.vtk`` is looked up in each MESH_DIR in
turn. Prints a Markdown table and writes ``RUNS_DIR/table.{json,md}``.
"""

import json
import re
import sys
from pathlib import Path
from urllib.parse import unquote

HEADER = (
    "| unit | area km² | box nodes | check points | triangles | vertices "
    "| final-check inserted (rounds) | peak footprint GB | wall s | worst angle "
    "| < 10° | max degree | .vtk MB | format |\n"
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"
)
SUMMED = ("area_km2", "box_nodes", "check_points", "triangles", "vertices")
SUMMED += ("final_check_inserted", "wall_s", "vtk_mb")


def cell(md: str, label: str) -> str:
    return re.search(rf"^\| {re.escape(label)} \| ([^|]+) \|", md, re.M).group(1).strip()


def find_mesh(dirs: list[Path], name: str) -> Path:
    return next(d / name for d in dirs if (d / name).exists())


def row(runs: Path, unit: dict, tol: str, dirs: list[Path]) -> dict:
    p = unit["prefix"]
    md = (runs / f"{p}.stats.md").read_text()
    log = (runs / f"{p}.log").read_text(errors="replace")
    vtk = find_mesh(dirs, f"sub_basin_{p}_anadem_tol{tol}m.vtk")
    raw = vtk.read_bytes()[:400_000].decode("latin-1")
    head = unquote(raw)
    m = re.search(r"checked against (\d+) source nodes: (\d+) inserted in (\d+) rounds", head)
    angle = re.search(r"^\| minimum angle \| [^|]+ \| [^|]+ \| ([^|]+) \| ([^|]+) \|", md, re.M)
    degree = re.search(r"^\| vertex degree \(triangles\) \| [^|]+ \| [^|]+ \| (\d+) \|", md, re.M)
    peak = re.search(r"(\d+)\s+peak memory footprint", log)
    return {
        "unit": p,
        "area_km2": round(unit["area_km2"], 1),
        "box_nodes": unit["grid_nodes"],
        "check_points": int(m.group(1)),
        "final_check_inserted": int(m.group(2)),
        "final_check_rounds": int(m.group(3)),
        "vertices": int(cell(md, "output vertices")),
        "triangles": int(cell(md, "output triangles")),
        "under_10deg": angle.group(1).strip(),
        "worst_angle": angle.group(2).strip(),
        "max_degree": int(degree.group(1)),
        "achieved_m": cell(md, f"{tol} m"),
        "peak_gb": int(peak.group(1)) / 1e9,
        "wall_s": float(re.search(r"([\d.]+) real", log).group(1)),
        "vtk_mb": vtk.stat().st_size / 1e6,
        "format": raw.splitlines()[2].strip().lower(),
        "mesh": str(vtk),
    }


def line(r: dict) -> str:
    cells = [
        r["unit"],
        f"{r['area_km2']:,.0f}",
        f"{r['box_nodes']:,}",
        f"{r['check_points']:,}",
        f"{r['triangles']:,}",
        f"{r['vertices']:,}",
        f"{r['final_check_inserted']:,} ({r['final_check_rounds']})",
        f"{r['peak_gb']:.2f}",
        f"{r['wall_s']:.1f}",
        r["worst_angle"],
        r["under_10deg"],
        str(r["max_degree"]),
        f"{r['vtk_mb']:.1f}",
        r["format"],
    ]
    return "| " + " | ".join(cells) + " |"


def main() -> None:
    runs, units_json, tol = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    dirs = [Path(a) for a in sys.argv[4:]]
    units = json.loads(units_json.read_text())
    rows = [row(runs, u, tol, dirs) for u in units if (runs / f"{u['prefix']}.stats.md").exists()]
    t = {k: sum(r[k] for r in rows) for k in SUMMED}
    total = [
        "**total**",
        f"{t['area_km2']:,.0f}",
        f"{t['box_nodes']:,}",
        f"{t['check_points']:,}",
        f"**{t['triangles']:,}**",
        f"{t['vertices']:,}",
        f"{t['final_check_inserted']:,}",
        f"max {max(r['peak_gb'] for r in rows):.2f}",
        f"{t['wall_s']:.1f}",
        "",
        "",
        "",
        f"{t['vtk_mb']:.1f}",
        "",
    ]
    text = "\n".join([HEADER, *map(line, rows), "| " + " | ".join(total) + " |"])
    (runs / "table.json").write_text(json.dumps(rows, indent=1))
    (runs / "table.md").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
