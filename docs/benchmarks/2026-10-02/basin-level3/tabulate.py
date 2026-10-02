"""The per-unit table of the level-3 runs. A measurement script.

    python tabulate.py RUNS_DIR UNITS_JSON MESH_DIR_761 MESH_DIR

Reads, per unit: area and box (target grid) nodes from UNITS_JSON; output
vertices, triangles, worst angle and max degree from ``<p>.stats.md``; the
final check's points, insertions and rounds from the mesh header's elevation
sentence (URL-encoded in the text field); wall time and peak footprint from
``/usr/bin/time -l`` in ``<p>.log``. Unit 761's mesh is in MESH_DIR_761 (made
earlier, same software); the rest in MESH_DIR. Prints a Markdown table and
writes ``RUNS_DIR/table.json``.
"""

import json
import re
import sys
from pathlib import Path
from urllib.parse import unquote


def cell(md: str, label: str) -> str:
    return re.search(rf"^\| {re.escape(label)} \| ([^|]+) \|", md, re.M).group(1).strip()


runs, units_json, mesh761, meshes = (Path(a) for a in sys.argv[1:5])
rows = []
for u in json.loads(units_json.read_text()):
    p = u["prefix"]
    md = (runs / f"{p}.stats.md").read_text()
    log = (runs / f"{p}.log").read_text(errors="replace")
    vtk = (mesh761 if p == "761" else meshes) / f"sub_basin_{p}_anadem_tol20m.vtk"
    head = unquote(vtk.read_bytes()[:400_000].decode("latin-1"))
    m = re.search(r"checked against (\d+) source nodes: (\d+) inserted in (\d+) rounds", head)
    quality = re.search(r"^\| minimum angle \| [^|]+ \| [^|]+ \| ([^|]+) \| ([^|]+) \|", md, re.M)
    degree = re.search(r"^\| vertex degree \(triangles\) \| [^|]+ \| [^|]+ \| (\d+) \|", md, re.M)
    rows.append({
        "unit": p, "area_km2": round(u["area_km2"], 1), "box_nodes": u["grid_nodes"],
        "check_points": int(m.group(1)), "final_check_inserted": int(m.group(2)),
        "final_check_rounds": int(m.group(3)),
        "vertices": int(cell(md, "output vertices")), "triangles": int(cell(md, "output triangles")),
        "under_10deg": quality.group(1).strip(), "worst_angle": quality.group(2).strip(),
        "max_degree": int(degree.group(1)),
        "peak_gb": int(re.search(r"(\d+)\s+peak memory footprint", log).group(1)) / 1e9,
        "wall_s": float(re.search(r"([\d.]+) real", log).group(1)),
        "vtk_mb": vtk.stat().st_size / 1e6,
    })  # fmt: skip
(runs / "table.json").write_text(json.dumps(rows, indent=1))
print("| unit | area km² | box nodes | check points | triangles | vertices | final-check inserted (rounds) | peak footprint GB | wall s | worst angle | < 10° | max degree | .vtk MB |")  # fmt: skip
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for r in rows:
    print(f"| {r['unit']} | {r['area_km2']:,.0f} | {r['box_nodes']:,} | {r['check_points']:,} | "
          f"{r['triangles']:,} | {r['vertices']:,} | {r['final_check_inserted']:,} ({r['final_check_rounds']}) | "
          f"{r['peak_gb']:.2f} | {r['wall_s']:.1f} | {r['worst_angle']} | {r['under_10deg']} | "
          f"{r['max_degree']} | {r['vtk_mb']:.1f} |")  # fmt: skip
t = {k: sum(r[k] for r in rows) for k in ("area_km2", "box_nodes", "check_points", "triangles",
                                          "vertices", "final_check_inserted", "wall_s", "vtk_mb")}  # fmt: skip
print(f"| **total** | {t['area_km2']:,.0f} | {t['box_nodes']:,} | {t['check_points']:,} | "
      f"**{t['triangles']:,}** | {t['vertices']:,} | {t['final_check_inserted']:,} | "
      f"max {max(r['peak_gb'] for r in rows):.2f} | {t['wall_s']:.1f} | | | | {t['vtk_mb']:.1f} |")  # fmt: skip
