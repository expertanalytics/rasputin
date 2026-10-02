"""Re-check of the 15c-2 acceptance's independent final check at 5 m on the
Velhas piece, after the crash fix: one ``--ascii`` run per CRS, then
``../15c-2-acceptance/run_geo.py``'s ``independent_check`` on it, and the
15 m-shift control on the first.

A measurement script, not production code; nothing imports it.

    python recheck_5m.py RASPUTIN CACHE OUTLINE WINDOW_TIF SCRATCH OUT_JSON CRS_FILE [CRS_FILE ...]

Each CRS_FILE holds one ``--out-crs`` value (a code or the pasted WKT).
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import urllib.parse
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "15c-2-acceptance"))
import run_geo  # noqa: E402


def main() -> None:
    rasputin, cache, outline, window, scratch, out_json = sys.argv[1:7]
    Path(scratch).mkdir(parents=True, exist_ok=True)
    results = []
    for k, crs_file in enumerate(sys.argv[7:]):
        crs = Path(crs_file).read_text().strip()
        vtk = Path(scratch) / f"velhas_t5_{k}.vtk"
        cmd = [rasputin, "mesh", "--dem", "anadem-v1", "--cache", cache, "--domain", outline,
               "--out-crs", crs, "--tolerance", "5", "--ascii", "--out", str(vtk)]  # fmt: skip
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit(p.stderr[-2000:])
        head = urllib.parse.unquote(vtk.read_text(errors="replace")[:20000])
        h = float(re.search(r"square ([\d.]+) m grid", head).group(1))
        m = re.search(r"checked against (\d+) source nodes: (\d+) inserted in (\d+) rounds, "
                      r"max error ([^ ]+) m", head)  # fmt: skip
        rec = {"crs_file": Path(crs_file).name, "h": h,
               "phase2": dict(zip(("source_nodes", "inserted", "rounds", "max_error"), m.groups())),
               "check": run_geo.independent_check(vtk, Path(window), Path(outline), crs, h, 5.0)}  # fmt: skip
        if k == 0:
            rec["control_shift15"] = run_geo.independent_check(vtk, Path(window), Path(outline), crs, h, 5.0, 15.0)
        results.append(rec)
        Path(out_json).write_text(json.dumps(results, indent=1))
        c = rec["check"]
        print(crs_file, "interior over", c["interior"]["over_tolerance"], "of", c["interior"]["nodes"],
              "strip over", c["strip"]["over_tolerance"], "of", c["strip"]["nodes"], flush=True)  # fmt: skip


if __name__ == "__main__":
    main()
