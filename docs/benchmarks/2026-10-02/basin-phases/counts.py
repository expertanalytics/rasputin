"""The counts a run printed, for `analyse.py`. A measurement script.

    python counts.py RUNS_DIR TAG MESH.vtk

Reads the CLI's report line from ``TAG.log`` (final triangles) and the
elevation sentence from the mesh file's header (phase 2's check points and
insertions), and writes ``TAG.counts.json``. Phase 1's triangle count is not
printed; it is taken as the final count less two per phase-2 insertion (an
insertion into a triangle or an edge adds two).
"""

import json
import re
import sys
from pathlib import Path

runs, tag, vtk = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
log = (runs / f"{tag}.log").read_text(errors="replace")
final = int(re.search(r"(\d+) triangles, achieved", log).group(1))
head = vtk.read_bytes()[:200_000].decode("latin-1")
m = re.search(r"checked against (\d+) source nodes: (\d+) inserted in (\d+) rounds", head)
points, inserted2, rounds2 = (int(g) for g in m.groups())
out = {"final_triangles": final, "check_points": points, "phase2_inserted": inserted2,
       "phase2_rounds": rounds2, "phase1_triangles": final - 2 * inserted2}  # fmt: skip
(runs / f"{tag}.counts.json").write_text(json.dumps(out, indent=1))
print(tag, out)
