"""One catchment run for 20c-2's acceptance (@perf, 2026-10-07; copied from ../20c-1/fix-69f37d1c/drive.py, tools path one level up).

``drive.py --pkg PKG --json OUT -- <rasputin mesh argv>``

Loads ``tin_engine`` from PKG (a ``tools/bench.py`` build-bench/pkg tree; the
editable finder is dropped, as bench.py's ``_child`` does), prints where
``tin_engine`` and ``_core`` came from, runs the CLI, and afterwards (outside
every timed phase: the ``--stats`` file is already written) measures the
trimmed mesh exactly: triangle count, the count under 1 degree, the worst
angle at full precision, the constrained Delaunay check (bench.py's), and a
hash of the arrays and of the output file from ``POINTS`` on. It also keeps
the ``--stats`` phase times unrounded (``PhaseClock.add``; the file rounds to
1 ms).
"""

import hashlib
import json
import sys
from pathlib import Path

argv = sys.argv[1:]
split = argv.index("--")
own, rasputin = argv[:split], argv[split + 1 :]
pkg = own[own.index("--pkg") + 1]
out_json = Path(own[own.index("--json") + 1])
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)

import numpy as np  # noqa: E402
import tin_engine  # noqa: E402
import tin_engine._core  # noqa: E402
import tin_engine.cli as cli  # noqa: E402
import tin_engine.stats  # noqa: E402

for module in (tin_engine, tin_engine._core):
    if not str(module.__file__).startswith(pkg):
        raise SystemExit(f"drive: {module.__file__} is not under {pkg}")
print(f"tin_engine.__file__ = {tin_engine.__file__}", file=sys.stderr)
print(f"_core.__file__ = {tin_engine._core.__file__}", file=sys.stderr)

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tools"))
import bench  # noqa: E402  (its _delaunay only)

seen: dict = {"phases": {}}
real_add = tin_engine.stats.PhaseClock.add


def add(self, name, seconds):  # the --stats phases at full precision (the file rounds to 1 ms)
    seen["phases"][name] = seen["phases"].get(name, 0.0) + seconds
    real_add(self, name, seconds)


tin_engine.stats.PhaseClock.add = add
real_trim = vars(cli)["trim"]


def trim(*args, **kwargs):  # the last trim call is the output mesh
    out = real_trim(*args, **kwargs)
    seen["trimmed"] = out
    return out


vars(cli)["trim"] = trim
cli.app(args=rasputin, prog_name="rasputin", standalone_mode=False)

t = seen["trimmed"]
xy = np.asarray(t.vertices, dtype=np.float64)
tris = np.asarray(t.triangles, dtype=np.int64).reshape(-1, 3)
cons = np.asarray(t.edges).astype(np.int64).reshape(-1, 2)
p = xy[:, :2][tris]
corners = []
for k in range(3):  # as tin_engine.stats.quality
    e = p[:, (k + 1) % 3] - p[:, k]
    f = p[:, (k + 2) % 3] - p[:, k]
    cross = np.abs(e[:, 0] * f[:, 1] - e[:, 1] * f[:, 0])
    corners.append(np.arctan2(cross, np.einsum("ij,ij->i", e, f)))
angles = np.degrees(np.min(corners, axis=0))
degree = np.bincount(tris.ravel(), minlength=len(xy))
checked, ambiguous, violations = bench._delaunay(xy[:, :2], tris, cons)
arrays = hashlib.sha256()
for a in (xy, tris, cons):
    arrays.update(np.ascontiguousarray(a).tobytes())
vtk = Path(rasputin[rasputin.index("--out") + 1]).read_bytes()
result = dict(
    tin_engine=tin_engine.__file__, core=tin_engine._core.__file__,
    hardening=getattr(tin_engine._core, "hardening", "none"),
    vertices=len(xy), triangles=len(tris), constraint_edges=len(cons),
    under_1=int(np.count_nonzero(angles < 1.0)), under_10=int(np.count_nonzero(angles < 10.0)),
    worst_angle=float(angles.min()), angle_median=float(np.median(angles)),
    max_degree=int(degree.max()),
    delaunay_checked=checked, delaunay_ambiguous=ambiguous, delaunay_violations=violations,
    arrays_sha256=arrays.hexdigest(),
    phases_s=seen["phases"],
    vtk_sha256=hashlib.sha256(vtk[vtk.index(b"\nPOINTS ") + 1 :]).hexdigest(),
)  # fmt: skip
out_json.write_text(json.dumps(result, indent=1) + "\n")
print("DRIVE " + json.dumps(result), file=sys.stderr)
