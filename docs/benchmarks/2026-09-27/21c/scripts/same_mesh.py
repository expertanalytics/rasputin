"""same_mesh.py A.vtk B.vtk -- are two ASCII VTK meshes the same triangulation up
to vertex and triangle numbering? Compares the sets of triangles as sets of
coordinate triples (bit-exact floats, as written)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[5] / "tools"))
import bench  # noqa: E402


def canon(path):
    m = bench.read_vtk_ascii(Path(path))
    pts = [tuple(p) for p in m.points.tolist()]
    return {frozenset(pts[i] for i in t) for t in m.triangles.tolist()}, len(m.triangles)


a, na = canon(sys.argv[1])
b, nb = canon(sys.argv[2])
print(f"triangles {na} / {nb}; distinct {len(a)} / {len(b)}; only in A {len(a - b)}; only in B {len(b - a)}; "
      f"{'SAME' if a == b else 'DIFFERENT'}")
