"""compare.py <master.vtk> <branch.vtk>: vertices, triangle set, codes."""
import sys
import numpy as np
from tin_engine.io.vtk_read import read_mesh

def load(p):
    m = read_mesh(open(p, "rb").read())
    v = np.asarray(m.vertices); t = np.asarray(m.triangles)
    code = np.asarray(m.land_cover_code) if m.land_cover_code is not None else None
    filled = getattr(m, "land_cover_filled", None)
    filled = np.asarray(filled).astype(bool) if filled is not None else None
    return v, t, code, filled

def keys(v, t):
    # vertex key: exact coordinates; triangle key: sorted corner vertex keys
    vk = {tuple(x): i for i, x in enumerate(map(tuple, v))}
    order = np.lexsort(v.T[::-1])
    rank = np.empty(len(v), int); rank[order] = np.arange(len(v))
    tk = np.sort(rank[t], axis=1)
    return vk, order, tk

vm, tm, cm, fm = load(sys.argv[1]); vb, tb, cb, fb = load(sys.argv[2])
print(f"vertices master {len(vm)} branch {len(vb)}; triangles master {len(tm)} branch {len(tb)}")
same_v = len(vm) == len(vb) and np.array_equal(vm[np.lexsort(vm.T[::-1])], vb[np.lexsort(vb.T[::-1])])
print("vertex set equal:", same_v, "; same order:", len(vm) == len(vb) and np.array_equal(vm, vb))
if not same_v:
    sys.exit(0)
_, _, km = keys(vm, tm); _, _, kb = keys(vb, tb)
om = np.lexsort(km.T[::-1]); ob = np.lexsort(kb.T[::-1])
same_t = np.array_equal(km[om], kb[ob])
print("triangle set equal (corners sorted):", same_t)
if not same_t or cm is None or cb is None:
    sys.exit(0)
cm, cb = cm[om], cb[ob]
fb = fb[ob] if fb is not None else np.zeros(len(cb), bool)
diff = cm != cb
print(f"master code 0: {(cm == 0).sum()}; branch code 0: {(cb == 0).sum()}")
print(f"codes differ: {diff.sum()}; where master != 0: {(diff & (cm != 0)).sum()}")
print(f"branch marked filled: {fb.sum()}; differ outside marked: {(diff & ~fb).sum()}; "
      f"marked with master 0: {(fb & (cm == 0)).sum()}; marked with a different nonzero master code: {(fb & diff & (cm != 0)).sum()}")
if (fb & diff).any():
    pairs = {}
    for a, b in zip(cm[fb & diff], cb[fb & diff]):
        pairs[(int(a), int(b))] = pairs.get((int(a), int(b)), 0) + 1
    print("marked and changed, (master, branch): count", pairs)
