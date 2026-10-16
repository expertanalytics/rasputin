import pickle, numpy as np, shapely
simp = pickle.load(open("band50.pkl", "rb"))[0]; src = pickle.load(open("clip_nosimp.pkl", "rb"))[0]
def table(ps):
    rows = []
    for pi, p in enumerate(ps):
        for q in shapely.get_parts(p):
            for r in shapely.get_rings(q):
                xy = shapely.get_coordinates(r)[:-1]
                rows += [(x, y, pi, k, len(xy)) for k, (x, y) in enumerate(xy)]
    return np.array(rows)
S, R = table(simp), table(src)
for cx, cy in [(145339.07, 6720531.47), (167389.45, 6667384.31), (146262.30, 6713029.95), (226633.86, 6578900.82), (152264.0, 6714980.8)]:
    print(f"near ({cx},{cy}):")
    for name, A in (("simplified", S), ("source", R)):
        d = np.hypot(A[:, 0] - cx, A[:, 1] - cy); sel = np.argsort(d)[:4]
        for i in sel:
            if d[i] < 2: print(f"   {name}: ({A[i,0]:.4f},{A[i,1]:.4f}) poly {int(A[i,2])} ring-vertex {int(A[i,3])}/{int(A[i,4])}")
