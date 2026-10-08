exec(open("angle2.py").read().split("print(\"under")[0])
lines = [shapely.LineString(P[e]) for e in m.edges]
ltree = shapely.STRtree(lines)
for t in np.argsort(mn)[:8]:
    tri = T[t]
    out = []
    for x in tri:
        pt = shapely.Point(P[x]); j = ltree.query_nearest(pt, return_distance=True)
        dl = float(j[1][0]); ds = min(shapely.distance(pt, rings[k]) for k in tree.query(pt, predicate="dwithin", distance=50)) if len(tree.query(pt, predicate="dwithin", distance=50)) else 99
        isv = (m.edges == x).any()
        out.append(f"({P[x][0]:.3f},{P[x][1]:.3f}) lineVertex={bool(isv)} d_line={dl:.4f} d_cover={ds:.4f}")
    print(f"{mn[t]:.4f}:", " | ".join(out))
