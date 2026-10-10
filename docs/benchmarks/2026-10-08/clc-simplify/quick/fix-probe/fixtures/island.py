from fx import *
import shapely
A, D = (100.0, 0.0), (100.0, 100.0)
def winner(B, C):
    pl = placements(A, B, C, D); devs = [deviation([A, B, C, D], e) for _, e in pl]
    return pl[0 if devs[0] <= devs[1] else 1][1], pl[1 if devs[0] <= devs[1] else 0][1]
def island_ok(B, C, E, isl):
    old = shapely.LineString([A, B, C, D]); new = shapely.LineString([A, E, D])
    tri = shapely.Polygon(isl)
    swept = shapely.make_valid(shapely.Polygon([A, B, C, D, E]))
    sq = shapely.Polygon([(100, 0), (200, 0), (200, 100), (100, 100)]) if isl[0][0] > 100 else shapely.Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
    right = shapely.Polygon(rings(B, C)[1]); left = shapely.Polygon(rings(B, C)[0])
    host = right if right.contains(tri) else left if left.contains(tri) else None
    return tri.is_valid and host is not None and not tri.intersects(old) and not tri.intersects(new) and not any(swept.covers(shapely.Point(p)) for p in isl)
def make_vertex_island(E, P, Q, d, side):  # one island vertex d from segment P-Q's midpoint, on `side`
    m = ((P[0]+Q[0])/2, (P[1]+Q[1])/2); t = sub(Q, P); L = math.hypot(*t); n = (-t[1]/L*side, t[0]/L*side)
    v = add(m, mul(n, d)); return [v, add(v, add(mul(n, 6), mul(t, 4/L))), add(v, add(mul(n, 6), mul(t, -4/L)))]
def make_edge_island(E, n, d):
    c = add(E, mul(n, d)); p = (-n[1], n[0])
    return [add(c, mul(p, 3)), add(add(c, mul(n, 6)), (0, 0)), add(c, mul(p, -3))]
