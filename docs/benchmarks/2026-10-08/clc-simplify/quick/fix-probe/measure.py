import pickle, time, sys, numpy as np, shapely
import tin_engine.feature_input as fi
for name in sys.argv[1:]:
    polys, outline, d = pickle.load(open(name, "rb"))
    t = time.perf_counter(); r = fi.snap_to_outline(polys, outline, d); t = time.perf_counter() - t
    bnd = outline.boundary
    on = near = edges = verts = 0
    for p in polys:
        for q in shapely.get_parts(p):
            for ring in shapely.get_rings(q):
                xy = shapely.get_coordinates(ring)
                segs = shapely.linestrings(np.stack([xy[:-1], xy[1:]], 1))
                dist = shapely.hausdorff_distance(segs, bnd) if False else None
                dd = shapely.distance(segs, bnd)
                # an edge is on the outline if both its ends and its midpoint are within 1e-6 of it
                mid = shapely.points((xy[:-1] + xy[1:]) / 2)
                onm = (shapely.distance(shapely.points(xy[:-1]), bnd) < 1e-6) & (shapely.distance(shapely.points(xy[1:]), bnd) < 1e-6) & (shapely.distance(mid, bnd) < 1e-6)
                edges += len(segs); on += int(onm.sum()); near += int(((dd <= d) & ~onm).sum())
    print(f"{name}: polys {len(polys)} edges {edges} on-outline {on} near-not-on {near} snap {t:.2f}s area_changed {r.area_changed:.0f} m2")
