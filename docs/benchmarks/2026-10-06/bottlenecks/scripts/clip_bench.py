"""features clip, decomposed: the no-reprojection path of _Tally._take
(feature_input.py), timed per step with perf_counter, plus input sizes, plus
candidate rewrites checked for identical output (WKB bytes).
usage: clip_bench.py <catchment>"""
import sys, time, json
from collections import Counter
import numpy as np, shapely
from shapely.geometry import LineString, Polygon
import tin_engine.feature_input as fi
from tin_engine.domain import read_domain
from tin_engine.crs import parse_crs

c = sys.argv[1]
D = "/Users/skavhaug/projects/rasputin_data"
dom = read_domain(__import__("pathlib").Path(f"/Users/skavhaug/projects/rasputin_scratch/norway/{c}/{c}_outline_nve.geojson"))
dem = parse_crs("EPSG:25833")
shapely.prepare(dom.polygon)
t = time.perf_counter()
read = fi.read_source(__import__("pathlib").Path(f"{D}/corine2018_dtm10_utm33.gpkg"), "corine2018", "Code_18",
                      lambda own: fi.source_region(dom, dem, own).bounds)
t_read = time.perf_counter() - t
region = fi.source_region(dom, dem, dem)
shapely.prepare(region)
geoms = [shapely.force_2d(g) for _, g, _ in read.rows]
sizes = dict(features_read=len(geoms), domain_vertices=int(shapely.get_num_coordinates(dom.polygon)),
             region_vertices=int(shapely.get_num_coordinates(region)),
             vertices_read=int(sum(shapely.get_num_coordinates(g) for g in geoms)),
             rings_read=int(sum(shapely.get_num_interior_rings(p)+1 for g in geoms for p in shapely.get_parts(g))),
             read_seconds=round(t_read, 3))
# --- baseline, as _take does it, timed per step
T = Counter(); N = Counter(); out_base = []; vtx_kept = 0
for g in geoms:
    t0 = time.perf_counter(); kept = list(fi.pre_clip(g, region)); t1 = time.perf_counter()
    T["pre_clip"] += t1 - t0
    lines = []
    for line in kept:
        ta = time.perf_counter()
        covered = dom.polygon.covers(line)
        tb = time.perf_counter()
        parts = shapely.get_parts(shapely.intersection(line, dom.polygon))
        tc = time.perf_counter()
        T["covers_probe(not in prod)"] += tb - ta
        T["intersection: chain covered"if covered else "intersection: chain cut"] += tc - tb
        N["intersection: chain covered" if covered else "intersection: chain cut"] += 1
        lines += [p for p in parts if isinstance(p, LineString) and p.length > 0]
    t2 = time.perf_counter()
    if lines:
        N["kept"] += 1; vtx_kept += sum(len(l.coords) for l in lines)
    elif g.geom_type in ("Polygon","MultiPolygon") and g.intersects(dom.polygon.point_on_surface()):
        N["kept_covering"] += 1
    else:
        N["outside"] += 1
    T["after"] += time.perf_counter() - t2
    N["chains"] += len(kept)
    N["vertices_in_chains"] += sum(len(k.coords) for k in kept)
    out_base.append(b"".join(shapely.to_wkb(lines)) if lines else b"")
sizes.update(N); sizes["vertices_kept"] = vtx_kept
# pre_clip internals, uninstrumented
P = Counter(); allfalse = alltrue = 0
for g in geoms:
    for part in shapely.get_parts(g):
        rings = [(r, True) for r in (part.exterior, *part.interiors)] if isinstance(part, Polygon) else [(part, False)]
        for ring, closed in rings:
            t0 = time.perf_counter(); xy = shapely.get_coordinates(ring)
            seg = shapely.linestrings(np.stack([xy[:-1], xy[1:]], axis=1)); t1 = time.perf_counter()
            keep = np.asarray(shapely.intersects(seg, region)); t2 = time.perf_counter()
            fi._runs(xy, keep, closed); t3 = time.perf_counter()
            P["coords+linestrings"] += t1 - t0; P["intersects(segments, region)"] += t2 - t1; P["_runs"] += t3 - t2
            allfalse += not keep.any(); alltrue += keep.all()
sizes.update(rings_no_edge_kept=allfalse, rings_all_edges_kept=alltrue,
             segments_tested=int(sizes["vertices_read"] - sizes["rings_read"]))
# --- candidate A: pre_clip with an early exit for rings with no kept edge, and endpoint prefilter
def runs_fast(xy, keep, closed):
    if not keep.any():
        return []
    return fi._runs(xy, keep, closed)
def pre_clip_A(geometry, region):
    out = []
    for part in shapely.get_parts(geometry):
        rings = [(r, True) for r in (part.exterior, *part.interiors)] if isinstance(part, Polygon) else [(part, False)]
        for ring, closed in rings:
            xy = shapely.get_coordinates(ring)
            if len(xy) < 2: continue
            inside = shapely.intersects_xy(region, xy[:, 0], xy[:, 1])
            keep = inside[:-1] | inside[1:]
            rest = np.flatnonzero(~keep)
            if len(rest):  # exact test only where neither endpoint is in the region
                seg = shapely.linestrings(np.stack([xy[rest], xy[rest + 1]], axis=1))
                keep[rest] = shapely.intersects(seg, region)
            out += runs_fast(xy, keep, closed)
    return tuple(out)
# candidate B: A plus a whole-geometry bbox/region test before pre_clip
def run_variant(pc, skip_covered, whole_test):
    res = []; t0 = time.perf_counter()
    for g in geoms:
        if whole_test and not region.intersects(g):
            kept = []
        else:
            kept = list(pc(g, region))
        lines = []
        for line in kept:
            if skip_covered and dom.polygon.covers(line):
                lines.append(line); continue
            lines += [p for p in shapely.get_parts(shapely.intersection(line, dom.polygon)) if isinstance(p, LineString) and p.length > 0]
        if not lines and g.geom_type in ("Polygon","MultiPolygon"):
            g.intersects(dom.polygon.point_on_surface())
        res.append(b"".join(shapely.to_wkb(lines)) if lines else b"")
    return time.perf_counter() - t0, res
variants = {}
for name, args in {"baseline (prod pre_clip)": (fi.pre_clip, False, False),
                   "A: endpoint prefilter + early exit": (pre_clip_A, False, False),
                   "A + whole-feature region test": (pre_clip_A, False, True),
                   "A + skip intersection when domain covers chain": (pre_clip_A, True, False),
                   "A + whole test + skip covered": (pre_clip_A, True, True)}.items():
    secs = []
    for _ in range(3):
        s, res = run_variant(*args); secs.append(s)
    variants[name] = dict(median_s=round(sorted(secs)[1], 3), runs=[round(x, 3) for x in secs],
                          identical_wkb=res == out_base,
                          differing_features=sum(a != b for a, b in zip(res, out_base)))
print(json.dumps(dict(catchment=c, sizes=sizes,
      steps_s={k: round(v, 3) for k, v in T.items()},
      pre_clip_internals_s={k: round(v, 3) for k, v in P.items()}, variants=variants), indent=1, default=int))
