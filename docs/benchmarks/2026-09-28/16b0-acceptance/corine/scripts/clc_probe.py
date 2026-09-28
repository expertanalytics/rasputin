"""16b probe: CORINE 2018 over a DTM10 tile, GDAL-free. Scratch only."""
import sqlite3, struct, sys, time, collections
from pathlib import Path
import numpy as np, shapely, shapely.wkb
from shapely.geometry import box, Polygon, MultiPolygon
from tin_engine.io.repository import TiffDemRepository
from tin_engine.crs import reprojector

GPKG = Path("../rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg")
T = "U2018_CLC2018_V2020_20u1"

def gpkg_wkb(blob):
    assert blob[:2] == b"GP", blob[:2]
    flags = blob[3]
    env = (flags >> 1) & 7
    n = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]
    assert not (flags >> 5) & 1, "extended"
    return bytes(blob[8 + n:]), flags

def rect_of(meta):
    x0 = meta.x_min; x1 = meta.x_min + (meta.cols - 1) * meta.delta_x
    y1 = meta.y_max; y0 = meta.y_max - (meta.rows - 1) * meta.delta_y
    return x0, y0, x1, y1

def probe(tile_path, sub=None):
    repo = TiffDemRepository([tile_path]); fp = repo.footprints()[0]; meta = fp.meta
    x0, y0, x1, y1 = rect_of(meta)
    if sub: x0, y0, x1, y1 = sub
    rect = box(x0, y0, x1, y1)
    # rectangle boundary densified to 1 km, into 3035, bbox
    ring = np.asarray(rect.segmentize(1000.0).exterior.coords)
    to3035 = reprojector("EPSG:25833", "EPSG:3035"); back = reprojector("EPSG:3035", "EPSG:25833")
    r = to3035(ring); bx0, by0 = r.min(0); bx1, by1 = r.max(0)
    c = sqlite3.connect(f"file:{GPKG}?mode=ro", uri=True)
    t = time.perf_counter()
    rows = c.execute(f"select f.OBJECTID, f.Code_18, f.Shape from {T} f join rtree_{T}_Shape r on f.OBJECTID = r.id "
                     "where r.maxx >= ? and r.minx <= ? and r.maxy >= ? and r.miny <= ?", (bx0, bx1, by0, by1)).fetchall()
    tq = time.perf_counter() - t
    t = time.perf_counter()
    feats = []; flagset = collections.Counter(); raw_vertices = 0; invalid = 0; parts = 0; holes = 0
    for oid, code, blob in rows:
        wkb, fl = gpkg_wkb(blob); flagset[fl] += 1
        g = shapely.wkb.loads(wkb)
        raw_vertices += shapely.get_num_coordinates(g)
        polys = list(g.geoms) if isinstance(g, MultiPolygon) else [g]
        parts += len(polys); holes += sum(len(p.interiors) for p in polys)
        moved = []
        for p in polys:
            rs = [back(np.asarray(rr.coords)) for rr in (p.exterior, *p.interiors)]
            moved.append(Polygon(rs[0], rs[1:]))
        mg = MultiPolygon(moved)
        if not mg.is_valid: invalid += 1
        feats.append((oid, code, mg))
    tp = time.perf_counter() - t
    # clip
    clipped = []
    for oid, code, g in feats:
        k = g.intersection(rect)
        if not k.is_empty: clipped.append((oid, code, k))
    print(f"{tile_path.name} rect {x0:.0f},{y0:.0f} .. {x1:.0f},{y1:.0f}")
    print(f"rtree query {tq*1e3:.0f} ms, {len(rows)} candidate rows; parse+reproject {tp*1e3:.0f} ms; blob flags {dict(flagset)}")
    print(f"raw: {raw_vertices} vertices, {parts} parts, {holes} holes, invalid after reprojection {invalid}")
    print(f"after clip: {len(clipped)} features; classes {collections.Counter(c for _, c, _ in clipped).most_common()}")
    # segments of clipped geometry
    segs = []; vcount = 0
    geomtypes = collections.Counter()
    for oid, code, k in clipped:
        geomtypes[k.geom_type] += 1
        for p in getattr(k, "geoms", [k]):
            if p.geom_type != "Polygon": continue
            for rr in (p.exterior, *p.interiors):
                a = np.asarray(rr.coords); vcount += len(a) - 1
                for i in range(len(a) - 1): segs.append((tuple(a[i]), tuple(a[i+1]), code))
    print(f"clipped geometry types {dict(geomtypes)}; ring vertices {vcount}; segments {len(segs)}")
    L = np.array([np.hypot(b[0]-a[0], b[1]-a[1]) for a, b, _ in segs])
    print(f"segment length m: min {L.min():.3f} p1 {np.percentile(L,1):.2f} median {np.median(L):.1f} max {L.max():.0f}; "
          f"<1m {np.sum(L<1)}, <{meta.delta_x:g}m (a cell) {np.sum(L<meta.delta_x)}")
    # sharing: undirected key counts
    key = collections.Counter(); dirk = collections.Counter(); samecode = 0
    codes_of = collections.defaultdict(set)
    for a, b, code in segs:
        k = (min(a, b), max(a, b)); key[k] += 1; dirk[(a, b)] += 1; codes_of[k].add(code)
    mult = collections.Counter(key.values())
    same_dir = sum(1 for k, v in dirk.items() if v > 1)
    both = sum(1 for k, v in key.items() if v == 2)
    samecode = sum(1 for k, v in key.items() if v == 2 and len(codes_of[k]) == 1)
    print(f"undirected segment multiplicity {dict(mult)}; same-direction duplicates {same_dir}; shared by 2 with same code {samecode}")
    # on the rectangle border
    onb = sum(1 for (a, b) in key if (a[0] == b[0] and a[0] in (x0, x1)) or (a[1] == b[1] and a[1] in (y0, y1)))
    print(f"unique segments {len(key)}; on the clip rectangle's border {onb}")
    # T-junction check: vertices of one polygon lying in the interior of another's segment (not shared vertices)
    verts = set(); [verts.update([a, b]) for a, b, _ in segs]
    lines = shapely.linestrings([[a, b] for (a, b) in key if key[(a, b)] == 1])
    tree = shapely.STRtree(lines)
    pts = shapely.points(np.array(list(verts)))
    idx = tree.query(pts, predicate="dwithin", distance=1e-3)
    tj = 0; tjd = []
    for pi, li in idx.T:
        ln = lines[li]; p = pts[pi]
        c0 = shapely.get_coordinates(ln)
        pc = tuple(shapely.get_coordinates(p)[0])
        if pc == tuple(c0[0]) or pc == tuple(c0[1]): continue
        tj += 1; tjd.append(ln.distance(p))
    print(f"vertices within 1 mm of an unshared segment's interior (T-junction or near-miss): {tj}" + (f", max dist {max(tjd):.2e} m" if tjd else ""))
    # are vertices on DEM nodes?
    xy = np.array(list(verts))
    col = (xy[:,0]-meta.x_min)/meta.delta_x; row = (meta.y_max-xy[:,1])/meta.delta_y
    print(f"distinct vertices {len(xy)}; exactly on a node {int(np.sum((col==np.round(col))&(row==np.round(row))))}")
    # nearest vertex-to-vertex distance (distinct)
    P = shapely.points(xy); tr = shapely.STRtree(P)
    ii = tr.query(P, predicate="dwithin", distance=1.0)
    ii = ii[:, ii[0] < ii[1]]
    dd = np.hypot(*(xy[ii[0]] - xy[ii[1]]).T) if ii.size else np.array([])
    print(f"distinct vertex pairs closer than 1 m: {ii.shape[1]}" + (f", min {dd.min():.3f} m, under 1 cm {int(np.sum(dd<0.01))}" if dd.size else ""))
    U = shapely.union_all([k for _, _, k in clipped])
    print(f"coverage of the rectangle: {U.area/rect.area*100:.3f} %; uncovered area {(rect.area-U.area)/1e6:.3f} km2")
    return meta, clipped, rect

if __name__ == "__main__":
    probe(Path(sys.argv[1]))
