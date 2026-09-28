"""Where the 'features read' time goes for the catchment: the widened R-tree
query (ids only), the same with the geometry blobs fetched, the decode, and the
candidates' vertex counts; and the plain-box query for comparison. Uses the
branch's own source_region, layer_info and decode_geometry."""
import sqlite3, sys, time
from pathlib import Path
import shapely
from tin_engine.crs import parse_crs
from tin_engine.domain import read_domain
from tin_engine.io.geopackage import decode_geometry, layer_info, _q
from tin_engine.feature_input import source_region

dom = read_domain(Path(sys.argv[1])).to_crs("EPSG:25833")
for path, table in (sys.argv[2], sys.argv[3]), (sys.argv[4], sys.argv[5]):
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    layer = layer_info(conn, table)
    box = source_region(dom, parse_crs("EPSG:25833"), layer.crs).bounds
    scale = 1.0 if parse_crs(layer.crs).is_geographic else 100_000.0
    s = "(r.maxx - r.minx + r.maxy - r.miny)"; w = f"({s} * max(1.0, {s} / :scale))"
    where = (f" WHERE r.maxx + {w} >= :minx AND r.minx - {w} <= :maxx AND r.maxy + {w} >= :miny AND r.miny - {w} <= :maxy")
    plain = " WHERE r.maxx >= :minx AND r.minx <= :maxx AND r.maxy >= :miny AND r.miny <= :maxy"
    p = dict(zip(("minx", "miny", "maxx", "maxy"), box), scale=scale)
    rt = _q(layer.rtree)
    n_rtree = conn.execute(f"SELECT count(*) FROM {rt}").fetchone()[0]
    t = time.perf_counter(); ids = [r[0] for r in conn.execute(f"SELECT r.id FROM {rt} r {where}", p)]; t_ids = time.perf_counter() - t
    t = time.perf_counter(); pids = [r[0] for r in conn.execute(f"SELECT r.id FROM {rt} r {plain}", p)]; t_plain = time.perf_counter() - t
    t = time.perf_counter()
    blobs = conn.execute(f"SELECT t.{_q(layer.pk)}, t.{_q(layer.column)} FROM {_q(layer.table)} t JOIN {rt} r ON t.{_q(layer.pk)} = r.id {where} ORDER BY t.{_q(layer.pk)}", p).fetchall()
    t_fetch = time.perf_counter() - t
    t = time.perf_counter(); geoms = [decode_geometry(b)[1] for _, b in blobs]; t_dec = time.perf_counter() - t
    nv = shapely.get_num_coordinates(geoms)
    big = sorted(zip(nv, (f for f, _ in blobs)), reverse=True)[:5]
    far = sorted(set(ids) - set(pids))
    nv_far = int(sum(n for n, (f, _) in zip(nv, blobs) if f in set(far)))
    print(f"{Path(path).name}:{table} ({layer.crs}): R-tree rows {n_rtree}; query box {tuple(round(b) for b in box)}")
    print(f"  widened query: {len(ids)} rows, ids only {t_ids:.3f} s; plain box: {len(pids)} rows, {t_plain*1000:.1f} ms")
    print(f"  widened query with blobs: {t_fetch:.3f} s, {sum(len(b) for _, b in blobs)/1e6:.1f} MB; decode {t_dec:.3f} s")
    print(f"  candidate vertices {int(nv.sum())}; in the {len(far)} rows only the widening adds: {nv_far}; largest (vertices, fid): {big}")
