"""The branch's candidate query (JOIN, planned as SCAN t) against the same rows
fetched by primary key from the R-tree's ids (IN subquery): time, rows, plan.
Diagnostic only; the branch's code is not changed."""
import sqlite3, sys, time
EU = ("/Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg", "U2018_CLC2018_V2020_20u1", "OBJECTID", "Shape")
NO = ("/Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg", "corine2018", "fid", "geom")
BOX = {"U2018_CLC2018_V2020_20u1": (4390346, 4082024, 4410933, 4103082), "corine2018": (290173, 6637964, 310429, 6659036)}
for path, table, pk, col in (EU, NO):
    c = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    rt = f'"rtree_{table}_{col}"'
    s = "(r.maxx - r.minx + r.maxy - r.miny)"; w = f"({s} * max(1.0, {s} / :scale))"
    cond = f"r.maxx + {w} >= :minx AND r.minx - {w} <= :maxx AND r.maxy + {w} >= :miny AND r.miny - {w} <= :maxy"
    p = dict(zip(("minx", "miny", "maxx", "maxy"), BOX[table]), scale=1e5)
    join = f'SELECT t."{pk}", t."{col}", t."Code_18" FROM "{table}" t JOIN {rt} r ON t."{pk}" = r.id WHERE {cond} ORDER BY t."{pk}"'
    sub = f'SELECT t."{pk}", t."{col}", t."Code_18" FROM "{table}" t WHERE t."{pk}" IN (SELECT r.id FROM {rt} r WHERE {cond}) ORDER BY t."{pk}"'
    for name, q in (("join (as built)", join), ("IN subquery", sub)):
        plan = "; ".join(r[3] for r in c.execute("EXPLAIN QUERY PLAN " + q, p))
        t = time.perf_counter(); rows = c.execute(q, p).fetchall(); dt = time.perf_counter() - t
        print(f"{table} {name}: {len(rows)} rows, {sum(len(r[1]) for r in rows)/1e6:.1f} MB, {dt:.3f} s; plan: {plan}")
