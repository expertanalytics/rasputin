"""Whether two mesh GeoPackages hold the same rows: every user table (not the
gpkg_* metadata), row by row, columns as stored; and, when not the same, the
largest difference of each real column. Run with the worktree's venv:
    .venv/bin/python docs/benchmarks/2026-10-26/106a/same_output.py A.gpkg B.gpkg
"""

import sqlite3
import sys


def tables(path: str) -> dict[str, list[tuple[object, ...]]]:
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    names = [
        n
        for (n,) in con.execute("select name from sqlite_master where type='table'")
        if not n.startswith(("gpkg_", "rtree_", "sqlite_"))
    ]
    out = {}
    for n in sorted(names):
        cols = [c[1] for c in con.execute(f'pragma table_info("{n}")')]
        out[n] = [cols, *con.execute(f'select * from "{n}" order by 1').fetchall()]
    return out


a, b = tables(sys.argv[1]), tables(sys.argv[2])
print("tables:", sorted(a) == sorted(b), sorted(a))
for name in sorted(set(a) & set(b)):
    ra, rb = a[name], b[name]
    if ra == rb:
        print(f"{name}: {len(ra) - 1} rows, identical")
        continue
    print(f"{name}: {len(ra) - 1} vs {len(rb) - 1} rows, columns equal: {ra[0] == rb[0]}")
    for i, col in enumerate(ra[0]):
        pairs = [(x[i], y[i]) for x, y in zip(ra[1:], rb[1:], strict=False)]
        diff = sum(x != y for x, y in pairs)
        if not diff:
            continue
        num = [abs(x - y) for x, y in pairs if isinstance(x, float) and isinstance(y, float)]
        print(f"  {col}: {diff} rows differ" + (f", max |diff| {max(num):.3g}" if num else ""))
