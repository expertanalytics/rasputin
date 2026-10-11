# ruff: noqa: E501 -- a one-off probe; its SQL and report lines are long
"""Increment 36 probe: how long a GeoPackage of a mesh takes to write, and how
big it is, for the candidate writer (one transaction, batched executemany,
the R-tree built after the rows) against Ola's stopgap converter.

    .venv/bin/python docs/increments/36-probes/gpkg_write_probe.py MESH.vtk OUT.gpkg
        [--stopgap PATH_TO_vtk_to_gpkg.py] [--no-text] [--no-rtree] [--flat]
        [--rtree select|executemany] [--rtree-order fid|morton] [--cache-mb N]
        [--batch 100000] [--stream] [--per-batch] [--twice] [--vtk-time]

Reads a binary legacy VTK written by `rasputin mesh` (points, triangles,
constraint lines, the `land_cover_code` cell array, the `crs` field) with
numpy, then writes the `triangles` layer with the attributes section 4 of
`docs/increments/36-geopackage-output.md` names, a `lines` layer for the
constraint edges, and the R-tree, timing each phase. `--twice` writes the
file a second time and compares the two files' SHA-256 (determinism).
`--stopgap` runs the converter on the same mesh for the baseline figures.
Writes nothing but OUT.gpkg (and OUT.2.gpkg with `--twice`).
"""

import argparse
import hashlib
import os
import sqlite3
import struct
import subprocess
import sys
import time

import numpy as np

_DT = {
    "double": ">f8",
    "float": ">f4",
    "int": ">i4",
    "unsigned_int": ">u4",
    "unsigned_char": ">u1",
    "vtktypeint64": ">i8",
}


def read_vtk(path):
    """Points (n,3), triangles (t,3), edges (e,2), per-triangle code or None, crs text."""
    with open(path, "rb") as f:
        b = f.read()
    assert b.split(b"\n", 3)[2].strip() == b"BINARY", "binary VTK only"
    pos, out, n_lines, tris = 0, {}, 0, None

    def line():
        nonlocal pos
        e = b.index(b"\n", pos)
        s = b[pos:e].decode("latin-1")
        pos = e + 1
        return s

    def arr(n, typ):
        nonlocal pos
        a = np.frombuffer(b, _DT[typ], n, pos)
        pos += a.nbytes
        if pos < len(b) and b[pos : pos + 1] == b"\n":
            pos += 1
        return a

    line()
    line()
    line()
    while pos < len(b):
        h = line().split()
        if not h or h[0] == "DATASET":
            continue
        k = h[0]
        if k == "FIELD":
            for _ in range(int(h[2])):
                name, nc, nt, typ = line().split()
                n = int(nc) * int(nt)
                if typ == "string":
                    vals = []
                    for _ in range(n):
                        ln = b[pos]
                        if ln >> 6 == 3:
                            length, hd = ln & 0x3F, 1
                        elif ln >> 6 == 2:
                            length, hd = struct.unpack(">H", b[pos : pos + 2])[0] & 0x3FFF, 2
                        elif ln >> 6 == 1:
                            length, hd = struct.unpack(">I", b[pos : pos + 4])[0] & 0x3FFFFFFF, 4
                        else:
                            length, hd = struct.unpack(">Q", b[pos : pos + 8])[0], 8
                        vals.append(b[pos + hd : pos + hd + length].decode())
                        pos += hd + length
                    if b[pos : pos + 1] == b"\n":
                        pos += 1
                    out["field:" + name] = vals
                else:
                    out["field:" + name] = arr(n, typ)
        elif k == "POINTS":
            out["points"] = arr(3 * int(h[1]), h[2]).reshape(-1, 3).astype(float)
        elif k in ("LINES", "POLYGONS"):
            n, size = int(h[1]), int(h[2])
            a = arr(size, "int")
            if k == "LINES":
                n_lines = n
                out["edges"] = a.reshape(n, 3)[:, 1:].astype(np.int64)
            else:
                tris = a.reshape(n, 4)
                assert (tris[:, 0] == 3).all()
                out["tris"] = tris[:, 1:].astype(np.int64)
        elif k in ("POINT_DATA", "CELL_DATA"):
            out["section"] = k
        elif k == "SCALARS":
            line()
            n = len(out["points"]) if out["section"] == "POINT_DATA" else n_lines + len(out["tris"])
            out[out["section"] + ":" + h[1]] = arr(n, h[2])
        else:
            raise ValueError(f"unexpected VTK section {h!r}")
    code = out.get("CELL_DATA:land_cover_code")
    if code is not None:
        code = np.asarray(code)[n_lines:]
    return (
        out["points"],
        out["tris"],
        out.get("edges", np.zeros((0, 2), np.int64)),
        code,
        out.get("field:crs", [None])[0],
    )


NAMES = {
    0: ("No polygon; constraint lines", "#3c3c3c"),
    311: ("Broad-leaved forest", "#4f8f3f"),
    312: ("Coniferous forest", "#1f5a2e"),
    313: ("Mixed forest", "#357438"),
    512: ("Water bodies", "#3e7bb6"),
    322: ("Moors and heathland", "#a7c47f"),
}


def polygon_blobs(xyz, srs, flat):
    """GeoPackage blobs for PolygonZ (or Polygon with --flat) triangles, xyz (m,3,3)."""
    m, nd = len(xyz), (2 if flat else 3)
    env = np.stack(
        [xyz[:, :, 0].min(1), xyz[:, :, 0].max(1), xyz[:, :, 1].min(1), xyz[:, :, 1].max(1)], 1
    )
    dt = np.dtype(
        [
            ("gp", "S2"),
            ("v", "u1"),
            ("f", "u1"),
            ("srs", "<i4"),
            ("env", "<f8", 4),
            ("bo", "u1"),
            ("t", "<u4"),
            ("nr", "<u4"),
            ("np", "<u4"),
            ("pts", "<f8", 4 * nd),
        ]
    )
    rec = np.zeros(m, dt)
    rec["gp"] = b"GP"
    rec["f"] = 0b0011
    rec["srs"] = srs
    rec["env"] = env
    rec["bo"] = 1
    rec["t"] = 3 if flat else 1003
    rec["nr"] = 1
    rec["np"] = 4
    rec["pts"] = np.concatenate([xyz, xyz[:, :1]], 1)[:, :, :nd].reshape(m, 4 * nd)
    return rec.tobytes(), dt.itemsize, env


def line_blobs(xyz, srs):
    """LineStringZ blobs, xyz (e,2,3)."""
    e = len(xyz)
    env = np.stack(
        [xyz[:, :, 0].min(1), xyz[:, :, 0].max(1), xyz[:, :, 1].min(1), xyz[:, :, 1].max(1)], 1
    )
    dt = np.dtype(
        [
            ("gp", "S2"),
            ("v", "u1"),
            ("f", "u1"),
            ("srs", "<i4"),
            ("env", "<f8", 4),
            ("bo", "u1"),
            ("t", "<u4"),
            ("np", "<u4"),
            ("pts", "<f8", 6),
        ]
    )
    rec = np.zeros(e, dt)
    rec["gp"] = b"GP"
    rec["f"] = 0b0011
    rec["srs"] = srs
    rec["env"] = env
    rec["bo"] = 1
    rec["t"] = 1002
    rec["np"] = 2
    rec["pts"] = xyz.reshape(e, 6)
    return rec.tobytes(), dt.itemsize, env


def chunks(raw, w, m):
    view = memoryview(raw)
    return [view[i * w : (i + 1) * w] for i in range(m)]


def write(out, pts, tris, edges, code, srs, a):
    t = {}
    t0 = time.perf_counter()
    xyz = pts[tris]
    e1, e2 = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
    area = 0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
    nz = np.cross(e1, e2)
    slope = np.degrees(np.arctan2(np.hypot(nz[:, 0], nz[:, 1]), np.abs(nz[:, 2])))
    z = xyz[:, :, 2]
    if a.per_batch:  # no whole-mesh blob array: envelopes are filled per batch below
        raw, w, env = b"", 0, np.zeros((len(tris), 4))
        del xyz
    else:
        raw, w, env = polygon_blobs(xyz, srs, a.flat)
    blobs = [] if a.stream else chunks(raw, w, len(tris))
    cc = code if code is not None else np.zeros(len(tris), int)
    t["blobs"] = time.perf_counter() - t0

    if os.path.exists(out):
        os.remove(out)
    con = sqlite3.connect(out, isolation_level=None)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA application_id=1196444487")
    con.execute("PRAGMA user_version=10400")
    con.executescript("""
    CREATE TABLE gpkg_spatial_ref_sys (srs_name TEXT NOT NULL, srs_id INTEGER PRIMARY KEY, organization TEXT NOT NULL,
      organization_coordsys_id INTEGER NOT NULL, definition TEXT NOT NULL, description TEXT);
    CREATE TABLE gpkg_contents (table_name TEXT PRIMARY KEY, data_type TEXT NOT NULL, identifier TEXT UNIQUE, description TEXT DEFAULT '',
      last_change DATETIME NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')), min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE,
      srs_id INTEGER, CONSTRAINT fk_gc_r_srs_id FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys(srs_id));
    CREATE TABLE gpkg_geometry_columns (table_name TEXT NOT NULL, column_name TEXT NOT NULL, geometry_type_name TEXT NOT NULL,
      srs_id INTEGER NOT NULL, z TINYINT NOT NULL, m TINYINT NOT NULL, CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name));
    CREATE TABLE gpkg_extensions (table_name TEXT, column_name TEXT, extension_name TEXT NOT NULL, definition TEXT NOT NULL,
      scope TEXT NOT NULL, CONSTRAINT ge_tce UNIQUE (table_name, column_name, extension_name));
    INSERT INTO gpkg_spatial_ref_sys VALUES ('Undefined cartesian SRS', -1, 'NONE', -1, 'undefined', NULL),
      ('Undefined geographic SRS', 0, 'NONE', 0, 'undefined', NULL),
      ('WGS 84', 4326, 'EPSG', 4326, 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]', NULL);
    """)
    con.execute(
        "INSERT INTO gpkg_spatial_ref_sys VALUES (?,?,?,?,?,?)",
        (f"EPSG:{srs}", srs, "EPSG", srs, f"EPSG:{srs} (probe: WKT omitted)", None),
    )
    t1 = time.perf_counter()
    text = "" if a.no_text else ", land_cover TEXT, colour TEXT"
    con.execute(
        f"CREATE TABLE triangles (fid INTEGER PRIMARY KEY, geom POLYGON, code INTEGER{text}, "
        "z_min REAL, z_max REAL, z_mean REAL, area_m2 REAL, slope_deg REAL)"
    )
    if a.no_text:
        rows = zip(
            range(1, len(tris) + 1),
            blobs,
            cc.tolist(),
            z.min(1).tolist(),
            z.max(1).tolist(),
            z.mean(1).tolist(),
            area.tolist(),
            slope.tolist(),
            strict=False,
        )
        sql = "INSERT INTO triangles VALUES (?,?,?,?,?,?,?,?)"
    else:
        names = [NAMES.get(int(c), (str(c), "#ff00ff"))[0] for c in cc]
        cols = [NAMES.get(int(c), (str(c), "#ff00ff"))[1] for c in cc]
        rows = zip(
            range(1, len(tris) + 1),
            blobs,
            cc.tolist(),
            names,
            cols,
            z.min(1).tolist(),
            z.max(1).tolist(),
            z.mean(1).tolist(),
            area.tolist(),
            slope.tolist(),
            strict=False,
        )
        sql = "INSERT INTO triangles VALUES (?,?,?,?,?,?,?,?,?,?)"
    con.execute("BEGIN")
    if a.stream:  # one batch of row tuples alive at a time, built from array slices
        view = memoryview(raw)
        for i in range(0, len(tris), a.batch):
            j = min(i + a.batch, len(tris))
            sl = slice(i, j)
            if a.per_batch:  # the blobs of this batch only, from the vertex and triangle arrays
                raw_b, w, env[sl] = polygon_blobs(pts[tris[sl]], srs, a.flat)
                view = memoryview(raw_b)
                i0 = i
            else:
                i0 = 0
            cols_ = [
                range(i + 1, j + 1),
                [view[(k - i0) * w : (k - i0 + 1) * w] for k in range(i, j)],
                cc[sl].tolist(),
            ]
            if not a.no_text:
                cols_ += [
                    [NAMES.get(int(c), (str(c), "#ff00ff"))[0] for c in cc[sl]],
                    [NAMES.get(int(c), (str(c), "#ff00ff"))[1] for c in cc[sl]],
                ]
            cols_ += [
                z[sl].min(1).tolist(),
                z[sl].max(1).tolist(),
                z[sl].mean(1).tolist(),
                area[sl].tolist(),
                slope[sl].tolist(),
            ]
            con.executemany(sql, zip(*cols_, strict=True))
    else:
        rows = list(rows)
        for i in range(0, len(rows), a.batch):
            con.executemany(sql, rows[i : i + a.batch])
    con.execute("COMMIT")
    t["insert"] = time.perf_counter() - t1
    t2 = time.perf_counter()
    if not a.no_rtree:
        con.execute(
            "CREATE VIRTUAL TABLE rtree_triangles_geom USING rtree(id, minx, maxx, miny, maxy)"
        )
        if a.cache_mb:
            con.execute(f"PRAGMA cache_size=-{a.cache_mb * 1024}")
        order = np.arange(len(tris))
        if a.rtree_order == "morton":  # the R-tree fed in a space-filling-curve order
            q = (
                (env[:, [0, 2]] - env[:, [0, 2]].min(0))
                / max(np.ptp(env[:, [0, 2]], 0).max(), 1e-9)
                * 65535
            ).astype(np.uint64)
            key = np.zeros(len(tris), np.uint64)
            for bit in range(16):
                key |= ((q[:, 0] >> bit) & 1) << (2 * bit) | ((q[:, 1] >> bit) & 1) << (2 * bit + 1)
            order = np.argsort(key, kind="stable")

        def env_batches():  # (id, minx, maxx, miny, maxy) tuples, one batch alive at a time
            for i in range(0, len(tris), a.batch):
                ids = (order[i : i + a.batch] + 1).tolist()
                yield zip(ids, *env[order[i : i + a.batch]].T.tolist(), strict=True)

        con.execute("BEGIN")
        if a.rtree == "executemany":
            for batch in env_batches():
                con.executemany("INSERT INTO rtree_triangles_geom VALUES (?,?,?,?,?)", batch)
        else:
            con.execute(
                "CREATE TEMP TABLE env (id INTEGER PRIMARY KEY, minx REAL, maxx REAL, miny REAL, maxy REAL)"
            )
            for batch in env_batches():
                con.executemany("INSERT INTO env VALUES (?,?,?,?,?)", batch)
            con.execute("INSERT INTO rtree_triangles_geom SELECT * FROM env ORDER BY rowid")
            con.execute("DROP TABLE env")
        con.execute("COMMIT")
        con.execute(
            "INSERT INTO gpkg_extensions VALUES ('triangles','geom','gpkg_rtree_index',"
            "'http://www.geopackage.org/spec120/#extension_rtree','write-only')"
        )
    t["rtree"] = time.perf_counter() - t2
    bb = (env[:, 0].min(), env[:, 2].min(), env[:, 1].max(), env[:, 3].max())
    con.execute(
        "INSERT INTO gpkg_contents VALUES ('triangles','features','triangles','', '1970-01-01T00:00:00.000Z',?,?,?,?,?)",
        (*bb, srs),
    )
    con.execute(
        "INSERT INTO gpkg_geometry_columns VALUES ('triangles','geom','POLYGON',?,?,0)",
        (srs, 0 if a.flat else 1),
    )
    t3 = time.perf_counter()
    if len(edges):
        exyz = pts[edges]
        raw, w, env2 = line_blobs(exyz, srs)
        con.execute("CREATE TABLE lines (fid INTEGER PRIMARY KEY, geom LINESTRING, length_m REAL)")
        con.execute("BEGIN")
        con.executemany(
            "INSERT INTO lines VALUES (?,?,?)",
            zip(
                range(1, len(edges) + 1),
                chunks(raw, w, len(edges)),
                np.hypot(*(exyz[:, 1, :2] - exyz[:, 0, :2]).T).tolist(),
                strict=False,
            ),
        )
        con.execute("COMMIT")
        bb = (env2[:, 0].min(), env2[:, 2].min(), env2[:, 1].max(), env2[:, 3].max())
        con.execute(
            "INSERT INTO gpkg_contents VALUES ('lines','features','lines','', '1970-01-01T00:00:00.000Z',?,?,?,?,?)",
            (*bb, srs),
        )
        con.execute(
            "INSERT INTO gpkg_geometry_columns VALUES ('lines','geom','LINESTRING',?,1,0)", (srs,)
        )
    t["lines"] = time.perf_counter() - t3
    con.close()
    t["total"] = time.perf_counter() - t0
    return t


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("vtk")
    ap.add_argument("out")
    ap.add_argument("--stopgap")
    ap.add_argument("--no-text", action="store_true")
    ap.add_argument("--no-rtree", action="store_true")
    ap.add_argument("--flat", action="store_true")
    ap.add_argument("--rtree", default="select", choices=["select", "executemany"])
    ap.add_argument("--batch", type=int, default=100_000)
    ap.add_argument("--twice", action="store_true")
    ap.add_argument(
        "--stream",
        action="store_true",
        help="build each batch's row tuples from array slices (bounded memory)",
    )
    ap.add_argument(
        "--per-batch",
        action="store_true",
        help="with --stream: build the blobs per batch too (no whole-mesh blob array)",
    )
    ap.add_argument(
        "--vtk-time",
        action="store_true",
        help="also time write_vtk(binary=True) on the same arrays",
    )
    ap.add_argument("--rtree-order", default="fid", choices=["fid", "morton"])
    ap.add_argument(
        "--cache-mb",
        type=int,
        default=0,
        help="PRAGMA cache_size for the R-tree build, MB (0: SQLite's default)",
    )
    a = ap.parse_args()
    t0 = time.perf_counter()
    pts, tris, edges, code, crs = read_vtk(a.vtk)
    print(
        f"read {a.vtk}: {len(pts)} points, {len(tris)} triangles, {len(edges)} constraint edges, "
        f"crs {crs}, {time.perf_counter() - t0:.2f} s; vtk {os.path.getsize(a.vtk) / 1e6:.1f} MB"
    )
    srs = int(crs.split(":")[1]) if crs and crs.upper().startswith("EPSG:") else -1
    flags = " ".join(
        f
        for f in (
            "--no-text" * a.no_text,
            "--no-rtree" * a.no_rtree,
            "--flat" * a.flat,
            f"--rtree {a.rtree}",
            f"--rtree-order {a.rtree_order}",
            f"--cache-mb {a.cache_mb}" * bool(a.cache_mb),
            "--stream" * a.stream,
            "--per-batch" * a.per_batch,
            f"--batch {a.batch}",
        )
        if f
    )
    t = write(a.out, pts, tris, edges, code, srs, a)
    size = os.path.getsize(a.out) / 1e6
    print(
        f"candidate [{flags}]: "
        + ", ".join(f"{k} {v:.2f} s" for k, v in t.items())
        + f"; {size:.1f} MB"
    )
    if a.twice:
        out2 = a.out[:-5] + ".2.gpkg"
        write(out2, pts, tris, edges, code, srs, a)
        same = sha(a.out) == sha(out2)
        print(f"determinism: second write {'identical' if same else 'DIFFERS'} (sha256)")
        os.remove(out2)
    if a.vtk_time:  # the binary VTK write of the same mesh, for the multiple
        from tin_engine.features import DEFAULT_VOCABULARY
        from tin_engine.io.vtk_legacy import write_vtk

        t1 = time.perf_counter()
        data = write_vtk(
            pts,
            triangles=tris,
            edges=edges,
            edge_masks=np.zeros(len(edges), np.uint32),
            vocabulary=DEFAULT_VOCABULARY,
            binary=True,
            triangle_codes=code,
            land_cover_codes="probe" if code is not None else "",
        )
        t2 = time.perf_counter()
        with open(a.out[:-5] + ".probe.vtk", "wb") as f:
            f.write(data)
        print(
            f"binary vtk: encode {t2 - t1:.2f} s, disk {time.perf_counter() - t2:.2f} s, {len(data) / 1e6:.1f} MB"
        )
        os.remove(a.out[:-5] + ".probe.vtk")
    if a.stopgap:
        out3 = a.out[:-5] + ".stopgap.gpkg"
        t1 = time.perf_counter()
        r = subprocess.run([sys.executable, a.stopgap, a.vtk, out3], capture_output=True, text=True)
        print(
            f"stopgap: {time.perf_counter() - t1:.2f} s wall (its read included), {os.path.getsize(out3) / 1e6:.1f} MB; "
            f"last line: {r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-200:]}"
        )
        os.remove(out3)


if __name__ == "__main__":
    main()
