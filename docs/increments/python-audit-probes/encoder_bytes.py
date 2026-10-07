"""The bytes every mesh encoder writes, and every refusal it gives
(python-audit-pr-d.md, "The probe").

Run against one revision's `tin_engine`, imported from a scratch copy:

    python3 tools/scratch_copy.py <rev> <dir>
    PYTHONPATH=<that value> .venv/bin/python \\
        docs/increments/python-audit-probes/encoder_bytes.py <dir>

scratch_copy.py prints one line on stdout, the pytest command for the copy
(`cd <dir> && PYTHONPATH=<value> <python> -c ... tests/python/`); take
<value> from that line only (stderr may carry a `scratch_copy:` warning).
The probe asserts that `tin_engine` was imported from `<dir>/src_python`,
in its own process and in each `rasputin` it starts, so a run cannot
measure the work tree by mistake. It writes its inputs to a temporary
directory and removes it.

One line per case: `case | outcome`. The outcome is the sha256 of each
file written (the temporary directory's path replaced by `<dir>` first, so
two runs agree), or `refused: <message>`, or `UNCAUGHT <class>: <message>`.
For a command, a refusal is its last stderr line and its exit code.

Part 1 calls the writers directly, with arguments whose meaning PR D does
not change: `write_ply` (surface and edge files, text and binary, with and
without codes), `write_vtk` (text and binary, codes, fields), on an empty
mesh, one triangle, two triangles with constraint edges, a NaN z, and
coordinates at a UTM easting; with comments and fields that are plain,
non-ASCII, carry each control character class, `%` and spaces, and with
codes at and past the int32 bounds, of the wrong length, and NaN in a
float array; and `crs_label` on CRSs with and without an EPSG code, one
named in non-ASCII.

Part 2 runs `rasputin mesh` in a child process: a gallery fixture with
`--flat` to `.vtk` and `.ply` (with `--out-edges`), text and binary, with
`--crs` plain, non-ASCII, with `\\r`, and empty; the refusals for an unknown
suffix and for `--out-edges` beside `.vtk`; and a DEM run (a hand-written
21 x 17 GeoTIFF, named in non-ASCII) with CORINE-coded `--features`, also
named in non-ASCII, to `.vtk` and `.ply` with `--out-edges`, so the land-cover
codes, their comment and the escaped file names are in the bytes.
"""

from __future__ import annotations

import hashlib
import json
import os
import struct
import subprocess
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path

import numpy as np

X, Y = 500_000.0, 6_600_000.0
ROWS, COLS = 17, 21
#: Shown for a command's files that hold them, so a hash is seen to cover them.
MARKS = (b"land_cover_code", b"d\\xe9m.tif")
#: The `--record` entries that hold a file name, escaped to ASCII.
RECORD = ("dem_source", "dem_tiles", "features", "land_cover_codes")
CHILD = (
    "import sys, tin_engine; from pathlib import Path\n"
    "assert Path(tin_engine.__file__).is_relative_to(Path(sys.argv[1]) / 'src_python'), "
    "tin_engine.__file__\n"
    "from tin_engine.cli import app; app(sys.argv[2:], prog_name='rasputin')\n"
)


def outcome(call: Callable[..., bytes], *args: object, **kwargs: object) -> str:
    """The sha256 of `call(*args, **kwargs)`'s bytes, or what it raised."""
    try:
        data = call(*args, **kwargs)
    except ValueError as exc:
        return f"refused: {exc}"
    except Exception as exc:  # a crash is what the probe reports
        return f"UNCAUGHT {type(exc).__name__}: {exc}"
    return hashlib.sha256(data).hexdigest()[:16]


def meshes() -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """name -> (vertices, triangles, edges, masks)."""
    v4 = np.array([[0, 0, 1.5], [1, 0, 2.25], [1, 1, 3.0], [0, 1, 0.1]], dtype=float)
    two = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.uint32)
    e2 = np.array([[0, 1], [1, 2]], dtype=np.uint32)
    nan = v4.copy()
    nan[2, 2] = np.nan
    utm = v4 * 10.0 + np.array([X + 0.001, Y - 0.003, 1000.0])
    return {
        "empty": (np.zeros((0, 3)), np.zeros((0, 3), np.uint32), np.zeros((0, 2), np.uint32),
                  np.zeros(0, np.uint32)),
        "one triangle": (v4[:3], two[:1], np.zeros((0, 2), np.uint32), np.zeros(0, np.uint32)),
        "two with edges": (v4, two, e2, np.array([1, 3], np.uint32)),
        "nan z": (nan, two, e2, np.array([0, 2], np.uint32)),
        "utm": (utm, two, e2, np.array([1, 0], np.uint32)),
    }  # fmt: skip


TEXTS = {
    "plain": "crs EPSG:25833",
    "non-ascii": "crs ETRS89 \N{DEGREE SIGN}N",
    "cr": "crs x\rcomment forged",
    "tab": "crs a\tb",
    "del": "crs a\x7fb",
    "nul": "crs a\x00b",
    "percent space": "crs a%20b c",
    "empty": "",
}
CODES = {
    "none": None,
    "small": [1, 312],
    "int32 bounds": [2**31 - 1, -(2**31)],
    "past int32": [2**31, 0],
    "below int32": [-(2**31) - 1, 0],
    "float nan": [float("nan"), 1.0],
    "one short": [1],
}


def part1() -> list[str]:
    from tin_engine.crs import crs_label
    from tin_engine.features import DEFAULT_VOCABULARY as VOC
    from tin_engine.io.ply import write_ply
    from tin_engine.io.vtk_legacy import write_vtk

    out = []
    for m, (v, t, e, k) in meshes().items():
        for binary in (False, True):
            mode = "binary" if binary else "text"
            out.append(f"ply surface | {m} | {mode} | "
                       + outcome(write_ply, v, faces=t, ascii=not binary))  # fmt: skip
            out.append(f"ply edges | {m} | {mode} | " + outcome(
                write_ply, v, edges=e, edge_properties=k, ascii=not binary,
                vocabulary=VOC))  # fmt: skip
            out.append(f"vtk | {m} | {mode} | " + outcome(
                write_vtk, v, triangles=t, edges=e, edge_masks=k, vocabulary=VOC,
                binary=binary))  # fmt: skip
    v, t, e, k = meshes()["two with edges"]
    for name, text in TEXTS.items():
        for binary in (False, True):
            mode = "binary" if binary else "text"
            out.append(f"ply comment | {name} | {mode} | " + outcome(
                write_ply, v, faces=t, ascii=not binary, comments=[text]))  # fmt: skip
            out.append(f"vtk field | {name} | {mode} | " + outcome(
                write_vtk, v, triangles=t, edges=e, edge_masks=k, vocabulary=VOC, binary=binary,
                fields=[("crs", text)]))  # fmt: skip
    for name, codes in CODES.items():
        for binary in (False, True):
            mode = "binary" if binary else "text"
            c = None if codes is None else np.array(codes)
            out.append(f"ply codes | {name} | {mode} | " + outcome(
                write_ply, v, faces=t, ascii=not binary, face_codes=c,
                comments=["land_cover_codes x"]))  # fmt: skip
            out.append(f"vtk codes | {name} | {mode} | " + outcome(
                write_vtk, v, triangles=t, edges=e, edge_masks=k, vocabulary=VOC, binary=binary,
                triangle_codes=c, land_cover_codes="x 0 = none"))  # fmt: skip
    out.append("ply edges | unnamed bit | " + outcome(
        write_ply, v, edges=e, edge_properties=np.array([1, 1 << 30]), vocabulary=VOC))  # fmt: skip
    out.append("vtk | unnamed bit | " + outcome(
        write_vtk, v, triangles=t, edges=e, edge_masks=np.array([1, 1 << 30]),
        vocabulary=VOC))  # fmt: skip
    out.append("vtk field | reserved name | " + outcome(
        write_vtk, v, triangles=t, edges=e, edge_masks=k, vocabulary=VOC,
        fields=[("elevation", "x")]))  # fmt: skip
    wkt = ('ENGCRS["local",EDATUM["d"],CS[Cartesian,2],AXIS["x",east],AXIS["y",north],'
           'LENGTHUNIT["metre",1]]')  # fmt: skip
    for name, crs in (("EPSG:25833", "EPSG:25833"), ("OGC:CRS84", "OGC:CRS84"),
                      ("PROJ string", "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"),
                      ("WKT ascii", wkt),
                      ("WKT non-ascii", wkt.replace("local", "l\u00f8kal"))):  # fmt: skip
        try:
            out.append(f"crs_label | {name} | {crs_label(crs)}")
        except Exception as exc:  # a crash is what the probe reports
            out.append(f"crs_label | {name} | {type(exc).__name__}: {exc}")
    return out


def tiff(z: np.ndarray) -> bytes:
    """A little-endian GeoTIFF, float32, one strip, EPSG:25833, pixel is point."""
    rows, cols = z.shape
    keys = [1, 1, 0, 3, 1024, 0, 1, 1, 1025, 0, 1, 2, 3072, 0, 1, 25833]
    tags: list[tuple[int, int, list[float] | list[int]]] = [
        (256, 4, [cols]), (257, 4, [rows]), (258, 3, [32]), (259, 3, [1]), (262, 3, [1]),
        (273, 4, [0]), (277, 3, [1]), (278, 4, [rows]), (279, 4, [z.nbytes]), (284, 3, [1]),
        (339, 3, [3]), (33550, 12, [10.0, 5.0, 0.0]), (33922, 12, [0, 0, 0, X, Y, 0]),
        (34735, 3, keys),
    ]  # fmt: skip
    fmt = {3: "H", 4: "I", 12: "d"}
    extra_at = 8 + 2 + 12 * len(tags) + 4
    extra = b""
    entries = b""
    blobs = []
    for tag, kind, values in tags:
        blob = struct.pack(f"<{len(values)}{fmt[kind]}", *values)
        blobs.append((tag, kind, len(values), blob))
    image_at = extra_at + sum(len(b) for *_, b in blobs if len(b) > 4)
    for tag, kind, count, blob in blobs:
        if tag == 273:
            blob = struct.pack("<I", image_at)
        if len(blob) <= 4:
            entries += struct.pack("<HHI", tag, kind, count) + blob.ljust(4, b"\0")
        else:
            entries += struct.pack("<HHII", tag, kind, count, extra_at + len(extra))
            extra += blob
    ifd = struct.pack("<H", len(tags)) + entries + b"\0\0\0\0"
    return b"II*\0" + struct.pack("<I", 8) + ifd + extra + z.astype("<f4").tobytes()


def part2(copy: Path, tmp: Path) -> list[str]:
    z = np.random.default_rng(16).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)
    dem = tmp / "dém.tif"
    dem.write_bytes(tiff(z))
    crs = {"type": "name", "properties": {"name": "EPSG:25833"}}
    square = [[X + 12.3, Y - 73.3], [X + 187.7, Y - 72.9], [X + 186.1, Y - 6.7],
              [X + 13.9, Y - 7.1], [X + 12.3, Y - 73.3]]  # fmt: skip
    (tmp / "domain.geojson").write_text(json.dumps(
        {"type": "Polygon", "crs": crs, "coordinates": [square]}))  # fmt: skip

    def rect(x0: float, y0: float, x1: float, y1: float) -> list[list[list[float]]]:
        return [[[X + x0, Y + y0], [X + x1, Y + y0], [X + x1, Y + y1], [X + x0, Y + y1],
                 [X + x0, Y + y0]]]  # fmt: skip

    features = [
        {"type": "Feature", "properties": {"Code_18": code}, "geometry": {"type": "Polygon",
         "coordinates": rect(*box)}}
        for code, box in (("312", (-50.3, -150.1, 100.2, 50.4)),
                          ("512", (100.2, -150.1, 300.2, 50.4)))
    ]  # fmt: skip
    cover = tmp / "lånd.geojson"
    cover.write_text(json.dumps({"type": "FeatureCollection", "crs": crs, "features": features}))
    dem_args = ["--dem", str(dem), "--domain", str(tmp / "domain.geojson"), "--tolerance", "1",
                "--features", str(cover), "--features-map", "corine"]  # fmt: skip
    cases: list[tuple[str, list[str], list[str]]] = []
    for mode in ("--ascii", "--binary"):
        for crs_name, crs_text in (("plain", "EPSG:25833"), ("non-ascii", "ETRS89 °N"),
                                   ("cr", "x\ry"), ("empty", "")):  # fmt: skip
            cases.append((f"mesh river --flat {mode} --crs {crs_name} | .vtk",
                          ["mesh", "river", "--flat", mode, "--crs", crs_text, "--out", "a.vtk"],
                          ["a.vtk"]))  # fmt: skip
            cases.append((f"mesh river --flat {mode} --crs {crs_name} | .ply --out-edges",
                          ["mesh", "river", "--flat", mode, "--crs", crs_text, "--out", "a.ply",
                           "--out-edges", "e.ply"], ["a.ply", "e.ply"]))  # fmt: skip
        cases.append((f"mesh --dem {mode} | .vtk", ["mesh", *dem_args, mode, "--out", "d.vtk"],
                      ["d.vtk"]))  # fmt: skip
        cases.append((f"mesh --dem {mode} | .ply --out-edges",
                      ["mesh", *dem_args, mode, "--out", "d.ply", "--out-edges", "f.ply"],
                      ["d.ply", "f.ply"]))  # fmt: skip
        cases.append((f"mesh --dem {mode} | .ply --record",
                      ["mesh", *dem_args, mode, "--out", "g.ply", "--record", "r.json"],
                      ["g.ply"]))  # fmt: skip
    cases.append(("mesh river --flat | .obj", ["mesh", "river", "--flat", "--out", "a.obj"], []))
    cases.append(("mesh river --flat | .vtk --out-edges",
                  ["mesh", "river", "--flat", "--out", "a.vtk", "--out-edges", "e.ply"], []))
    out = []
    for name, args, files in cases:
        for f in files:
            (tmp / f).unlink(missing_ok=True)
        run = subprocess.run([sys.executable, "-c", CHILD, str(copy), *args], cwd=tmp,
                             capture_output=True, check=False)  # fmt: skip
        err = run.stderr.decode(errors="replace").replace(str(tmp), "<dir>")
        if "Traceback" in err:
            out.append(f"{name} | UNCAUGHT {err.strip().splitlines()[-1]}")
        elif run.returncode != 0:
            last = [ln.strip("│ ") for ln in err.splitlines() if ln.strip("│╰╯─ ")]
            out.append(f"{name} | refused (exit {run.returncode}): {last[-1]}")
        else:
            data = [(tmp / f).read_bytes().replace(str(tmp).encode(), b"<dir>") for f in files]
            marks = [m for m in MARKS if any(m in b for b in data)]
            out.append(f"{name} | " + " ".join(hashlib.sha256(b).hexdigest()[:16] for b in data)
                       + (f" | holds {b', '.join(marks).decode()}" if marks else ""))  # fmt: skip
            if "--record" in args:  # the escaped names; the rest of it holds timings
                record = json.loads((tmp / "r.json").read_text(encoding="ascii"))
                out.append(f"{name} | record " + json.dumps({k: record.get(k) for k in RECORD}))
    return out


def main() -> None:
    copy = Path(sys.argv[1]).resolve()
    import tin_engine

    assert Path(tin_engine.__file__).is_relative_to(copy / "src_python"), tin_engine.__file__
    os.environ["COLUMNS"] = "200"  # one refusal per line in the command's error box
    with tempfile.TemporaryDirectory() as name:
        lines = part1() + part2(copy, Path(name).resolve())
    print("\n".join(lines))


if __name__ == "__main__":
    main()
