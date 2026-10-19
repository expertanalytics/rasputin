"""23c-4's short runs at a cell size binary cannot hold (@perf, 2026-10-19).

python cellsize_runs.py PKG TESTS RECORD.jsonl SCRATCH

PKG: a tree's Release package (`build-bench/pkg`, built by bench.py), imported
ahead of any installed tin_engine. TESTS: the 23c-4 tree's `tests/python`, whose
`pieces_fixtures` holds the oracles (`problems`, with `near=True` judging seams
within 1e-6 m of a lattice line at a decimal spacing); the same oracles for
both trees. Writes one JSON line per run: the scene, the spacing, `--pieces`,
the exit code, the wall time of `rasputin mesh` in-process, the oracles' time
and findings, and a sha256 over every file in `<out>.pieces`.

Scenes: `square`, the design's (`hills()` at 41 x 41 nodes, `SQUARE`,
`--tolerance 1`); `big`, the same hills at 201 x 201 nodes and `SQUARE` scaled
by 5 about the tie point (whole half-decimetres, so on the 1 mm snap grid).
"""

import hashlib
import json
import sys
import time
from pathlib import Path

pkg, tests, record, scratch = sys.argv[1:5]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path[:0] = [pkg, tests]

import pieces_fixtures as pf  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402

import tin_engine  # noqa: E402

BIG_RING = [(pf.TIE_X + 5 * (x - pf.TIE_X), pf.TIE_Y + 5 * (y - pf.TIE_Y)) for x, y in pf.SQUARE]
SCENES = {"square": (41, pf.SQUARE, ("4", "16")), "big": (201, BIG_RING, ("4", "16", "64"))}
SPACINGS = (9.7, 9.9, 10.0)


def digest(d: Path | None) -> str | None:
    if d is None or not (d / "index.json").exists():
        return None
    h = hashlib.sha256()
    for f in sorted(d.iterdir()):
        h.update(f.name.encode() + b"\0" + f.read_bytes())
    return h.hexdigest()


with open(record, "w") as out:
    out.write(json.dumps({"tin_engine": tin_engine.__file__}) + "\n")
    for scene, (n, ring, counts) in SCENES.items():
        array = pf.hills(n)
        for h in SPACINGS:
            for pieces in counts:
                where = Path(scratch) / f"{scene}-{h}-{pieces}"
                dem = pf.write_dem(where / "dem.tif", array, h=h)
                dom = pf.write_domain(where / "domain.geojson", ring)
                args = ("--dem", str(dem), "--domain", str(dom), "--tolerance", "1",
                        "--pieces", pieces)  # fmt: skip
                t0 = time.perf_counter()
                done = pf.run(where / "run", *args)
                wall = time.perf_counter() - t0
                row = {"scene": scene, "h": h, "pieces": pieces, "exit": done.code,
                       "wall_s": round(wall, 3), "pieces_sha256": digest(done.dir)}  # fmt: skip
                if done.code != 0:
                    row["message"] = " ".join(done.output.split())[-300:]
                else:
                    pf.read_index_and_pieces(done)
                    t1 = time.perf_counter()
                    found = pf.problems(done, array, Polygon(ring), 1.0, h, near=h != 10.0)
                    row |= {"files": len(done.pieces), "problems": found,
                            "check_s": round(time.perf_counter() - t1, 3),
                            "triangles": sum(len(p.triangles) for p in done.pieces)}  # fmt: skip
                out.write(json.dumps(row) + "\n")
                out.flush()
                print(scene, h, pieces, row["exit"], row["wall_s"], len(row.get("problems", [])))
