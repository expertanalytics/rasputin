#!/usr/bin/env python3
"""Increment 137 section 8: wall time and peak resident memory of three
commands, with `--dem-heights` off on master and head and on at head.

    python3 run_137.py --repo <head worktree> --scratch <dir> --repeats 3

`<scratch>/pkg-base` and `<scratch>/pkg-head` are symlink trees of each
tree's `src_python/tin_engine` plus its Release `_core` (the README says how
they are made). Each sample is one `bench.py _child --pkg` process under
`/usr/bin/time -l`; repeats are interleaved over (variant, case). Writes
`samples.jsonl` and `meshes.json` next to this file; meshes go to
`<scratch>/out` and stay out of the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def cases(repo: Path) -> dict[str, list[str]]:
    fixture = str(repo / "tests/fixtures/dem_archive/7908_3_10m_z33.tif")
    velhas = repo / "tests/fixtures/velhas"
    quarter = str(repo / "docs/benchmarks/2026-09-26/quarter.geojson")
    return {
        "tile": ["--dem", fixture, "--tolerance", "10"],
        "cut4": ["--dem", fixture, "--domain", quarter, "--tolerance", "10", "--pieces", "4"],
        "velhas": ["--dem", str(velhas / "anadem_velhas.tif"),
                   "--domain", str(velhas / "catchment.geojson"),
                   "--out-crs", "EPSG:31983", "--tolerance", "5"],
    }  # fmt: skip


VARIANTS = (("base-off", "pkg-base", False), ("head-off", "pkg-head", False),
            ("head-on", "pkg-head", True))  # fmt: skip


def one(repo: Path, scratch: Path, variant: str, pkg: str, on: bool, case: str,
        args: list[str], r: int) -> dict[str, object]:  # fmt: skip
    out = scratch / "out" / f"{case}-{variant}.vtk"
    argv = ["/usr/bin/time", "-l", sys.executable, str(repo / "tools/bench.py"), "_child",
            "--pkg", str(scratch / pkg), "--threads", "0", "--", "mesh", *args,
            *(["--dem-heights"] if on else []), "--out", str(out)]  # fmt: skip
    done = subprocess.run(argv, capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise SystemExit(f"{variant} {case} r={r} failed:\n{done.stderr[-2000:]}")
    bench = json.loads(next(x[6:] for x in done.stderr.splitlines() if x.startswith("BENCH ")))
    real = float(re.search(r"([\d.]+) real", done.stderr).group(1))  # type: ignore[union-attr]
    rss = int(re.search(r"(\d+)\s+maximum resident set size", done.stderr).group(1))  # type: ignore[union-attr]
    return {"variant": variant, "case": case, "repeat": r, "wall_s": real,
            "app_s": bench["app_s"], "peak_rss_mb": rss / 2**20, "phases": bench.get("phases", {}),
            "hardening": bench.get("hardening")}  # fmt: skip


def digest(path: Path) -> str:
    """A single file's bytes, or every piece file's name and bytes (the
    names inside `<out>.pieces/` are the same for every variant)."""
    if path.exists():
        return hashlib.sha256(path.read_bytes()).hexdigest()[:16]
    h = hashlib.sha256()
    for f in sorted(path.parent.glob(path.name + ".pieces/*")):
        h.update(f.name.encode() + f.read_bytes())
    return h.hexdigest()[:16]


def shape(out: Path) -> dict[str, object]:
    """Per case and file: triangles, whether the flag left vertices and
    triangles as they are without it, and the on file's mean class count."""
    import numpy as np

    from tin_engine.io.vtk_read import read_mesh

    found: dict[str, object] = {}
    for case in ("tile", "cut4", "velhas"):
        for off in sorted(out.glob(f"{case}-head-off.vtk*")):
            pairs = [(off, Path(str(off).replace("head-off", "head-on")))] if off.is_file() else [
                (f, out / f"{case}-head-on.vtk.pieces" / f.name) for f in sorted(off.glob("*.vtk"))]
            for f_off, f_on in pairs:
                m0, m1 = read_mesh(f_off.read_bytes()), read_mesh(f_on.read_bytes())
                t = m1.dem_heights
                assert t is not None
                found[f"{case}/{f_off.name}"] = {
                    "triangles": len(m1.triangles),
                    "same_mesh": bool(np.array_equal(m0.vertices, m1.vertices)
                                      and np.array_equal(m0.triangles, m1.triangles)),
                    "mean_classes": float(np.diff(t.offsets.astype(np.int64)).mean()),
                    "bytes_off": f_off.stat().st_size, "bytes_on": f_on.stat().st_size}  # fmt: skip
    return found


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--scratch", type=Path, required=True)
    p.add_argument("--repeats", type=int, default=3)
    a = p.parse_args()
    (a.scratch / "out").mkdir(parents=True, exist_ok=True)
    samples = []
    with (HERE / "samples.jsonl").open("w") as fh:
        for r in range(a.repeats):
            for case, args in cases(a.repo).items():
                for variant, pkg, on in VARIANTS:
                    s = one(a.repo, a.scratch, variant, pkg, on, case, args, r)
                    samples.append(s)
                    fh.write(json.dumps(s) + "\n")
                    fh.flush()
                    print(f"{case:7} {variant:9} r={r} wall {s['wall_s']:.2f} s "
                          f"rss {s['peak_rss_mb']:.0f} MB", flush=True)  # fmt: skip
    meshes = {f"{c}-{v}": digest(a.scratch / "out" / f"{c}-{v}.vtk")
              for c in cases(a.repo) for v, _, _ in VARIANTS}  # fmt: skip
    meshes.update(shape(a.scratch / "out"))
    (HERE / "meshes.json").write_text(json.dumps(meshes, indent=1) + "\n")
    for case in cases(a.repo):
        for variant, _, _ in VARIANTS:
            g = [s for s in samples if s["case"] == case and s["variant"] == variant]
            print(case, variant, "median wall", statistics.median(float(s["wall_s"]) for s in g),
                  "median rss", statistics.median(float(s["peak_rss_mb"]) for s in g))  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
