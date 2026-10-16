"""Increment 30c's safety net: what `assemble` and `decode_window` return, as hashes.

`docs/increments/30c-dem-read-speed.md`, section 6 ("The byte-identical
gate"). Two modes, each printing one line per result; run a mode before and
after the change, with the interpreter whose `tin_engine` is the tree to
measure (the first line prints that path, and the run stops if it is not an
installed package). Every line of the gate's base file that starts with
`fixture ` or `mesh ` must appear unchanged in the run
(`grep -E '^(fixture|mesh) '`); the only other such lines allowed are those
of tests the red suite added, which section 6 lists. The rest is pytest's and
the CLI's own output.

    fixtures   run the suites that reach the DEM code in-process, recording
               every call of `assemble` (from `tin_engine.mosaic`, from
               `dem_input`, `catchment` or a test) and every call of
               `decode_window` (from `io.cog` or `io.repository`), keyed by
               test id and call number.
    mesh C...  run `rasputin mesh` on each named catchment as @perf's profile
               did (`docs/benchmarks/2026-10-06/bottlenecks/README.md`,
               "Method"), and hash the assembled DEM, its seam report and the
               `.vtk` written.

An `assemble` line: the canvas's shape and dtype, `sha`, a hash of its meta
and its bytes (NaN payloads included), and `seams`, every seam's two names,
node count, largest and median, with the floats as `repr`. A
`decode_window` line: the window, dtype, and a hash of the tile's meta and
bytes. A call that raises leaves no line.

From the repository root (`tests/python` on the path for `fixtures`, and
pytest importable; the suites' `addopts` are cleared, so pytest-cov is not
needed):

    PYTHONPATH=tests/python python docs/increments/30c-probes/dem_bytes.py fixtures
    python docs/increments/30c-probes/dem_bytes.py mesh numedalslagen skiensvassdraget

The data paths are `../rasputin_data` and `../rasputin_scratch/norway`,
next to the main checkout (`RASPUTIN_DATA`, `RASPUTIN_SCRATCH` override them).
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

import tin_engine
import tin_engine.catchment as catchment
import tin_engine.cli as cli
import tin_engine.dem_input as dem_input
import tin_engine.io.cog as cog
import tin_engine.io.repository as repository
import tin_engine.mosaic as mosaic

HERE = Path(__file__).resolve()
GIT = ["git", "-C", str(HERE.parent), "rev-parse", "--path-format=absolute", "--git-common-dir"]
MAIN = Path(subprocess.run(GIT, capture_output=True, text=True, check=True).stdout.strip()).parent
DATA = Path(os.environ.get("RASPUTIN_DATA", MAIN.parent / "rasputin_data"))
SCRATCH = Path(os.environ.get("RASPUTIN_SCRATCH", MAIN.parent / "rasputin_scratch" / "norway"))
SUITES = [
    "test_mosaic.py",
    "test_mosaic_windowed.py",
    "test_dem_input.py",
    "test_dem_input_domain.py",
    "test_io_cog.py",
    "test_io_repository.py",
    "test_io_cache.py",
    "test_cli_mesh_dem.py",
    "test_cli_mesh_mosaic.py",
    "test_cli_mesh_domain.py",
    "test_cli_mesh_domain_crs.py",
    "test_cli_mesh_cache.py",
    "test_cli_mesh_geographic.py",
    "test_cli_mesh_plain_output.py",
]
#: Where each function is looked up at call time: (module, attribute).
ASSEMBLE = [(mosaic, "assemble"), (dem_input, "assemble"), (catchment, "assemble")]
DECODE = [(cog, "decode_window"), (repository, "decode_window")]


def tile_sha(tile: Any) -> str:
    array = np.ascontiguousarray(np.asarray(tile.array))
    h = hashlib.sha256(repr(tile.meta).encode() + b"|" + str(array.dtype).encode() + b"|")
    h.update(array.tobytes())
    return h.hexdigest()[:16]


def assemble_line(out: Any) -> str:
    array = np.asarray(out.tile.array)
    seams = ";".join(
        f"{s.first}/{s.second}/{s.nodes}/{s.largest!r}/{s.median!r}" for s in out.seams
    )
    shape = "x".join(map(str, array.shape))
    return f"shape={shape} dtype={array.dtype} sha={tile_sha(out.tile)} seams=[{seams}]"


def decode_line(window: Any, out: Any) -> str:
    w = f"{window.row0},{window.col0},{window.rows},{window.cols}"
    return f"window={w} dtype={np.asarray(out.array).dtype} sha={tile_sha(out)}"


def patched(record: Any) -> list[tuple[Any, str, Any]]:
    """Wrap `assemble` and `decode_window` wherever they are looked up; return
    what to restore. `record(kind, text)` gets one line per call that returns."""
    original_assemble, original_decode = mosaic.assemble, cog.decode_window
    inside = [0]

    def assemble(*args: Any, **kwargs: Any) -> Any:
        inside[0] += 1
        try:
            out = original_assemble(*args, **kwargs)
        finally:
            inside[0] -= 1
        record("assemble", assemble_line(out))
        return out

    def decode_window(source: Any, meta: Any, dtype: Any, window: Any, **kwargs: Any) -> Any:
        out = original_decode(source, meta, dtype, window, **kwargs)
        record("decode_window" + (" in assemble" if inside[0] else ""), decode_line(window, out))
        return out

    saved = [(m, a, getattr(m, a)) for m, a in ASSEMBLE + DECODE]
    for m, a in ASSEMBLE:
        setattr(m, a, assemble)
    for m, a in DECODE:
        setattr(m, a, decode_window)
    return saved


def fixtures() -> int:
    import pytest

    calls: dict[str, int] = {}
    lines: list[str] = []

    def record(kind: str, text: str) -> None:
        test = os.environ.get("PYTEST_CURRENT_TEST", "?").rsplit(" (", 1)[0]
        calls[test] = calls.get(test, 0) + 1
        lines.append(f"fixture {test} #{calls[test]} {kind} {text}")

    saved = patched(record)
    paths = [str(Path.cwd() / "tests" / "python" / s) for s in SUITES]
    status = int(pytest.main(["-q", "-o", "addopts=", "-p", "no:cacheprovider", *paths]))
    for m, a, f in saved:
        setattr(m, a, f)
    print(*lines, sep="\n")
    print(f"fixtures: pytest exit {status}, {len(lines)} calls recorded")
    return status


def mesh(names: list[str]) -> None:
    for name in names:
        seen: list[str] = []
        saved = patched(
            lambda kind, text, seen=seen: seen.append(text) if kind == "assemble" else None
        )
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / f"{name}.vtk"
            cli.app(["mesh", "--dem", str(DATA / "DTM10_UTM33_20260925"),
                     "--domain", str(SCRATCH / name / f"{name}_outline_nve.geojson"),
                     "--features", str(DATA / "corine2018_dtm10_utm33.gpkg"),
                     "--features-layer", "corine2018", "--features-map", "corine",
                     "--tolerance", "10", "--binary", "--out", str(target)],
                    standalone_mode=False)  # fmt: skip
            vtk = hashlib.sha256(target.read_bytes()).hexdigest()[:16]
        for m, a, f in saved:
            setattr(m, a, f)
        assert len(seen) == 1, seen
        print(f"mesh {name} {seen[0]} vtk={vtk}", flush=True)


def main(argv: list[str]) -> int:
    where = Path(tin_engine.__file__).parent
    print(f"tin_engine from {where}", flush=True)
    if "site-packages" not in where.parts:
        raise SystemExit(f"tin_engine is not an installed package here ({where}); see section 6")
    mode, rest = argv[0], argv[1:]
    if mode == "fixtures":
        return fixtures()
    if mode == "mesh":
        mesh(rest)
        return 0
    raise SystemExit(f"unknown mode {mode!r}: fixtures or mesh")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
