"""Increment 30a's safety net: what `label_triangles` and `regions` return, as hashes.

`docs/increments/30a-landcover-speed.md`, section 6 ("The byte-identical
gate"). Three modes, each printing one line per result; run a mode before
and after the change, with the interpreter whose `tin_engine` is the tree to
measure (the first line prints that path). The lines that start with
`fixture ` or `mesh ` must be identical to `base_a7154ec.txt`'s
(`grep -E '^(fixture|mesh) '`); the rest is pytest's and the CLI's own output.

    fixtures   run the suites that call the land-cover code
               (`test_landcover.py`, `test_cli_mesh_landcover.py`,
               `test_cli_mesh_multi_features.py`) in-process, recording every
               call of `label_triangles` (from the module or from the CLI) and
               every direct call of `regions`, keyed by test id and call number.
    mesh C...  run `rasputin mesh` on each named catchment as @perf's profile
               did (`docs/benchmarks/2026-10-06/bottlenecks/README.md`,
               "Method"), and hash the codes, the four counts of the stderr
               line, and the `.vtk` written. `--save DIR` also keeps
               `label_triangles`' inputs in DIR (`<C>.npz`), for `replay`.
    replay DIR C...  call `label_triangles` on inputs `mesh --save` kept, and
               print the same codes-and-counts line as `mesh` (no `.vtk`
               hash), with the call's own seconds on stderr.

From the repository root (`tests/python` on the path for `fixtures`, and
pytest importable; the suites' `addopts` are cleared, so pytest-cov is not
needed):

    PYTHONPATH=tests/python python docs/increments/30a-probes/landcover_bytes.py fixtures
    python docs/increments/30a-probes/landcover_bytes.py mesh numedalslagen skiensvassdraget

The data paths are `../rasputin_data` and `../rasputin_scratch/norway`,
next to the main checkout (`RASPUTIN_DATA`, `RASPUTIN_SCRATCH` override them).
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import numpy as np
import shapely

import tin_engine
import tin_engine.cli as cli
import tin_engine.landcover as lc

HERE = Path(__file__).resolve()
GIT = ["git", "-C", str(HERE.parent), "rev-parse", "--path-format=absolute", "--git-common-dir"]
MAIN = Path(subprocess.run(GIT, capture_output=True, text=True, check=True).stdout.strip()).parent
DATA = Path(os.environ.get("RASPUTIN_DATA", MAIN.parent / "rasputin_data"))
SCRATCH = Path(os.environ.get("RASPUTIN_SCRATCH", MAIN.parent / "rasputin_scratch" / "norway"))


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()[:16]


def labels_line(labels: Any) -> str:
    codes = np.ascontiguousarray(labels.codes)
    return (
        f"triangles={len(codes)} dtype={codes.dtype} regions={labels.regions} "
        f"outside={labels.outside} overlapped={labels.overlapped} thin={labels.thin} "
        f"codes={sha(codes.tobytes())}"
    )


def fixtures() -> int:
    import pytest

    original_label, original_regions = lc.label_triangles, lc.regions
    calls: dict[str, int] = {}
    lines: list[str] = []
    inside = [0]

    def key(kind: str) -> str:
        test = os.environ.get("PYTEST_CURRENT_TEST", "?").rsplit(" (", 1)[0]
        calls[test] = calls.get(test, 0) + 1
        return f"fixture {test} #{calls[test]} {kind}"

    def label(*args: Any, **kwargs: Any) -> Any:
        inside[0] += 1
        try:
            out = original_label(*args, **kwargs)
        finally:
            inside[0] -= 1
        lines.append(f"{key('label_triangles')} {labels_line(out)}")
        return out

    def regions(*args: Any, **kwargs: Any) -> Any:
        out = original_regions(*args, **kwargs)
        if not inside[0]:
            ids = np.ascontiguousarray(out, dtype=np.int64)
            lines.append(f"{key('regions')} n={len(ids)} ids={sha(ids.tobytes())}")
        return out

    lc.label_triangles, lc.regions, cli.label_triangles = label, regions, label
    suites = ["test_landcover.py", "test_cli_mesh_landcover.py", "test_cli_mesh_multi_features.py"]
    paths = [str(Path.cwd() / "tests" / "python" / s) for s in suites]
    status = int(pytest.main(["-q", "-o", "addopts=", "-p", "no:cacheprovider", *paths]))
    print(*lines, sep="\n")
    print(f"fixtures: pytest exit {status}, {len(lines)} calls recorded")
    return status


def save_inputs(path: Path, vertices: Any, triangles: Any, edges: Any, **kw: Any) -> None:
    """The call's arrays, and the polygons as WKB joined end to end with offsets."""
    wkb = [shapely.to_wkb(p) for p, _ in kw["polygons"]]
    np.savez(path, vertices=np.asarray(vertices), triangles=np.asarray(triangles),
             edges=np.asarray(edges), margin=kw["margin"],
             codes=np.array([c for _, c in kw["polygons"]], dtype=np.int64),
             wkb=np.frombuffer(b"".join(wkb), dtype=np.uint8),
             ends=np.cumsum([len(w) for w in wkb], dtype=np.int64))  # fmt: skip


def load_inputs(path: Path) -> tuple[Any, Any, Any, list[tuple[Any, int]], float]:
    z = np.load(path)
    blob, ends = z["wkb"].tobytes(), z["ends"].tolist()
    starts = [0, *ends[:-1]]
    polygons = [(shapely.from_wkb(blob[i:j]), int(c))
                for i, j, c in zip(starts, ends, z["codes"], strict=True)]  # fmt: skip
    return z["vertices"], z["triangles"], z["edges"], polygons, float(z["margin"])


def mesh(names: list[str], save: Path | None) -> None:
    original = cli.label_triangles
    for name in names:
        seen: list[str] = []

        def grab(*args: Any, name: str = name, seen: list[str] = seen, **kw: Any) -> Any:
            if save is not None:
                save_inputs(save / f"{name}.npz", *args, **kw)
            out = original(*args, **kw)
            seen.append(labels_line(out))
            return out

        cli.label_triangles = grab
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / f"{name}.vtk"
            cli.app(["mesh", "--dem", str(DATA / "DTM10_UTM33_20260925"),
                     "--domain", str(SCRATCH / name / f"{name}_outline_nve.geojson"),
                     "--features", str(DATA / "corine2018_dtm10_utm33.gpkg"),
                     "--features-layer", "corine2018", "--features-map", "corine",
                     "--tolerance", "10", "--binary", "--out", str(target)],
                    standalone_mode=False)  # fmt: skip
            vtk = sha(target.read_bytes())
        assert len(seen) == 1, seen
        print(f"mesh {name} {seen[0]} vtk={vtk}", flush=True)
    cli.label_triangles = original


def replay(folder: Path, names: list[str]) -> None:
    for name in names:
        vertices, triangles, edges, polygons, margin = load_inputs(folder / f"{name}.npz")
        start = time.perf_counter()
        out = lc.label_triangles(vertices, triangles, edges, polygons=polygons, margin=margin)
        print(f"replay {name}: {time.perf_counter() - start:.3f} s", file=sys.stderr)
        print(f"mesh {name} {labels_line(out)}", flush=True)


def main(argv: list[str]) -> int:
    print(f"tin_engine from {Path(tin_engine.__file__).parent}", flush=True)
    mode, rest = argv[0], argv[1:]
    if mode == "fixtures":
        return fixtures()
    if mode == "mesh":
        save = None
        if rest[:1] == ["--save"]:
            save, rest = Path(rest[1]), rest[2:]
            save.mkdir(parents=True, exist_ok=True)
        mesh(rest, save)
        return 0
    if mode == "replay":
        replay(Path(rest[0]), rest[1:])
        return 0
    raise SystemExit(f"unknown mode {mode!r}: fixtures, mesh or replay")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
