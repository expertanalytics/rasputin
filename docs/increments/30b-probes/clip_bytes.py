"""Increment 30b's safety net: what `open_features` and `pre_clip` return, as hashes.

`docs/increments/30b-clip-speed.md`, section 6 ("The byte-identical gate").
Two modes, each printing one line per result; run a mode before and after the
change, with the interpreter whose `tin_engine` is the tree to measure (the
first line prints that path). The lines that start with `fixture ` or `mesh `
must be identical to the gate's base file, `base_5e2fbe0.txt`
(`grep -E '^(fixture|mesh) '`; `base_a7154ec.txt` is the earlier base, kept
as history); the rest is pytest's and the CLI's own output.

    fixtures   run the suites that reach the features code in-process,
               recording every call of `open_features` (from the module or
               from the CLI) and every direct call of `pre_clip` (one not made
               from inside `open_features`), keyed by test id and call number;
               then `open_features` on three hand-made sources (`cases`),
               keyed `probe::<case>`.
    mesh C...  run `rasputin mesh` on each named catchment as @perf's profile
               did (`docs/benchmarks/2026-10-06/bottlenecks/README.md`,
               "Method"), and hash the feature set (every feature's fid,
               mask, code, lines and polygon, in order, and the four counts)
               and the `.vtk` written.

A feature set's line: `features` kept, `outside`, `clipped`, `empty`,
`vertices` (in all kept lines), `lines` (their number) and `sha`, a hash of
every feature's fid, mask, code, each line's WKB and the polygon's WKB, in
order. A `pre_clip` line: the number of chains and a hash of their WKB.

From the repository root (`tests/python` on the path for `fixtures`, and
pytest importable; the suites' `addopts` are cleared, so pytest-cov is not
needed):

    PYTHONPATH=tests/python python docs/increments/30b-probes/clip_bytes.py fixtures
    python docs/increments/30b-probes/clip_bytes.py mesh numedalslagen skiensvassdraget

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

import shapely

import tin_engine
import tin_engine.cli as cli
import tin_engine.feature_input as fi

HERE = Path(__file__).resolve()
GIT = ["git", "-C", str(HERE.parent), "rev-parse", "--path-format=absolute", "--git-common-dir"]
MAIN = Path(subprocess.run(GIT, capture_output=True, text=True, check=True).stdout.strip()).parent
DATA = Path(os.environ.get("RASPUTIN_DATA", MAIN.parent / "rasputin_data"))
SCRATCH = Path(os.environ.get("RASPUTIN_SCRATCH", MAIN.parent / "rasputin_scratch" / "norway"))
SUITES = [
    "test_feature_input.py",
    "test_cli_mesh_features.py",
    "test_cli_mesh_multi_features.py",
    "test_cli_mesh_landcover.py",
    "test_cli_mesh_edge_strip.py",
    "test_cli_mesh_plain_output.py",
]


def wkb(geometry: Any) -> bytes:
    return b"" if geometry is None else bytes(shapely.to_wkb(geometry))


def set_line(found: Any) -> str:
    h = hashlib.sha256()
    vertices = lines = 0
    for f in found.features:
        h.update(f"{f.fid!r}|{f.mask}|{f.code}|{len(f.lines)}|".encode())
        for line in f.lines:
            h.update(wkb(line))
            vertices += len(line.coords)
        lines += len(f.lines)
        h.update(b"|" + wkb(f.polygon) + b"#")
    return (
        f"features={len(found.features)} outside={found.outside} clipped={found.clipped} "
        f"empty={found.empty} vertices={vertices} lines={lines} sha={h.hexdigest()[:16]}"
    )


def chains_line(chains: Any) -> str:
    h = hashlib.sha256(b"".join(wkb(c) + b"#" for c in chains))
    return f"chains={len(chains)} sha={h.hexdigest()[:16]}"


def fixtures() -> int:
    import pytest

    original_open, original_clip = fi.open_features, fi.pre_clip
    calls: dict[str, int] = {}
    lines: list[str] = []
    inside = [0]

    def key(kind: str) -> str:
        test = os.environ.get("PYTEST_CURRENT_TEST", "?").rsplit(" (", 1)[0]
        calls[test] = calls.get(test, 0) + 1
        return f"fixture {test} #{calls[test]} {kind}"

    def open_features(*args: Any, **kwargs: Any) -> Any:
        inside[0] += 1
        try:
            out = original_open(*args, **kwargs)
        finally:
            inside[0] -= 1
        lines.append(f"{key('open_features')} {set_line(out)}")
        return out

    def pre_clip(*args: Any, **kwargs: Any) -> Any:
        out = original_clip(*args, **kwargs)
        if not inside[0]:
            lines.append(f"{key('pre_clip')} {chains_line(out)}")
        return out

    fi.open_features, fi.pre_clip, cli.open_features = open_features, pre_clip, open_features
    paths = [str(Path.cwd() / "tests" / "python" / s) for s in SUITES]
    status = int(pytest.main(["-q", "-o", "addopts=", "-p", "no:cacheprovider", *paths]))
    fi.open_features, fi.pre_clip, cli.open_features = original_open, original_clip, original_open
    lines += cases()
    print(*lines, sep="\n")
    print(f"fixtures: pytest exit {status}, {len(lines)} calls recorded")
    return status


def cases() -> list[str]:
    """Hand-made sources that no suite at the base has (section 7, P5), each
    through `open_features` on the suites' 300 m box domain: an invalid
    polygon whose shell lies far outside the region and whose hole lies inside
    the domain (its shell's box misses the region's); the same with a shell
    whose box overlaps the region's; and an open line whose two ends lie
    outside the region and which crosses the domain."""
    from shapely.geometry import LineString, Polygon

    from feature_fixtures import Feat, at, domain_of, open_one, square, write_geojson

    box, hole = domain_of(square(0, 0, 300, 300)), [square(100, 100, 200, 200).exterior.coords]
    ell = [at(500, -200), at(600, -200), at(600, 600), at(-200, 600), at(-200, 500), at(500, 500)]
    sources = {
        "P5 hole inside, shell's box misses the region's": Polygon(
            square(2_000, 2_000, 2_100, 2_100).exterior.coords, hole
        ),
        "P5 hole inside, shell's box overlaps the region's": Polygon(ell, hole),
        "P5 line with both ends outside the region": LineString([at(-500, 150), at(800, 150)]),
    }
    out: list[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        for name, geometry in sources.items():
            feature = Feat("a", geometry, {"property": "road"})
            path = write_geojson(Path(tmp) / "f.geojson", [feature])
            out.append(f"fixture probe::{name} #1 open_features {set_line(open_one(path, box))}")
    return out


def mesh(names: list[str]) -> None:
    original = cli.open_features
    for name in names:
        seen: list[str] = []

        def grab(*args: Any, seen: list[str] = seen, **kw: Any) -> Any:
            out = original(*args, **kw)
            seen.append(set_line(out))
            return out

        cli.open_features = grab
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / f"{name}.vtk"
            cli.app(["mesh", "--dem", str(DATA / "DTM10_UTM33_20260925"),
                     "--domain", str(SCRATCH / name / f"{name}_outline_nve.geojson"),
                     "--features", str(DATA / "corine2018_dtm10_utm33.gpkg"),
                     "--features-layer", "corine2018", "--features-map", "corine",
                     "--tolerance", "10", "--binary", "--out", str(target)],
                    standalone_mode=False)  # fmt: skip
            vtk = hashlib.sha256(target.read_bytes()).hexdigest()[:16]
        assert len(seen) == 1, seen
        print(f"mesh {name} {seen[0]} vtk={vtk}", flush=True)
    cli.open_features = original


def main(argv: list[str]) -> int:
    print(f"tin_engine from {Path(tin_engine.__file__).parent}", flush=True)
    mode, rest = argv[0], argv[1:]
    if mode == "fixtures":
        return fixtures()
    if mode == "mesh":
        mesh(rest)
        return 0
    raise SystemExit(f"unknown mode {mode!r}: fixtures or mesh")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
