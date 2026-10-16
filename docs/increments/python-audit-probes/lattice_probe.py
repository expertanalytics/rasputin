"""What every function PR A (`audit-lattice`) edits returns, on every input
shape it is built for and a few it is not (`python-audit-pr-a.md`, "The
differential probe").

Run against one revision's `tin_engine`, imported from a scratch copy:

    python3 tools/scratch_copy.py <rev> <dir>
    PYTHONPATH=<that value> .venv/bin/python \\
        docs/increments/python-audit-probes/lattice_probe.py <dir>

scratch_copy.py prints one line on stdout, the pytest command for the copy
(`cd <dir> && PYTHONPATH=<value> <python> -c ... tests/python/`); take
<value> from that line only (stderr may carry a `scratch_copy:` warning).
scratch_copy.py copies `_core*.so` from the `.venv` of the worktree it is
run from. When that worktree has none (it warns `no built _core`), copy one
by hand into `<dir>/src_python/tin_engine/`, from a venv whose C++ is
unchanged against <rev>: `git diff --quiet <rev> <that venv's head> --
include src bindings CMakeLists.txt` exits 0 (the paths scratch_copy.py
itself compares).
The copy's sitecustomize drops the editable finder, and the probe asserts
that `tin_engine` was imported from `<dir>/src_python`, so a run cannot
measure the work tree by mistake. It reads the DEM fixtures under
`<dir>/tests/fixtures/dtm10/` and writes no file.

Imports: the standard library, numpy, shapely and `tin_engine` only
(`tools/check_prohibited_deps.py` does not scan `docs/`).

One line per (entry, case): `entry | case | outcome`. The outcome is the
value, written exactly (a float as its type and `float.hex`, an array as
dtype, shape, its NaN, +inf and -inf counts and a hash of its bytes, a
geometry as its type and a hash of its WKB, a model or dataclass field by
field, wall-clock fields left out, a lazy field read out), or
`raised <class>: <message>`, the
message with any memory address blanked. Entries are the functions at both
revisions, public or private, that the PR edits or that call what it edits;
`diff` of two revisions' outputs is every changed result.

Cases. Grids: a 10 m UTM grid, a non-dyadic one (origin 0.3, spacing 0.1,
where another spelling of a node lands an ulp off it), one row, one
column, one node, a geographic grid (spacing 1/1200 degree), a negative
spacing, and `None` and a dict where a `RasterMeta` belongs. Values: clean,
one NaN, one +inf, one -inf, +inf beside -inf, one sentinel, in float32 and
float64, with and without a sentinel. DEM fixtures: the two-tile seam set
and the two-lattice set, whole and windowed, with and without a box;
`open_dem` on the seam set with a box, two domains and a domain reprojected
to UTM 32. Catchments: a cone (points, a lake, a reach) and a valley whose
catchment grows its window.
"""


from __future__ import annotations

import asyncio
import dataclasses
import hashlib
import math
import re
import sys
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import shapely

ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ROOT / "src_python"))

import tin_engine  # noqa: E402

assert Path(tin_engine.__file__).resolve().is_relative_to(ROOT / "src_python"), tin_engine.__file__

from tin_engine import (  # noqa: E402
    burn,
    catchment,
    catchment_batch,
    cli,
    dem_input,
    domain,
    grid_domain,
    mosaic,
    reference,
    target_grid,
)
from tin_engine.fetch import plan as fetch_plan  # noqa: E402
from tin_engine.gauge import Reach  # noqa: E402
from tin_engine.io.models import DemTile, IndexWindow, RasterMeta  # noqa: E402
from tin_engine.io.repository import TiffDemRepository, TileFootprint  # noqa: E402

Bounds = mosaic.Bounds  # its home at the base; PR A keeps the name importable there
ADDRESS = re.compile(r"0x[0-9a-f]{6,}")
UTM = "EPSG:25833"


def show(v: Any) -> str:
    """`v` written exactly, as the module docstring says."""
    if isinstance(v, BaseException):
        text = ADDRESS.sub("0x?", str(v)).replace("\n", " / ")
        return f"raised {type(v).__name__}: {text}"
    if isinstance(v, np.ndarray):
        counts = ""
        if v.dtype.kind == "f":
            nan, pos, neg = (int(f(v).sum()) for f in (np.isnan, np.isposinf, np.isneginf))
            counts = f" nan={nan} +inf={pos} -inf={neg}"
        digest = hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest()[:16]
        return f"array {v.dtype} {v.shape}{counts} sha={digest}"
    if isinstance(v, bool | type(None) | str | int) and not isinstance(v, np.integer):
        return repr(v)
    if isinstance(v, float | np.floating):
        return f"{type(v).__name__}:{float(v).hex()}"
    if isinstance(v, np.integer):
        return f"{type(v).__name__}:{int(v)}"
    if isinstance(v, shapely.Geometry):
        return f"{v.geom_type} sha={hashlib.sha256(shapely.to_wkb(v)).hexdigest()[:16]}"
    if isinstance(v, tuple | list):
        return "(" + ", ".join(show(x) for x in v) + ")"
    if isinstance(v, dict):
        return "{" + ", ".join(f"{k}={show(x)}" for k, x in v.items()) + "}"
    if hasattr(type(v), "model_fields"):  # a pydantic model
        fields = {k: getattr(v, k) for k in type(v).model_fields}
        return f"{type(v).__name__}{show(fields)}"
    if dataclasses.is_dataclass(v) and not isinstance(v, type):
        fields = {f.name: getattr(v, f.name) for f in dataclasses.fields(v)}
        fields = {k: list(x) if hasattr(x, "__next__") else x  # a lazy field, read out
                  for k, x in fields.items() if "seconds" not in k}  # fmt: skip
        return f"{type(v).__name__}{show(fields)}"
    return f"<{type(v).__name__}>"


def line(entry: str, case: str, run: Callable[..., Any], *args: Any) -> None:
    """One line: `run(*args)`'s outcome, a generator's items read out."""
    try:
        out = run(*args)
        if hasattr(out, "__next__"):
            out = list(out)
        print(f"{entry} | {case} | {show(out)}")
    except Exception as exc:  # the probe records every outcome, a refusal included
        print(f"{entry} | {case} | {show(exc)}")


def meta(**over: Any) -> RasterMeta:
    base: dict[str, Any] = {
        "x_min": 500_000.0, "y_max": 6_600_000.0, "delta_x": 10.0, "delta_y": 10.0, "cols": 5,
        "rows": 6, "epsg": 25833, "nodata": -32767.0, "nodata_source": "tag",
        "pixel_is_area": False, "vertical_unit_assumed": False,
    }  # fmt: skip
    return RasterMeta(**{**base, **over})


def grids() -> dict[str, Any]:
    return {
        "utm 10 m 6x5": meta(),
        "non-dyadic 7x6": meta(x_min=0.3, y_max=0.7, delta_x=0.1, delta_y=0.1, cols=6, rows=7,
                               nodata=None, nodata_source="absent"),
        "one row": meta(rows=1),
        "one column": meta(cols=1),
        "one node": meta(rows=1, cols=1),
        "geographic": meta(x_min=10.0, y_max=60.0, delta_x=1 / 1200, delta_y=1 / 1200, cols=8,
                           rows=9, epsg=4326, geographic=True),
        "negative spacing": meta(delta_x=-10.0),
        "None": None,
        "dict": {"x_min": 0.0},
    }  # fmt: skip


def real(m: Any) -> RasterMeta:
    """`m`, or the 10 m grid for a wrong type: the probe's own arithmetic
    reads this, so the entry itself is what meets the wrong type."""
    return m if isinstance(m, RasterMeta) else meta()


def a_box(m: Any, grow: float) -> Bounds:
    """The node rectangle of `m` grown by `grow` cells, written out by hand."""
    g = real(m)
    x1 = g.x_min + (g.cols - 1) * g.delta_x
    y0 = g.y_max - (g.rows - 1) * g.delta_y
    gx, gy = abs(g.delta_x) * grow, g.delta_y * grow
    return Bounds(x_min=min(g.x_min, x1) - gx - 1e-9, y_min=y0 - gy - 1e-9,
                  x_max=max(g.x_min, x1) + gx + 1e-9, y_max=g.y_max + gy + 1e-9)  # fmt: skip


def polygon(m: Any, grow: float) -> Any:
    b = a_box(m, grow)
    return shapely.box(b.x_min, b.y_min, b.x_max, b.y_max)


def extent(m: Any, grow: float) -> Any:
    crs = m.crs if isinstance(m, RasterMeta) else UTM
    return domain.check_extent(domain.DomainPolygon(polygon=polygon(m, grow), crs=crs), m)


def past(m: Any, grow: float) -> Any:
    return dem_input._past(a_box(m, grow), m)


def inside(m: Any, grow: float) -> Any:
    return reference.nodes_inside([polygon(m, grow).buffer(0)], m)


def lake_seed(m: Any, grow: float) -> Any:
    return catchment._seed_mask(m, polygon(m, grow), (0.0, 0.0))


def seed(m: Any, off: float) -> Any:
    g = real(m)
    return catchment._seed_mask(m, None, (g.x_min + off * g.delta_x, g.y_max - off * g.delta_y))


def off_node(m: Any) -> Any:
    g = real(m)
    c, r = np.arange(g.cols, dtype=np.float64), np.arange(g.rows, dtype=np.float64)
    on = np.column_stack([g.x_min + c * g.delta_x, np.full(c.size, g.y_max - r[-1] * g.delta_y)])
    off = on + np.array([[g.delta_x / 3, 0.0]])
    return cli._off_node(np.vstack([on, off, np.nextafter(on, np.inf)]), m)


def plan_object(m: Any) -> Any:
    g = real(m)
    blocks = -(-g.rows // 2) * -(-g.cols // 2)
    offsets = list(range(0, 8 * blocks, 8))
    page = SimpleNamespace(chunks=(2, 2), imagewidth=g.cols, dataoffsets=offsets,
                           databytecounts=[8] * blocks)  # fmt: skip
    return fetch_plan.plan_object("o", "u", page, m, a_box(m, -0.25), ())


def source_box(m: Any) -> Any:
    return fetch_plan.source_box(fetch_plan.FetchRequest(source="s", box=a_box(m, 0.0)), m)


def reprojected(m: Any) -> Any:
    return dem_input._reprojected(dem_input.DemRequest(sources=(Path("absent.tif"),)), [m, m])


def region(m: Any) -> Any:
    grid = target_grid.TargetGrid(crs=UTM, spacing=5, row0=-1_320_000, col0=100_000, rows=11,
                                  cols=9)  # fmt: skip
    grown = shapely.box(500_000.0, 6_599_950.0, 500_040.0, 6_600_000.0)
    return target_grid.source_region(grid, m, grown)


def grid_entries() -> None:
    for name, m in grids().items():
        for stride in (1, 2, 0):
            case = f"{name}, stride {stride}"
            line("grid_domain.subsample", case, grid_domain.subsample, m, stride)
        for grow, label in ((-0.25, "inside"), (0.5, "outside")):
            line("domain.check_extent", f"{name}, domain {label}", extent, m, grow)
            line("dem_input._past", f"{name}, box {label}", past, m, grow)
            line("reference.nodes_inside", f"{name}, polygon {label}", inside, m, grow)
            line("catchment._seed_mask", f"{name}, lake {label}", lake_seed, m, grow)
        for point, off in (("first node", 0), ("mid-cell", 0.5), ("outside", -3)):
            line("catchment._seed_mask", f"{name}, no lake, point {point}", seed, m, off)
        line("cli._off_node", name, off_node, m)
        line("fetch.plan.plan_object", name, plan_object, m)
        line("fetch.plan.source_box", name, source_box, m)
        line("dem_input._reprojected", name, reprojected, m)
        line("target_grid.source_region", name, region, m)


VALUES = {
    "clean": {},
    "nan": {(2, 2): math.nan},
    "+inf": {(2, 2): math.inf},
    "-inf": {(2, 2): -math.inf},
    "+inf beside -inf": {(2, 2): math.inf, (2, 3): -math.inf},
    "sentinel": {(2, 2): -32767.0},
}


def tile_of(values: dict[tuple[int, int], float], dtype: Any, nodata: float | None,
            x_min: float = 500_000.0) -> DemTile:  # fmt: skip
    """A 6 x 5 tile of a tilted plane, `values` written over it."""
    r, c = np.indices((6, 5))
    z = (100.0 + r * 3.0 + c * 2.0).astype(dtype)
    for (i, j), v in values.items():
        z[i, j] = v
    source = "absent" if nodata is None else "tag"
    return DemTile(meta=meta(x_min=x_min, nodata=nodata, nodata_source=source), array=z)


class Memory:
    """A `DemRepository` over tiles held in memory."""

    def __init__(self, tiles: dict[str, DemTile]) -> None:
        self.tiles = tiles

    def footprints(self) -> tuple[Any, ...]:
        return tuple(TileFootprint(name=n, meta=t.meta, dtype=t.array.dtype)
                     for n, t in sorted(self.tiles.items()))  # fmt: skip

    def load(self, name: str) -> DemTile:
        return self.tiles[name]

    def load_window(self, name: str, w: IndexWindow) -> DemTile:
        """The window, its meta moved by whole cells (the oracle's spelling)."""
        t = self.tiles[name]
        m = t.meta.model_copy(update={"x_min": t.meta.x_min + w.col0 * t.meta.delta_x,
                                      "y_max": t.meta.y_max - w.row0 * t.meta.delta_y,
                                      "rows": w.rows, "cols": w.cols})  # fmt: skip
        return DemTile(meta=m, array=t.array[w.row0 : w.row0 + w.rows, w.col0 : w.col0 + w.cols])

    def check(self, plan: Any) -> None:
        return None


def assembled(repo: Any, bounds: Any, needed: Any, windowed: bool) -> Any:
    plan = mosaic.plan_mosaic(repo.footprints(), bounds, needed)
    repo.check(plan)
    load_window = repo.load_window if windowed else None
    return mosaic.assemble(plan, repo.load, needed, load_window=load_window)


def value_entries() -> None:
    reach = Reach(line=((500_000.0, 6_599_980.0), (500_040.0, 6_599_980.0)), at=20.0,
                  uncertainty=5.0, corridor=10.0)  # fmt: skip
    dom = domain.DomainPolygon(polygon=shapely.box(499_990.0, 6_599_940.0, 500_050.0, 6_600_010.0),
                               crs=UTM)  # fmt: skip
    box = Bounds(x_min=500_000.0, y_min=6_599_950.0, x_max=500_050.0, y_max=6_600_000.0)
    for vname, values in VALUES.items():
        for dtype in (np.float32, np.float64):
            for nodata in (-32767.0, None):
                case = f"{vname}, {np.dtype(dtype).name}, sentinel {nodata}"
                tile = tile_of(values, dtype, nodata)
                windows = target_grid.TileWindows(tile)
                for h, rows, cols, label in ((5, 11, 9, "half cells"), (10, 6, 5, "on nodes")):
                    grid = target_grid.TargetGrid(
                        crs=UTM, spacing=h, row0=-6_600_000 // h, col0=500_000 // h, rows=rows,
                        cols=cols)  # fmt: skip
                    for threads in (1, 2):
                        line("target_grid.resample", f"{case}, {label}, {threads} threads",
                             target_grid.resample, grid, windows, threads)  # fmt: skip
                    line("target_grid.check_point_blocks", f"{case}, {label}",
                         target_grid.check_point_blocks, grid, windows, dom, 1)  # fmt: skip
                line("burn.burn_reach", case, burn.burn_reach, tile, reach)
                east = {(i, j - 3): v for (i, j), v in values.items() if j >= 3}
                repo = Memory({"east": tile_of(east, dtype, nodata, x_min=500_030.0), "west": tile})
                for bounds, blabel in ((None, "no box"), (box, "box")):
                    label = f"{case}, {blabel}"
                    line("mosaic.assemble", label, assembled, repo, bounds, None, False)
                    windowed = f"{label}, windowed"
                    line("mosaic.assemble", windowed, assembled, repo, bounds, None, True)


def cone(spacing: float = 50.0, n: int = 161) -> DemTile:
    """A cone 8 km across at 50 m, its peak at the centre node."""
    r, c = np.indices((n, n)).astype(np.float64)
    z = 2000.0 - spacing * np.hypot(r - n // 2, c - n // 2) / 10.0
    return DemTile(meta=meta(rows=n, cols=n, delta_x=spacing, delta_y=spacing),
                   array=z.astype(np.float32))  # fmt: skip


def delineate(repo: Any, fields: dict[str, Any]) -> Any:
    return catchment.delineate(catchment.CatchmentRequest(seed_crs=UTM, **fields), repo)


def catchment_entries() -> None:
    peak = cone()
    m = peak.meta
    x, y = m.x_min + 80 * m.delta_x, m.y_max - 80 * m.delta_y
    repo = Memory({"cone": peak})
    lake = shapely.Point(x + 600.0, y).buffer(120.0)
    reach = Reach(line=((x + 300.0, y + 10.0), (x + 1500.0, y + 10.0)), at=600.0, uncertainty=50.0)
    cases: dict[str, dict[str, Any]] = {
        "point on a slope": {"seed": (x + 500.0, y + 3.0)},
        "point on a node": {"seed": (x + 500.0, y)},
        "lake": {"seed": (x + 600.0, y), "lakes": (lake,), "lakes_crs": UTM},
        "reach": {"seed": (x + 900.0, y + 10.0), "reach": reach},
        "point outside": {"seed": (m.x_min - 5000.0, m.y_max)},
        "outline tolerance 0": {"seed": (x + 500.0, y + 3.0), "outline_tolerance": 0.0},
    }
    for name, fields in cases.items():
        line("catchment.delineate", name, delineate, repo, fields)
    # A valley falling west: the catchment of a point near its foot runs past
    # the first window, so the window grows (`_joined`).
    r, c = np.indices((161, 161)).astype(np.float64)
    z = 100.0 + 2.5 * c + 10.0 * np.abs(r - 80)
    valley = DemTile(meta=meta(rows=161, cols=161, delta_x=50.0, delta_y=50.0),
                     array=z.astype(np.float32))  # fmt: skip
    foot = {"seed": (m.x_min + 20 * m.delta_x, y)}
    line("catchment.delineate", "valley, the window grows", delineate, Memory({"v": valley}), foot)


def fixture_entries() -> None:
    for folder in ("seam", "lattices"):
        repo = TiffDemRepository.from_directory(ROOT / "tests" / "fixtures" / "dtm10" / folder)
        foot = repo.footprints()
        first = foot[0].meta
        x1 = first.x_min + (first.cols - 1) * first.delta_x
        boxes: dict[str, Any] = {
            "no box": None,
            "inside the first": Bounds(x_min=first.x_min + 100.5, y_min=first.y_max - 900.0,
                                       x_max=x1 - 103.0, y_max=first.y_max - 37.0),
            "across both": Bounds(x_min=min(f.meta.x_min for f in foot) + 15.0,
                                  y_min=min(f.meta.y_max for f in foot) - 1200.0,
                                  x_max=max(f.meta.x_min for f in foot) + 1400.0,
                                  y_max=max(f.meta.y_max for f in foot) - 5.0),
            "meets no tile": Bounds(x_min=0.0, y_min=0.0, x_max=1.0, y_max=1.0),
        }  # fmt: skip
        for bname, b in boxes.items():
            needed = None if b is None else shapely.box(b.x_min, b.y_min, b.x_max, b.y_max)
            shrunk = None if needed is None else needed.buffer(-2.0)
            for nlabel, nd in (("", None), (", needed", shrunk)):
                case = f"{folder}, {bname}{nlabel}"
                line("mosaic.plan_mosaic", case, mosaic.plan_mosaic, foot, b, nd)
                line("mosaic.assemble", f"{case}, whole", assembled, repo, b, nd, False)
                line("mosaic.assemble", f"{case}, windowed", assembled, repo, b, nd, True)


def opened(folder: str, fields: dict[str, Any]) -> Any:
    """`open_dem` on a DEM fixture folder: the plan, the mosaic or the target
    grid, and on the reprojected path its check points, read out."""
    source = ROOT / "tests" / "fixtures" / "dtm10" / folder
    return dem_input.open_dem(dem_input.DemRequest(sources=(source,), **fields))


def open_entries() -> None:
    corner = (49_800.0, 6_467_100.0)  # inside the seam set's first tile
    square = shapely.box(corner[0], corner[1], corner[0] + 900.5, corner[1] + 700.25)
    tilted = shapely.Polygon([(49_810.0, 6_467_120.0), (50_690.0, 6_467_300.0),
                              (50_500.0, 6_467_790.0), (49_830.0, 6_467_640.0)])  # fmt: skip
    cases: dict[str, dict[str, Any]] = {
        "box": {"bounds": Bounds(x_min=corner[0], y_min=corner[1], x_max=corner[0] + 1200.0,
                                 y_max=corner[1] + 800.0)},
        "domain, square": {"domain": domain.DomainPolygon(polygon=square, crs=UTM)},
        "domain, tilted": {"domain": domain.DomainPolygon(polygon=tilted, crs=UTM)},
        "domain, reprojected to UTM 32": {"domain": domain.DomainPolygon(polygon=tilted, crs=UTM),
                                          "target_crs": "EPSG:25832"},
    }  # fmt: skip
    for name, fields in cases.items():
        line("dem_input.open_dem", f"seam, {name}", opened, "seam", fields)


class Sink:
    def catchment(self, station: Any, result: Any) -> None: ...

    def row(self, row: Any) -> None: ...


def batch(repo: Any) -> Any:
    request = catchment_batch.BatchRequest()
    return asyncio.run(catchment_batch.run_batch(request, repo, [], UTM, [], UTM, None, Sink()))


if __name__ == "__main__":
    grid_entries()
    value_entries()
    catchment_entries()
    fixture_entries()
    open_entries()
    # No station: run_batch builds its tile boxes and reads nothing else.
    line("catchment_batch.run_batch", "no station", batch, Memory({"cone": cone()}))
