"""F2 and F7: planning a fetch, pure (23a-2, `tin_engine.fetch.plan`).

`docs/increments/23-basin-scale.md`, "Planning (`fetch/plan.py`, pure)" and
23a-2's F2 and F7. The oracles share no code with the planner:

- **The source box** must contain the domain's box in the frame, grown by
  `margin` times the source's north-south spacing in metres, sampled densely
  and moved into the source CRS with pyproj here, then grown by two source
  cells (the containment is checked against the source box shrunk by those
  two cells). The spacing of a geographic source is taken as
  `dy * 111,000 m` (the planner may use any figure at least that large).
- **The blocks** of a box are those with a node strictly within one cell of
  it (the box snapped outward to the lattice, as 23a-1's windows are),
  decided node by node from the header's affine numbers and each block's
  pixel rectangle (`cog_fixtures.block_rectangles`).
- **The ranges** are checked against the missing blocks' spans alone.

API pinned in `fetch_fixtures`' docstring. HOW THIS FILE GOES RED: every test
reaches `tin_engine.fetch.plan` through the `plan` fixture, so each errors at
setup with `ModuleNotFoundError` until 23a-2 lands, and the rest collects.
"""

from __future__ import annotations

import importlib
import io
import itertools
from collections.abc import Sequence
from types import ModuleType
from typing import Any

import numpy as np
import pyproj
import pytest
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient

from cog_fixtures import block_rectangles
from crs_fixtures import proj4_of, refuse_point_moves
from fetch_fixtures import (
    GAP,
    MAX_RANGE,
    MIB,
    PROJECTED_CRS,
    block_span,
    geographic,
    page_of,
    projected,
    wide,
    with_long_header,
)
from geotiff_fixtures import TIE_X, TIE_Y
from tin_engine.domain import DomainPolygon
from tin_engine.io.geotiff import read_page
from tin_engine.mosaic import Bounds

METRES_PER_DEGREE = 111_000.0


@pytest.fixture(scope="module")
def plan() -> ModuleType:
    return importlib.import_module("tin_engine.fetch.plan")


@pytest.fixture(scope="module")
def fetch_error() -> type[Exception]:
    error: type[Exception] = importlib.import_module("tin_engine.fetch.http").FetchError
    return error


def header(data: bytes, *, geographic: bool = False) -> tuple[Any, Any, Any]:
    """(meta, dtype, page) by the reader; the geographic flag is 23a-2's."""
    extra = {"geographic": True} if geographic else {}
    return read_page(io.BytesIO(data), nodata=None, **extra)


def box_of(x0: float, y0: float, x1: float, y1: float) -> Bounds:
    return Bounds(x_min=x0, y_min=y0, x_max=x1, y_max=y1)


def domain(ring: Sequence[tuple[float, float]], crs: str) -> DomainPolygon:
    return DomainPolygon(polygon=orient(Polygon(ring), 1.0), crs=crs)


def brute_force_blocks(page: Any, meta: Any, box: Bounds) -> tuple[int, ...]:
    """Blocks with a node strictly within one cell of `box`, node by node."""
    xs = meta.x_min + meta.delta_x * np.arange(meta.cols)
    ys = meta.y_max - meta.delta_y * np.arange(meta.rows)
    near_x = (xs > box.x_min - meta.delta_x) & (xs < box.x_max + meta.delta_x)
    near_y = (ys > box.y_min - meta.delta_y) & (ys < box.y_max + meta.delta_y)
    return tuple(
        i
        for i, (r0, r1, c0, c1) in enumerate(block_rectangles(page))
        if near_y[r0:r1].any() and near_x[c0:c1].any()
    )


def dense(box: Bounds, grow: float, n: int = 41) -> np.ndarray:
    """An n x n grid of points over `box` grown by `grow`, edges included."""
    xs = np.linspace(box.x_min - grow, box.x_max + grow, n)
    ys = np.linspace(box.y_min - grow, box.y_max + grow, n)
    gx, gy = np.meshgrid(xs, ys)
    return np.column_stack([gx.ravel(), gy.ravel()])


def inside(points: np.ndarray, box: Bounds, shrink_x: float, shrink_y: float) -> bool:
    x, y = points[:, 0], points[:, 1]
    return bool(
        np.all(x >= box.x_min + shrink_x)
        and np.all(x <= box.x_max - shrink_x)
        and np.all(y >= box.y_min + shrink_y)
        and np.all(y <= box.y_max - shrink_y)
    )


# --------------------------------------------------------------------------
# The request
# --------------------------------------------------------------------------


class TestTheRequest:
    def test_exactly_one_of_domain_and_box(self, plan: ModuleType) -> None:
        ring = [(TIE_X + 50, TIE_Y - 50), (TIE_X + 200, TIE_Y - 50), (TIE_X + 100, TIE_Y - 150)]
        with pytest.raises(ValueError, match="domain"):
            plan.FetchRequest(
                source="s",
                domain=domain(ring, PROJECTED_CRS),
                box=box_of(TIE_X, TIE_Y - 100, TIE_X + 100, TIE_Y),
            )
        with pytest.raises(ValueError, match="domain"):
            plan.FetchRequest(source="s")

    def test_the_defaults(self, plan: ModuleType) -> None:
        request = plan.FetchRequest(source="s", box=box_of(0, 0, 1, 1))
        assert (request.margin, request.connections) == (4, 8)
        assert (request.out_crs, request.dry_run, request.refresh) == (None, False, False)


# --------------------------------------------------------------------------
# F2: the source box
# --------------------------------------------------------------------------


class TestF2TheSourceBox:
    def test_in_the_source_crs_the_box_grows_by_margin_and_two_cells(
        self, plan: ModuleType
    ) -> None:
        meta = header(projected())[0]
        box = box_of(TIE_X + 101.3, TIE_Y - 187.7, TIE_X + 333.9, TIE_Y - 61.1)
        request = plan.FetchRequest(source="s", box=box, margin=4)
        got = plan.source_box(request, meta)
        grow = 4 * meta.delta_y
        points = dense(box, grow)
        assert inside(points, got, 2 * meta.delta_x * 0.999, 2 * meta.delta_y * 0.999)

    @pytest.mark.parametrize("margin", [0, 4, 9])
    def test_a_reprojected_domain_box_grown_by_margin_lies_inside(
        self, plan: ModuleType, margin: int
    ) -> None:
        """The frame is EPSG:3035; the source is EPSG:25833."""
        meta = header(projected())[0]
        to_frame = pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:3035", always_xy=True)
        corners = [
            (TIE_X + 120, TIE_Y - 200),
            (TIE_X + 600, TIE_Y - 180),
            (TIE_X + 300, TIE_Y - 30),
        ]
        ring = [to_frame.transform(x, y) for x, y in corners]
        given = domain(ring, "EPSG:3035")
        request = plan.FetchRequest(source="s", domain=given, out_crs="EPSG:3035", margin=margin)
        got = plan.source_box(request, meta)
        frame_box = box_of(*given.polygon.bounds)
        points = dense(frame_box, margin * meta.delta_y)
        back = pyproj.Transformer.from_crs("EPSG:3035", PROJECTED_CRS, always_xy=True)
        moved = np.column_stack(back.transform(points[:, 0], points[:, 1]))
        assert inside(moved, got, 2 * meta.delta_x * 0.999, 2 * meta.delta_y * 0.999)

    def test_out_crs_spelt_as_the_sources_proj_string_is_no_out_crs(
        self, plan: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Audit PR B, red test 10: `out_crs` as the PROJ string of the
        source's EPSG code (25833) is the source's CRS (`crs.same_crs`), so the box is
        the one with no `out_crs`, and no bounds are transformed."""
        meta = header(projected())[0]
        assert meta.crs == f"EPSG:{meta.epsg}" == PROJECTED_CRS
        box = box_of(TIE_X + 101.3, TIE_Y - 187.7, TIE_X + 333.9, TIE_Y - 61.1)
        bare = plan.source_box(plan.FetchRequest(source="s", box=box, margin=4), meta)
        refuse_point_moves(monkeypatch)
        spelt = plan.FetchRequest(source="s", box=box, margin=4, out_crs=proj4_of(meta.epsg))
        assert plan.source_box(spelt, meta) == bare

    def test_a_geographic_source_grows_by_its_spacing_in_metres(self, plan: ModuleType) -> None:
        """Source EPSG:4326 at 0.001°; the frame is UTM 33N, in metres."""
        meta = header(geographic(), geographic=True)[0]
        to_frame = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32633", always_xy=True)
        corners = [(15.005, 59.47), (15.03, 59.48), (15.02, 59.495)]
        ring = [to_frame.transform(lon, lat) for lon, lat in corners]
        given = domain(ring, "EPSG:32633")
        request = plan.FetchRequest(source="s", domain=given, out_crs="EPSG:32633", margin=4)
        got = plan.source_box(request, meta)
        points = dense(box_of(*given.polygon.bounds), 4 * meta.delta_y * METRES_PER_DEGREE)
        back = pyproj.Transformer.from_crs("EPSG:32633", "EPSG:4326", always_xy=True)
        moved = np.column_stack(back.transform(points[:, 0], points[:, 1]))
        assert inside(moved, got, 2 * meta.delta_x * 0.999, 2 * meta.delta_y * 0.999)

    def test_a_box_crossing_the_antimeridian_is_refused(
        self, plan: ModuleType, fetch_error: type[Exception]
    ) -> None:
        meta = header(geographic(), geographic=True)[0]
        request = plan.FetchRequest(source="s", box=box_of(179.99, 10, 180.5, 11))
        with pytest.raises(fetch_error):
            plan.source_box(request, meta)


# --------------------------------------------------------------------------
# F2: the blocks
# --------------------------------------------------------------------------


def boxes_over(meta: Any, count: int, seed: int) -> list[Bounds]:
    """Boxes off the lattice: inside, across the edges, and thin."""
    rng = np.random.default_rng(seed)
    width, height = (meta.cols - 1) * meta.delta_x, (meta.rows - 1) * meta.delta_y
    out = []
    for _ in range(count):
        x0, x1 = sorted(rng.uniform(-0.2, 1.2, 2) * width + meta.x_min)
        y1, y0 = sorted(meta.y_max - rng.uniform(-0.2, 1.2, 2) * height, reverse=True)
        if x1 - x0 < 0.3 * meta.delta_x or y1 - y0 < 0.3 * meta.delta_y:
            continue
        overlaps = (
            x1 > meta.x_min
            and x0 < meta.x_min + width
            and y0 < meta.y_max
            and y1 > meta.y_max - height
        )
        if overlaps:
            out.append(box_of(x0, y0, x1, y1))
    out.append(
        box_of(
            meta.x_min + 3.37 * meta.delta_x,
            meta.y_max - 3.61 * meta.delta_y,
            meta.x_min + 3.41 * meta.delta_x,
            meta.y_max - 3.55 * meta.delta_y,
        )
    )
    return out


class TestF2TheBlocks:
    @pytest.mark.parametrize("which", ["projected", "wide", "geographic"])
    def test_the_blocks_equal_a_brute_force_over_every_block(
        self, plan: ModuleType, which: str
    ) -> None:
        data = {"projected": projected, "wide": wide, "geographic": geographic}[which]()
        meta, _, page = header(data, geographic=which == "geographic")
        boxes = boxes_over(meta, 40, seed=len(which))
        assert len(boxes) > 20
        for box in boxes:
            got = plan.plan_object("obj", "http://x/obj.tif", page, meta, box, ())
            assert tuple(got.blocks) == brute_force_blocks(page, meta, box), box

    def test_the_brute_force_can_fail(self) -> None:
        """A box one cell wide meets fewer blocks than one four blocks wide."""
        meta, _, page = header(wide())
        thin = box_of(
            meta.x_min + 5.5 * meta.delta_x,
            meta.y_max - 30,
            meta.x_min + 6.5 * meta.delta_x,
            meta.y_max,
        )
        broad = box_of(
            meta.x_min + 5.5 * meta.delta_x,
            meta.y_max - 30,
            meta.x_min + 70 * meta.delta_x,
            meta.y_max,
        )
        assert len(brute_force_blocks(page, meta, thin)) < len(
            brute_force_blocks(page, meta, broad)
        )

    def test_a_box_outside_the_raster_is_refused(
        self, plan: ModuleType, fetch_error: type[Exception]
    ) -> None:
        meta, _, page = header(projected())
        far = box_of(TIE_X - 5000, TIE_Y + 1000, TIE_X - 4000, TIE_Y + 2000)
        with pytest.raises(fetch_error):
            plan.plan_object("obj", "http://x/obj.tif", page, meta, far, ())


# --------------------------------------------------------------------------
# F2: the ranges
# --------------------------------------------------------------------------


def check_ranges(ranges: Sequence[tuple[int, int]], spans: Sequence[tuple[int, int]]) -> None:
    """Every span in exactly one range; each range starts and ends on a span
    it holds; no gap between held spans over 64 KiB; no range over 8 MiB
    unless it holds one span."""
    ranges = [tuple(r) for r in ranges]
    assert ranges == sorted(ranges)
    held: dict[tuple[int, int], list[tuple[int, int]]] = {r: [] for r in ranges}
    for span in spans:
        owners = [r for r in ranges if r[0] <= span[0] and span[1] <= r[1]]
        assert len(owners) == 1, (span, owners)
        held[owners[0]].append(span)
    for (start, stop), inner in held.items():
        inner.sort()
        assert inner, (start, stop)
        assert (start, stop) == (inner[0][0], max(s[1] for s in inner))
        for (_, end), (begin, _) in itertools.pairwise(inner):
            assert begin - end <= GAP, (start, stop)
        assert stop - start <= MAX_RANGE or len(inner) == 1, (start, stop)


class TestF2TheRanges:
    def test_ranges_cover_exactly_the_missing_blocks(self, plan: ModuleType) -> None:
        meta, _, page = header(wide())
        box = box_of(meta.x_min + 100.5, meta.y_max - 200.5, meta.x_min + 900.5, meta.y_max - 3.3)
        rng = np.random.default_rng(5)
        everything = plan.plan_object("obj", "http://x/obj.tif", page, meta, box, ())
        present = tuple(int(i) for i in rng.choice(everything.blocks, 5, replace=False))
        got = plan.plan_object("obj", "http://x/obj.tif", page, meta, box, present)
        assert got.blocks == everything.blocks
        missing = [i for i in got.blocks if i not in present]
        check_ranges(got.ranges, [block_span(page, i) for i in missing])
        assert got.bytes == sum(b - a for a, b in got.ranges)
        assert len(got.ranges) > 1  # block rows lie 100 KiB apart

    def test_the_rule_can_fail(self) -> None:
        with pytest.raises(AssertionError):
            check_ranges([(0, 10), (10, 20)], [(0, 20)])
        with pytest.raises(AssertionError):
            check_ranges([(0, 201 + GAP)], [(0, 100), (101 + GAP, 201 + GAP)])

    def test_coalesce_at_exactly_64_kib_and_one_byte_more(self, plan: ModuleType) -> None:
        assert plan.coalesce([(0, 100), (100 + GAP, 200 + GAP)]) == ((0, 200 + GAP),)
        assert plan.coalesce([(0, 100), (101 + GAP, 200 + GAP)]) == (
            (0, 100),
            (101 + GAP, 200 + GAP),
        )

    def test_coalesce_caps_a_range_at_8_mib_but_not_one_block(self, plan: ModuleType) -> None:
        spans = [(i * MIB, (i + 1) * MIB) for i in range(20)]
        got = plan.coalesce(spans)
        check_ranges(got, spans)
        assert len(got) == 3
        assert plan.coalesce([(0, 9 * MIB)]) == ((0, 9 * MIB),)

    def test_coalesce_property(self, plan: ModuleType) -> None:
        rng = np.random.default_rng(17)
        for _ in range(200):
            spans, at = [], int(rng.integers(0, 1000))
            for _ in range(int(rng.integers(1, 40))):
                size = int(rng.choice([1, 300, 70_000, MIB, 3 * MIB, 9 * MIB]))
                spans.append((at, at + size))
                at += size + int(rng.choice([0, 1, GAP - 1, GAP, GAP + 1, 5 * MIB]))
            order = rng.permutation(len(spans))
            check_ranges(plan.coalesce([spans[i] for i in order]), spans)


# --------------------------------------------------------------------------
# F7: the header prefix
# --------------------------------------------------------------------------


class TestF7ThePrefix:
    def test_a_prefix_short_of_the_offsets_is_incomplete(self, plan: ModuleType) -> None:
        data = with_long_header(3 * MIB)
        assert plan.parse_prefix(data[:MIB], nodata=None) is None
        assert plan.parse_prefix(data[: 2 * MIB], nodata=None) is None
        meta, _, page = plan.parse_prefix(data[: 4 * MIB], nodata=None)
        assert tuple(page.dataoffsets) == tuple(page_of(data).dataoffsets)
        assert tuple(page.databytecounts) == tuple(page_of(data).databytecounts)
        assert (meta.rows, meta.cols) == (48, 64)

    def test_a_prefix_cutting_the_offset_array_is_incomplete(self, plan: ModuleType) -> None:
        data = with_long_header(3 * MIB)
        offsets_at = page_of(data).tags[324].valueoffset
        for cut in (offsets_at + 4, offsets_at + 20):
            assert plan.parse_prefix(data[:cut], nodata=None) is None, cut

    def test_a_prefix_cutting_the_first_ifd_is_incomplete(self, plan: ModuleType) -> None:
        assert plan.parse_prefix(projected()[:24], nodata=None) is None

    def test_a_geographic_header_is_read(self, plan: ModuleType) -> None:
        meta, _, _ = plan.parse_prefix(geographic()[:MIB], nodata=None)
        assert meta.epsg == 4326

    def test_read_page_reads_geographic_with_or_without_the_flag(self) -> None:
        """Amended by 15c-2 (D6): the mesh path's 2048 refusal (23a-1 W9) is
        lifted, so `read_page` reads a geographic header without the flag too,
        to the same meta; fetch's `geographic=True` call still reads it."""
        with_flag = header(geographic(), geographic=True)[0]
        assert (with_flag.epsg, with_flag.geographic, with_flag.crs) == (4326, True, "EPSG:4326")
        assert header(geographic())[0] == with_flag
