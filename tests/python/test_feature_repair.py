"""The land-cover stage of ``open_features``: repair, merge and simplify (20c-3).

``docs/increments/20c-soft-quality.md``, "Design of PR 20c-3" (steps 1 to 3,
"Where it sits", "Off") and "Tests ``@tester`` writes red first", 20c-3: RP1
to RP8, and "Steps 2 and 3". Not mutation-critical (the design's words: the
geometry is GEOS's; what 20c-3 owns is which polygons go in, with which
tolerance, in which order, and each of those is an assertion here).

Every fixture is a hand-made coverage of a few polygons in EPSG:25833 at the
UTM-shaped offset of ``feature_fixtures`` (``at``), under the ``corine`` map
(``Code_18``), so the stage applies; no data file except RP8's, which is the
committed extract (``tests/fixtures/corine``), not ``rasputin_data``.

Scale of the bounds: coordinates near 5e5 and 6.6e6 m, where a double
resolves about 1e-9 m; the fixtures span at most 200 m, except RP7's
10 km polygons. The 1e-6 m² area bound (RP1) was checked against shapely
2.2.0 / GEOS 3.14.1, where the moved area matched the wedge to 1e-8 m².

PINNED HERE, where the design leaves it open (listed for ``@architect``):

- ``FeatureRequest`` carries the four switches as ``repair_m: float``,
  ``merge_same_class: bool``, ``tolerance_m: float`` and
  ``outline_snap_m: float``, each defaulting to OFF (0.0, False, 0.0,
  0.0): the library default is today's path, and the CLI supplies the
  run defaults (0.05, on, 0, 5). So every older ``open_features`` test,
  which builds a request without them, keeps today's output.
- The label polygon after the stage is ``TerrainFeature.polygon``, in
  source order; with the merge off, one feature per input polygon.
- An overlap is given to the smaller polygon (16c's D2, "the smallest
  wins"), which is ``coverage_clean``'s ``merge_strategy="min_area"``; its
  default, ``longest_border``, empties a lake lying inside an unholed
  forest (checked on GEOS 3.14.1: area 0). See ``TestOverlap``.
- RP1's third polygon has a straight side through both slit corners, so the
  "two thin triangles the moved corner sweeps" on it have zero area here,
  and the moved area is the wedge's alone.
- A merged polygon's feature keeps the first fid of its class in source
  order; the order of the remaining features is not asserted.

RED at the commit that adds this file: ``FeatureRequest`` forbids extra
fields, so every request with a switch is refused by Pydantic.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon
from shapely.geometry.base import BaseGeometry

import feature_fixtures as ff
from feature_fixtures import UTM33, X0, Y0, Feat, at, domain_of, write_geojson
from gpkg_fixtures import EXTRACT
from test_cli_mesh_domain import quarter_circle

S = 0.05  # question 9's default repair tolerance, metres (ruling 9)
L = 79.13  # the slit's border length (M5, in local coordinates)
WIDTH = 0.01  # the slit's open end, 1 cm
WEDGE = 0.5 * L * WIDTH  # its area, 0.39565 m²
#: All four switches off: "the off switch the gates use".
OFF = {"repair_m": 0.0, "merge_same_class": False, "tolerance_m": 0.0, "outline_snap_m": 0.0}


@pytest.fixture(scope="module")
def fi() -> ModuleType:
    return ff.feature_input()


def coded(fid: Any, geometry: BaseGeometry, code: str) -> Feat:
    return Feat(fid, geometry, {"Code_18": code})


def poly(*points: tuple[float, float]) -> Polygon:
    return Polygon([at(x, y) for x, y in points])


def box(x0: float, y0: float, x1: float, y1: float) -> Polygon:
    return poly((x0, y0), (x1, y0), (x1, y1), (x0, y1))


#: Holds every fixture below at least 10 m from its outline, so the read
#: region (its hull grown by 100 m) clips nothing and no border is near it.
BIG = domain_of(box(-70, -70, 200, 70))


def request(fi: ModuleType, *paths: Path, map_name: str = "corine", **switches: Any) -> Any:
    sources = tuple(ff.source(p, map_name) for p in paths)
    return fi.FeatureRequest(sources=sources, **{**OFF, **switches})


def opened(
    fi: ModuleType, tmp_path: Path, features: list[Feat], domain: Any = BIG, **switches: Any
) -> Any:
    path = write_geojson(tmp_path / "cover.geojson", features)
    return fi.open_features(request(fi, path, **switches), domain, UTM33)


def labels(fs: Any) -> list[BaseGeometry]:
    return [f.polygon for f in fs.features]


def ring_edges(polygons: list[BaseGeometry]) -> np.ndarray:
    """Every ring edge's length, over all polygons."""
    out = []
    for p in polygons:
        for ring in shapely.get_rings(p):
            xy = shapely.get_coordinates(ring)
            out.append(np.hypot(*np.diff(xy, axis=0).T))
    return np.concatenate(out)


def holes(polygons: list[BaseGeometry]) -> list[Polygon]:
    union = shapely.union_all(polygons)
    return [Polygon(r) for part in shapely.get_parts(union) for r in part.interiors]


def all_xy(fs: Any) -> np.ndarray:
    """Every vertex of every label polygon and every line."""
    parts = [shapely.get_coordinates(g) for g in labels(fs) if g is not None]
    parts += [shapely.get_coordinates(line) for f in fs.features for line in f.lines]
    return np.concatenate(parts)


# ---------------------------------------------------------------- RP1, RP6


def slit(extra: tuple[float, float] | None = None) -> list[Polygon]:
    """M5's slit along the x axis: A (above) and B (below) share the far
    vertex F = (0, 0); A's border runs to (L, 1 cm), B's to (L, 0); C, 50 m
    tall on each side, closes the wedge with its 1 cm edge; D lies west of
    F, so F is inside the coverage and the wedge is a hole. ``extra`` is one
    more vertex on A's border (RP6)."""
    a = [(0.0, 0.0), *([extra] if extra else []), (L, WIDTH), (L, 50.0), (0.0, 50.0)]
    return [
        poly(*a),
        poly((0, 0), (0, -50), (L, -50), (L, 0)),
        poly((L, -50), (L + 50, -50), (L + 50, 50), (L, 50), (L, WIDTH), (L, 0)),
        poly((-50, -50), (0, -50), (0, 0), (0, 50), (-50, 50)),
    ]


def slit_features(extra: tuple[float, float] | None = None) -> list[Feat]:
    return [coded(i + 1, p, c) for i, (p, c) in enumerate(zip(slit(extra), CODES, strict=True))]


CODES = ("311", "312", "211", "231")


class TestTheSlit:
    """RP1, at S = 0.05 and at S = 0."""

    def test_rp1_the_premise_the_slit_is_a_hole(self) -> None:
        (hole,) = holes(slit())
        assert hole.area == pytest.approx(WEDGE, abs=1e-6)  # 1 mm²; see the docstring

    def test_rp1_at_5_cm_the_slit_is_repaired(self, fi: ModuleType, tmp_path: Path) -> None:
        fs = opened(fi, tmp_path, slit_features(), repair_m=S)
        assert [f.fid for f in fs.features] == [1, 2, 3, 4]
        got = labels(fs)
        assert holes(got) == []
        shared = shapely.intersection(got[0], got[1])
        assert shapely.length(shared) == pytest.approx(L, abs=1e-3)  # the two borders are one
        assert ring_edges(got).min() >= S  # no ring edge under 5 cm is left
        moved = sum(
            shapely.symmetric_difference(a, b).area for a, b in zip(slit(), got, strict=True)
        )
        assert moved == pytest.approx(WEDGE, abs=1e-6)  # m²; see the docstring

    def test_rp1_at_0_the_slit_stays(self, fi: ModuleType, tmp_path: Path) -> None:
        """S = 0 with the merge on, so the stage runs (the design's "S = 0
        still runs ``coverage_clean`` at zero"): the wedge is still a hole,
        and nothing moved."""
        fs = opened(fi, tmp_path, slit_features(), repair_m=0.0, merge_same_class=True)
        got = labels(fs)
        (hole,) = holes(got)
        assert hole.area == pytest.approx(WEDGE, abs=1e-6)  # 1 mm²; see the docstring
        assert all(shapely.equals(a, b) for a, b in zip(slit(), got, strict=True))


class TestEmptied:
    """Ruling G6 (e): RP1's degenerate case. C, 1 cm thick between the slit
    and a fifth polygon E east of it, is given wholly to its neighbours by
    the repair at S; it is counted ``empty``, not ``outside`` (it is not
    outside the domain). At S = 0 (the merge on, so the stage runs) C stays:
    the premise."""

    THIN = 0.01  # C's thickness, metres

    def features(self) -> list[Feat]:
        t = self.THIN
        polygons = [
            *slit()[:2],
            poly((L, -50), (L + t, -50), (L + t, 50), (L, 50), (L, WIDTH), (L, 0)),
            slit()[3],
            box(L + t, -50, L + 50, 50),
        ]
        return [
            coded(i + 1, p, c)
            for i, (p, c) in enumerate(zip(polygons, (*CODES, "121"), strict=True))
        ]

    def test_g6e_the_premise_at_0_c_stays(self, fi: ModuleType, tmp_path: Path) -> None:
        fs = opened(fi, tmp_path, self.features(), repair_m=0.0, merge_same_class=True)
        assert [f.fid for f in fs.features] == [1, 2, 3, 4, 5]
        assert (fs.outside, fs.empty) == (0, 0)

    def test_g6e_the_emptied_polygon_is_counted_empty(self, fi: ModuleType, tmp_path: Path) -> None:
        fs = opened(fi, tmp_path, self.features(), repair_m=S)
        assert [f.fid for f in fs.features] == [1, 2, 4, 5]  # the premise: C emptied
        assert (fs.outside, fs.empty) == (0, 1)


class TestTheOrder:
    """RP6: repair before simplifying. A vertex 4.5 mm off A's straight
    border, halfway along and away from the slit, with
    ``--features-tolerance 2``: repaired first, the two borders are one inner
    border and the simplification drops it; simplified first, the border is
    a side of the gap, which ``simplify_boundary=False`` keeps."""

    EXTRA = (L / 2, WIDTH / 2 + 0.0045)

    def nearest(self, xy: np.ndarray) -> float:
        return float(
            np.hypot(xy[:, 0] - (X0 + self.EXTRA[0]), xy[:, 1] - (Y0 + self.EXTRA[1])).min()
        )

    def test_rp6_the_premise_simplified_first_keeps_the_vertex(self) -> None:
        polygons = np.array(slit(self.EXTRA), dtype=object)
        simplified = shapely.coverage_simplify(polygons, 2.0, simplify_boundary=False)
        cleaned = shapely.coverage_clean(simplified, snapping_distance=S, gap_width=S)
        assert self.nearest(np.concatenate([shapely.get_coordinates(g) for g in cleaned])) < 0.01

    def test_rp6_repaired_first_drops_it(self, fi: ModuleType, tmp_path: Path) -> None:
        fs = opened(fi, tmp_path, slit_features(self.EXTRA), repair_m=S, tolerance_m=2.0)
        assert self.nearest(all_xy(fs)) > 1.0


# ---------------------------------------------------------------- RP2-RP5


def same(fs: Any, inputs: list[Polygon]) -> None:
    got = labels(fs)
    assert len(got) == len(inputs)
    for a, b in zip(inputs, got, strict=True):
        assert shapely.equals(a, b), (a.wkt, b.wkt)


class TestWhatStays:
    def test_rp2_a_real_gap_wider_than_s_stays(self, fi: ModuleType, tmp_path: Path) -> None:
        """A strip 2 m wide between two polygons, closed at both ends."""
        inputs = [
            box(0, -50, 100, 0),
            box(0, 2, 100, 50),
            box(-50, -50, 0, 50),
            box(100, -50, 150, 50),
        ]
        (before,) = holes(inputs)
        assert before.area == pytest.approx(200.0, rel=1e-9)  # the premise: a hole
        fs = opened(
            fi,
            tmp_path,
            [coded(i, p, c) for i, (p, c) in enumerate(zip(inputs, CODES, strict=True))],
            repair_m=S,
        )
        (after,) = holes(labels(fs))
        assert after.area == pytest.approx(before.area, rel=1e-9)
        same(fs, inputs)

    def test_rp3_nothing_to_repair(self, fi: ModuleType, tmp_path: Path) -> None:
        """A valid coverage with a bent border and no near miss under S."""
        inputs = [
            poly((0, 0), (40.3, 3.1), (80.7, -2.2), (80.7, 60), (0, 60)),
            poly((0, 0), (0, -50), (80.7, -50), (80.7, -2.2), (40.3, 3.1)),
        ]
        assert shapely.coverage_is_valid(np.array(inputs, dtype=object))  # the premise
        fs = opened(
            fi, tmp_path, [coded(1, inputs[0], "311"), coded(2, inputs[1], "312")], repair_m=S
        )
        same(fs, inputs)
        assert [shapely.get_num_coordinates(g) for g in labels(fs)] == [
            shapely.get_num_coordinates(g) for g in inputs
        ]

    def test_rp4_one_source_at_a_time(self, fi: ModuleType, tmp_path: Path) -> None:
        """Two sources, one polygon each, 1 cm apart: the gap stays."""
        a, b = box(0, 0, 50, 50), box(50.01, 0, 100, 50)
        joined = shapely.coverage_clean(
            np.array([a, b], dtype=object), snapping_distance=S, gap_width=S
        )
        assert shapely.distance(*joined) == 0.0  # the premise: in one source they would join
        pa = write_geojson(tmp_path / "a.geojson", [coded(1, a, "311")])
        pb = write_geojson(tmp_path / "b.geojson", [coded(2, b, "312")])
        fs = fi.open_features(request(fi, pa, pb, repair_m=S), BIG, UTM33)
        same(fs, [a, b])
        assert shapely.distance(*labels(fs)) == pytest.approx(0.01, abs=1e-6)

    def test_rp5_a_map_without_codes_is_not_land_cover(
        self, fi: ModuleType, tmp_path: Path
    ) -> None:
        """The slit under the ``property`` map, with the repair and the merge
        on: the same lines as with everything off, and no label polygons."""
        feats = [Feat(i + 1, p, {"property": "land_cover"}) for i, p in enumerate(slit())]
        path = write_geojson(tmp_path / "p.geojson", feats)
        on = fi.open_features(
            request(fi, path, map_name="property", repair_m=S, merge_same_class=True), BIG, UTM33
        )
        off = fi.open_features(request(fi, path, map_name="property"), BIG, UTM33)
        assert [f.fid for f in on.features] == [f.fid for f in off.features] == [1, 2, 3, 4]
        for f, g in zip(on.features, off.features, strict=True):
            assert f.polygon is None and g.polygon is None
            assert len(f.lines) == len(g.lines)
            assert all(shapely.equals_exact(x, y, 0) for x, y in zip(f.lines, g.lines, strict=True))


# ---------------------------------------------------------------- RP7


class TestTheClip:
    """RP7: two polygons 10 km wide, their shared border crossing the
    domain: each label polygon lies within the read region, together they
    cover the domain, and the lines inside it are today's."""

    WEST = box(-5000, -5000, 60, 5000)
    EAST = box(60, -5000, 5000, 5000)

    def test_rp7(self, fi: ModuleType, tmp_path: Path) -> None:
        domain = domain_of(box(0, 0, 120, 50))
        feats = [coded(1, self.WEST, "311"), coded(2, self.EAST, "312")]
        on = opened(fi, tmp_path, feats, domain=domain, repair_m=S)
        region = fi.source_region(domain, UTM33, UTM33)
        got = labels(on)
        assert len(got) == 2
        for g in got:
            assert region.buffer(1e-6).covers(g)  # 1 µm: the clip's own rounding
            assert g.area < 0.01 * self.WEST.area  # it was clipped
        assert shapely.union_all(got).covers(domain.polygon)
        off = opened(fi, tmp_path, feats, domain=domain)
        for f, g in zip(on.features, off.features, strict=True):
            assert f.fid == g.fid
            # Normalised: the repair may start or orient a ring differently.
            mine = sorted(shapely.normalize(x).wkb for x in f.lines)
            today = sorted(shapely.normalize(y).wkb for y in g.lines)
            assert mine == today and len(mine) == 1


# ---------------------------------------------------------------- RP8


#: RECORDED at 985b4a47 (no 20c-3 production change) by ``digest`` below,
#: over ``open_features`` of the committed extract on the quarter circle
#: with the ``corine`` map and today's request. Integers only (fid, mask,
#: code, vertices per line, the label polygon's type and vertex count), so
#: a last bit of pyproj does not move it. No commit may update it to agree
#: with new code.
TODAY = (60, "e977f021bebf6f52a1e0e36eccc09e01351b879cf0813d1c91906c9a7bee0957", 9, 19, 0)


def digest(fs: Any) -> tuple[int, str, int, int, int]:
    rows = [
        [
            f.fid,
            int(f.mask),
            f.code,
            [int(shapely.get_num_coordinates(line)) for line in f.lines],
            None
            if f.polygon is None
            else [f.polygon.geom_type, int(shapely.get_num_coordinates(f.polygon))],
        ]
        for f in fs.features
    ]
    text = json.dumps(rows).encode()
    return len(rows), hashlib.sha256(text).hexdigest(), fs.outside, fs.clipped, fs.empty


class TestAllOff:
    """RP8: the four switches off give today's ``FeatureSet``."""

    @pytest.fixture(scope="class")
    def domain(self) -> Any:
        return domain_of(Polygon(quarter_circle()))

    def test_rp8_all_off_is_today(
        self, fi: ModuleType, domain: Any, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        def refused(*_: Any, **__: Any) -> Any:
            raise AssertionError("the stage ran with every switch off")

        for name in ("coverage_clean", "coverage_union", "coverage_simplify"):
            monkeypatch.setattr(shapely, name, refused)
        off = fi.open_features(
            fi.FeatureRequest(sources=(ff.source(EXTRACT, "corine"),), **OFF), domain, UTM33
        )
        assert digest(off) == TODAY

    def test_rp8_the_request_defaults_are_off(self, fi: ModuleType, domain: Any) -> None:
        default = fi.open_features(
            fi.FeatureRequest(sources=(ff.source(EXTRACT, "corine"),)), domain, UTM33
        )
        off = fi.open_features(
            fi.FeatureRequest(sources=(ff.source(EXTRACT, "corine"),), **OFF), domain, UTM33
        )
        assert digest(default) == digest(off)
        for f, g in zip(default.features, off.features, strict=True):
            assert (f.fid, f.mask, f.code) == (g.fid, g.mask, g.code)
            assert all(shapely.equals_exact(x, y, 0) for x, y in zip(f.lines, g.lines, strict=True))
            assert (f.polygon is None) == (g.polygon is None)
            assert f.polygon is None or shapely.equals_exact(f.polygon, g.polygon, 0)


# ---------------------------------------------------------------- steps 2 and 3


class TestMerge:
    """Step 2: neighbours of one class are one polygon, keeping the first
    fid of the class in source order; the partition stays valid."""

    def test_two_of_one_class_merge(self, fi: ModuleType, tmp_path: Path) -> None:
        a, c, b = box(0, 0, 50, 50), box(0, -50, 100, 0), box(50, 0, 100, 50)
        feats = [coded(7, a, "311"), coded(3, c, "312"), coded(5, b, "311")]
        fs = opened(fi, tmp_path, feats, repair_m=S, merge_same_class=True)
        by_fid = {f.fid: f for f in fs.features}
        assert set(by_fid) == {7, 3}
        assert shapely.equals(by_fid[7].polygon, shapely.union(a, b))
        assert by_fid[7].code == 311
        shared = LineString([at(50, 0), at(50, 50)])
        lines = shapely.union_all(list(by_fid[7].lines))
        assert shapely.length(shapely.intersection(lines, shared)) == pytest.approx(0.0, abs=1e-9)
        assert shapely.coverage_is_valid(np.array(labels(fs), dtype=object))

    def test_off_keeps_both(self, fi: ModuleType, tmp_path: Path) -> None:
        a, b = box(0, 0, 50, 50), box(50, 0, 100, 50)
        fs = opened(fi, tmp_path, [coded(7, a, "311"), coded(5, b, "311")], repair_m=S)
        assert [f.fid for f in fs.features] == [7, 5]


class TestTolerance:
    """Step 3: ``coverage_simplify`` at the tolerance, also with the repair
    at 0 (the old design's coverage check and its error are gone)."""

    #: A 100 m border with a vertex 5 cm off its straight line: a triangle of
    #: 2.5 m², under the 4 m² (tolerance squared) that ``coverage_simplify``
    #: removes (its tolerance is about the square root of the area dropped).
    BENT = (50.0, 0.05)

    def inputs(self) -> list[Polygon]:
        return [
            poly((0, 0), self.BENT, (100, 0), (100, 50), (0, 50)),
            poly((0, 0), (0, -50), (100, -50), (100, 0), self.BENT),
            box(-50, -50, 0, 50),
        ]

    @pytest.mark.parametrize("repair", [0.0, S])
    def test_tolerance_2(self, fi: ModuleType, tmp_path: Path, repair: float) -> None:
        inputs = self.inputs()
        feats = [coded(i + 1, p, c) for i, (p, c) in enumerate(zip(inputs, CODES[:3], strict=True))]
        fs = opened(fi, tmp_path, feats, repair_m=repair, tolerance_m=2.0)
        got = labels(fs)
        assert shapely.coverage_is_valid(np.array(got, dtype=object))
        bent = np.array(at(*self.BENT))
        assert np.hypot(*(all_xy(fs) - bent).T).min() > 1.0  # simplified: the bend is gone
        for before, after in zip(inputs, got, strict=True):
            assert abs(after.area - before.area) < 2.0 * before.length


# ---------------------------------------------------------------- overlaps (pinned)


class TestOverlap:
    """PINNED (see the docstring): an overlap goes to the smaller polygon, as
    16c's D2 labels it. A lake inside an unholed forest keeps its area under
    the default repair; the forest loses it."""

    FOREST = box(0, 0, 140, 50)
    LAKE = box(50, 15, 90, 35)

    def test_the_lake_keeps_its_area(self, fi: ModuleType, tmp_path: Path) -> None:
        feats = [coded("f", self.FOREST, "312"), coded("l", self.LAKE, "512")]
        fs = opened(fi, tmp_path, feats, repair_m=S)
        by_fid = {f.fid: f.polygon for f in fs.features}
        assert by_fid["l"].area == pytest.approx(self.LAKE.area, rel=1e-9)
        assert by_fid["f"].area == pytest.approx(self.FOREST.area - self.LAKE.area, rel=1e-9)
