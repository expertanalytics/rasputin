"""``tolerance_field.line_segments``: the lines a tolerance follows, ready for C++.

Increment 33 (``docs/increments/33-feature-tolerance.md``, sections 3, 4.5, 8
and 9, test 11; invariant-critical, mutants M8 and M9). ``line_segments(spec,
window, mesh_crs)`` reads the file through ``feature_input.read_source``,
transforms every LineString and MultiLineString to ``mesh_crs``
(``always_xy``), simplifies each by ``margin_m`` (Douglas-Peucker, no topology
kept), keeps the whole segments whose box meets ``window`` grown by ``end_m +
margin_m``, and returns them as a float64 ``(k, 4)`` array ``x0 y0 x1 y1`` with
the margin they were simplified by. A line the simplification empties comes
back as a zero-length segment at its first vertex. Polygons and points, and a
NaN or infinite input coordinate, are refused in words naming the file.

The guarantee checked (G7): for every point, its distance to the returned
segments less the returned margin is never above its distance to the original
line, transformed and not simplified. Distances here are this file's own
(point to segment, clamped projection), not the product's.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
import shapely
from pydantic import ValidationError
from pyproj import Transformer
from shapely.geometry import LineString, MultiLineString, Point, Polygon, mapping
from shapely.geometry.base import BaseGeometry

import tin_engine.tolerance_field as tf
from gpkg_fixtures import Layer, Row, write_gpkg

UTM33 = "EPSG:25833"
UTM33_URN = "urn:ogc:def:crs:EPSG::25833"
UTM32 = "EPSG:25832"
X0, Y0 = 500_000.0, 6_600_000.0
#: A 1 km window in UTM33.
WINDOW = (X0, Y0, X0 + 1000.0, Y0 + 1000.0)


def write(path: Path, geometries: list[BaseGeometry], crs: str | None = UTM33_URN) -> Path:
    """A FeatureCollection of ``geometries``; ``crs`` None writes no ``crs``
    member, so the file is RFC 7946's EPSG:4326."""
    doc: dict[str, Any] = {
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature", "geometry": mapping(g), "properties": {}} for g in geometries
        ],
    }
    if crs is not None:
        doc["crs"] = {"type": "name", "properties": {"name": crs}}
    path.write_text(json.dumps(doc))
    return path


def spec(
    path: Path,
    *,
    crs: str | None = None,
    near: float = 1.0,
    start: float = 0.0,
    end: float = 3000.0,
    margin: float = 1.0,
) -> Any:
    return tf.ToleranceLines(
        path=path, crs=crs, near_m=near, start_m=start, end_m=end, margin_m=margin
    )


def segments_of(
    s: Any, window: tuple[float, float, float, float], mesh_crs: str
) -> tuple[npt.NDArray[np.float64], float]:
    segs, margin = tf.line_segments(s, window, mesh_crs)
    assert isinstance(segs, np.ndarray)
    assert segs.dtype == np.float64
    assert segs.ndim == 2
    assert segs.shape[1] == 4
    return segs, margin


def as_set(segs: npt.NDArray[np.float64]) -> set[tuple[float, ...]]:
    return {tuple(round(v, 6) for v in row) for row in segs}


def point_segment(p: npt.NDArray[np.float64], segs: npt.NDArray[np.float64]) -> np.ndarray:
    """Each point's distance to the nearest of ``segs`` (``(k, 4)``), clamped
    projection; a zero-length segment is a point."""
    a, b = segs[:, 0:2], segs[:, 2:4]
    u = b - a
    len2 = np.einsum("ij,ij->i", u, u)
    w = p[:, None, :] - a[None, :, :]
    t = np.divide(np.einsum("pij,ij->pi", w, u), len2, out=np.zeros(w.shape[:2]), where=len2 > 0)
    t = np.clip(t, 0.0, 1.0)
    foot = a[None, :, :] + t[:, :, None] * u[None, :, :]
    return np.asarray(np.hypot(*(p[:, None, :] - foot).transpose(2, 0, 1)).min(axis=1))


def line_rows(line: LineString) -> npt.NDArray[np.float64]:
    xy = np.asarray(line.coords, dtype=np.float64)
    return np.hstack([xy[:-1], xy[1:]])


def zigzag(amplitude: float, period: float, length: float) -> LineString:
    """A line along x from (X0 + 100, Y0 + 500), ``amplitude`` metres either side."""
    n = int(length / (period / 2))
    xs = X0 + 100.0 + np.arange(n + 1) * period / 2
    ys = Y0 + 500.0 + amplitude * np.where(np.arange(n + 1) % 2 == 0, 1.0, -1.0)
    return LineString(np.column_stack([xs, ys]))


def near_points(line: LineString, band: float, count: int, seed: int) -> npt.NDArray[np.float64]:
    """``count`` points within ``band`` metres of ``line``'s box, seeded."""
    rng = np.random.default_rng(seed)
    x0, y0, x1, y1 = line.bounds
    return np.column_stack(
        [rng.uniform(x0 - band, x1 + band, count), rng.uniform(y0 - band, y1 + band, count)]
    )


def to_lonlat(line: LineString, crs: str) -> LineString:
    back = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    x, y = back.transform(*np.asarray(line.coords).T)
    return LineString(np.column_stack([x, y]))


def to_crs(line: LineString, crs: str) -> LineString:
    there = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = there.transform(*np.asarray(line.coords).T)
    return LineString(np.column_stack([x, y]))


# ---------------------------------------------------------------- the CRS


class TestTransform:
    def test_a_line_in_epsg_4326_comes_back_in_the_mesh_crs(self, tmp_path: Path) -> None:
        """No ``crs`` member: RFC 7946's EPSG:4326. Both ends are the PROJ
        transform's, to 1e-6 m (coordinates up to 7e6 m, checked there)."""
        lonlat = LineString([(8.00, 60.60), (8.01, 60.61)])
        path = write(tmp_path / "line.geojson", [lonlat], crs=None)
        expected = to_crs(lonlat, UTM32)
        x0, y0, x1, y1 = expected.bounds
        segs, margin = segments_of(spec(path), (x0 - 10, y0 - 10, x1 + 10, y1 + 10), UTM32)
        assert margin == 1.0
        assert segs.shape == (1, 4)
        np.testing.assert_allclose(segs[0], line_rows(expected)[0], rtol=0, atol=1e-6)

    def test_the_crs_given_is_used_when_the_file_says_none(self, tmp_path: Path) -> None:
        lonlat = LineString([(8.00, 60.60), (8.01, 60.61)])
        path = write(tmp_path / "line.geojson", [lonlat], crs=None)
        expected = to_crs(lonlat, UTM32)
        x0, y0, x1, y1 = expected.bounds
        given = spec(path, crs="EPSG:4326")
        segs, _ = segments_of(given, (x0 - 10, y0 - 10, x1 + 10, y1 + 10), UTM32)
        np.testing.assert_allclose(segs[0], line_rows(expected)[0], rtol=0, atol=1e-6)


# ---------------------------------------------------------------- the shapes


class TestShapes:
    def test_a_multilinestring_is_split_into_its_segments(self, tmp_path: Path) -> None:
        a = LineString([(X0 + 100, Y0 + 100), (X0 + 300, Y0 + 100)])
        b = LineString([(X0 + 100, Y0 + 500), (X0 + 300, Y0 + 500), (X0 + 300, Y0 + 800)])
        path = write(tmp_path / "multi.geojson", [MultiLineString([a, b])])
        segs, _ = segments_of(spec(path), WINDOW, UTM33)
        assert as_set(segs) == as_set(np.vstack([line_rows(a), line_rows(b)]))

    def test_a_linestring_gives_one_row_per_kept_segment(self, tmp_path: Path) -> None:
        line = LineString([(X0 + 10, Y0 + 10), (X0 + 400, Y0 + 20), (X0 + 420, Y0 + 600)])
        path = write(tmp_path / "line.geojson", [line])
        segs, _ = segments_of(spec(path), WINDOW, UTM33)
        assert as_set(segs) == as_set(line_rows(line))

    @pytest.mark.parametrize(
        "geometry",
        [
            Polygon([(X0, Y0), (X0 + 10, Y0), (X0 + 10, Y0 + 10)]),
            Point(X0 + 5, Y0 + 5),
        ],
        ids=["polygon", "point"],
    )
    def test_a_file_with_no_lines_is_refused_in_words(
        self, tmp_path: Path, geometry: BaseGeometry
    ) -> None:
        path = write(tmp_path / "bergen_line.geojson", [geometry])
        with pytest.raises(
            ValueError,
            match=r"bergen_line\.geojson has no lines; polygons and points are not used here",
        ):
            tf.line_segments(spec(path), WINDOW, UTM33)

    @pytest.mark.parametrize(
        "other",
        [Polygon([(X0, Y0), (X0 + 10, Y0), (X0 + 10, Y0 + 10)]), Point(X0 + 5, Y0 + 5)],
        ids=["polygon", "point"],
    )
    def test_a_polygon_or_point_beside_lines_is_refused(
        self, tmp_path: Path, other: BaseGeometry
    ) -> None:
        line = LineString([(X0 + 100, Y0 + 100), (X0 + 300, Y0 + 100)])
        path = write(tmp_path / "mixed.geojson", [line, other])
        with pytest.raises(
            ValueError, match=r"mixed\.geojson.*polygons and points are not used here"
        ):
            tf.line_segments(spec(path), WINDOW, UTM33)


# ---------------------------------------------------------------- the selection


class TestSelection:
    """``E = 3000``, margin 1: the window grows by 3001 m (M9)."""

    @staticmethod
    def vertical(x: float) -> LineString:
        return LineString([(x, Y0 + 100), (x, Y0 + 900)])

    def test_a_line_2999_m_outside_is_kept_and_one_at_3002_m_dropped(self, tmp_path: Path) -> None:
        east_in, east_out = self.vertical(X0 + 1000 + 2999), self.vertical(X0 + 1000 + 3002)
        west_in, west_out = self.vertical(X0 - 2999), self.vertical(X0 - 3002)
        path = write(tmp_path / "lines.geojson", [east_in, east_out, west_in, west_out])
        segs, _ = segments_of(spec(path), WINDOW, UTM33)
        assert as_set(segs) == as_set(np.vstack([line_rows(east_in), line_rows(west_in)]))

    def test_a_segment_exactly_at_the_reach_is_kept(self, tmp_path: Path) -> None:
        """Section 9.1: the box test is closed. 504 001 = 501 000 + 3 001, exact."""
        at_reach = self.vertical(X0 + 1000 + 3001)
        path = write(tmp_path / "edge.geojson", [at_reach])
        segs, _ = segments_of(spec(path), WINDOW, UTM33)
        assert as_set(segs) == as_set(line_rows(at_reach))

    def test_a_geopackage_with_every_line_out_of_reach_is_an_empty_array(
        self, tmp_path: Path
    ) -> None:
        """Code review round 1, fix 1: the R-tree query with the grown window
        finds no row, which is zero segments (section 5: stderr says so), not
        a file with no lines."""
        far = self.vertical(X0 + 1000 + 10_000)
        path = write_gpkg(tmp_path / "far.gpkg", [Layer("lines", 25833, [Row(1, far)])])
        segs, margin = segments_of(spec(path), WINDOW, UTM33)
        assert segs.shape == (0, 4)
        assert margin == 1.0

    def test_a_geopackage_with_a_polygon_in_reach_and_lines_out_of_reach_is_refused(
        self, tmp_path: Path
    ) -> None:
        """Code review round 3: the R-tree returns only the polygon, so the file
        must be read whole to see its lines, and is then a mixed file (section
        9.1, pin 10), not one with no lines."""
        polygon = Polygon([(X0 + 100, Y0 + 100), (X0 + 200, Y0 + 100), (X0 + 200, Y0 + 200)])
        far = self.vertical(X0 + 1000 + 10_000)
        layer = Layer("shapes", 25833, [Row(1, polygon), Row(2, far)])
        path = write_gpkg(tmp_path / "mixed_far.gpkg", [layer])
        with pytest.raises(
            ValueError,
            match=r"mixed_far\.gpkg holds lines and other shapes; "
            r"polygons and points are not used here",
        ):
            tf.line_segments(spec(path), WINDOW, UTM33)

    def test_no_line_within_reach_is_an_empty_array(self, tmp_path: Path) -> None:
        path = write(tmp_path / "far.geojson", [self.vertical(X0 + 1000 + 3002)])
        segs, margin = segments_of(spec(path), WINDOW, UTM33)
        assert segs.shape == (0, 4)
        assert margin == 1.0

    def test_a_segment_reaching_into_the_window_is_kept_whole(self, tmp_path: Path) -> None:
        """Kept whole, never cut at the window (section 4.6: two pieces beside
        a seam hold the same segments)."""
        long = LineString([(X0 + 500, Y0 + 500), (X0 + 20_000, Y0 + 500)])
        path = write(tmp_path / "long.geojson", [long])
        segs, _ = segments_of(spec(path), WINDOW, UTM33)
        assert as_set(segs) == as_set(line_rows(long))


# ---------------------------------------------------------------- non-finite input


class TestNonFinite:
    """Code review round 1, fix 3 (Ola away; the main session took "refuse"):
    ``shapely.simplify`` drops a NaN vertex silently, so a NaN or infinite
    input coordinate is refused before it, naming the file."""

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
    def test_a_non_finite_coordinate_is_refused(self, tmp_path: Path, bad: float) -> None:
        line = LineString(
            [(X0 + 500, Y0 + 500), (X0 + 600, Y0 + 600), (bad, bad), (X0 + 700, Y0 + 500)]
        )
        layer = Layer("lines", 25833, [Row(1, line)], rtree=False)
        path = write_gpkg(tmp_path / "broken_line.gpkg", [layer])
        with pytest.raises(ValueError, match=r"broken_line\.gpkg"):
            tf.line_segments(spec(path), WINDOW, UTM33)


# ---------------------------------------------------------------- the margin (G7)


class TestMargin:
    """M8: ``d(segments) - margin <= d(original)`` at 1 000 seeded points
    within 3 m of the line's box. Slack 1e-9 m (coordinates to 7e6 m)."""

    @pytest.mark.parametrize("given_in", ["utm33", "epsg4326"])
    def test_the_distance_less_the_margin_is_never_above_the_original(
        self, tmp_path: Path, given_in: str
    ) -> None:
        # 0.45 m either side every 2.5 m: Douglas-Peucker at 1 m flattens it to
        # one segment along the upper peaks (GEOS 3.14: 81 vertices to 2), so a
        # lower peak is 0.9 m from the simplified line.
        utm = zigzag(0.45, 5.0, 200.0)
        if given_in == "utm33":
            path, original = write(tmp_path / "zig.geojson", [utm]), utm
        else:
            lonlat = to_lonlat(utm, UTM33)
            path, original = (
                write(tmp_path / "zig.geojson", [lonlat], crs=None),
                to_crs(lonlat, UTM33),
            )
        segs, margin = segments_of(spec(path, margin=1.0), WINDOW, UTM33)
        assert margin == 1.0
        assert 1 <= len(segs) < len(original.coords) - 1, "the line was not simplified"
        points = near_points(original, 3.0, 1000, seed=33)
        here = point_segment(points, segs)
        there = point_segment(points, line_rows(original))
        excess = here - margin - there
        assert excess.max() <= 1e-9, f"{(excess > 1e-9).sum()} points over, worst {excess.max()} m"
        # Not vacuous: without the margin the bound breaks here.
        assert (here - there).max() > 0.5

    def test_the_margin_returned_is_the_one_simplified_by(self, tmp_path: Path) -> None:
        """0.8 m either side: Douglas-Peucker keeps every vertex at 1 m and
        flattens the line at 2 m (GEOS 3.14), so a line simplified by more
        than the margin it returns breaks the bound by about 0.6 m."""
        utm = zigzag(0.8, 5.0, 200.0)
        path = write(tmp_path / "zig.geojson", [utm])
        segs, margin = segments_of(spec(path, margin=1.0), WINDOW, UTM33)
        points = near_points(utm, 3.0, 1000, seed=35)
        excess = point_segment(points, segs) - margin - point_segment(points, line_rows(utm))
        assert excess.max() <= 1e-9, f"worst {excess.max()} m"

    def test_a_line_simplification_empties_is_a_zero_length_segment(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """GEOS 3.14 never empties a line, so the fallback is reached only with
        ``shapely.simplify`` replaced (code review round 1, suggestion)."""
        monkeypatch.setattr(shapely, "simplify", lambda *_a, **_k: LineString())
        line = LineString([(X0 + 500, Y0 + 500), (X0 + 600, Y0 + 500)])
        path = write(tmp_path / "line.geojson", [line])
        segs, _ = segments_of(spec(path), WINDOW, UTM33)
        np.testing.assert_array_equal(segs, [[X0 + 500, Y0 + 500, X0 + 500, Y0 + 500]])

    def test_a_closed_loop_smaller_than_the_margin_is_kept(self, tmp_path: Path) -> None:
        loop = LineString(
            [(X0 + 500, Y0 + 500), (X0 + 500.5, Y0 + 500), (X0 + 500.5, Y0 + 500.5),
             (X0 + 500, Y0 + 500.5), (X0 + 500, Y0 + 500)]
        )  # fmt: skip
        path = write(tmp_path / "loop.geojson", [loop])
        segs, margin = segments_of(spec(path, margin=1.0), WINDOW, UTM33)
        assert len(segs) >= 1, "a loop smaller than the margin came back as nothing"
        points = near_points(loop, 3.0, 1000, seed=34)
        excess = point_segment(points, segs) - margin - point_segment(points, line_rows(loop))
        assert excess.max() <= 1e-9


# ---------------------------------------------------------------- increment 34


class TestToleranceSlope:
    """Increment 34, section 4.5: ``--tolerance-slope N START END`` as data,
    beside ``ToleranceLines``. The bounds are the CLI's to check (section 5)."""

    def test_its_fields(self) -> None:
        s = tf.ToleranceSlope(near_m=2.0, start_deg=25.0, end_deg=35.0)  # type: ignore[attr-defined]
        assert (s.near_m, s.start_deg, s.end_deg) == (2.0, 25.0, 35.0)
        assert set(type(s).model_fields) == {"near_m", "start_deg", "end_deg"}

    def test_it_is_frozen(self) -> None:
        s = tf.ToleranceSlope(near_m=2.0, start_deg=30.0, end_deg=30.0)  # type: ignore[attr-defined]
        with pytest.raises(ValidationError, match="frozen"):
            s.near_m = 1.0
