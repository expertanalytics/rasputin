"""A tolerance that varies with distance to named lines: the binding and ``rasputin mesh``.

Increment 33 (``docs/increments/33-feature-tolerance.md``, sections 4.5, 5 and
9, tests 12, 13 and 14; the invariant-critical suites are C++ and
``test_tolerance_field.py``). ``--tolerance-near FILE N``, ``--tolerance-ramp
START END`` and ``--tolerance-near-crs CRS`` hold the mesh to ``N`` metres on
the lines and to ``--tolerance`` (F) from ``END`` metres out, linearly between.

- test 12: each refusal of section 5, in its words; the stderr line for lines
  out of reach (and then the ``--tolerance F`` mesh, G2); the record's four new
  fields; ``max_error_near_lines_m <= N``; none of them without the flag.
- test 13: on the resampled path (``--out-crs``) a line given in the DEM's CRS
  is transformed to the target grid's CRS, so the mesh is denser along the line
  where the mesh's own coordinates put it.
- test 14: more vertices within 100 m of the line than the ``--tolerance F``
  mesh has there (no count against the ``--tolerance N`` mesh: greedy
  insertion is not monotone in the tolerance).
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from pyproj import Transformer
from shapely.geometry import LineString, Point, Polygon, mapping
from shapely.geometry.base import BaseGeometry

import tin_engine.cli as cli
import tin_engine.edge_strip as edge_strip
import tin_engine.final_check as final_check
from cli_driver import SQUARE, UTM33, geojson, invoke, refused, rough_dem
from geotiff_fixtures import TIE_X, TIE_Y
from gpkg_fixtures import Layer, Row, write_gpkg
from test_cli_constraint_feet_all_paths import Spy
from test_cli_mesh_geographic import TARGET, domain_4674, geographic_dem, lonlat_ring
from test_cli_mesh_stats import mesh
from test_core_refine_points import STUB, X_MIN, Y_MAX, H, Phase1, scattered, surface, world
from tin_engine import _core
from tin_engine._core import ChainRole
from tin_engine.cli import DEFAULT_SNAP_SPACING, _constraint_arrays, _engine
from vtkread import VtkFile, read_vtk

__all__ = ["domain_4674", "geographic_dem"]  # fixtures, used by name

FIELDS = ("tolerance_near_m", "tolerance_ramp_m", "tolerance_lines", "max_error_near_lines_m")

bumpy = rough_dem(20)

#: A line across SQUARE from north to south, 48 m east of its west side.
LINE_X = TIE_X + 60.0
LINE = LineString([(LINE_X, TIE_Y + 10.0), (LINE_X, TIE_Y - 90.0)])


def lines_file(path: Path, geometries: list[BaseGeometry], crs: str | None = UTM33) -> Path:
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


@pytest.fixture
def box(tmp_path: Path) -> Path:
    return geojson(tmp_path / "box.geojson", SQUARE)


@pytest.fixture
def line(tmp_path: Path) -> Path:
    return lines_file(tmp_path / "line.geojson", [LINE])


def run(tmp_path: Path, *args: str, out: str = "x") -> tuple[VtkFile, dict[str, Any], str]:
    """(mesh, record, output) of one ``rasputin mesh`` run."""
    vtk, rec = tmp_path / f"{out}.vtk", tmp_path / f"{out}.json"
    result = mesh(*args, "--out", str(vtk), "--record", str(rec))
    return read_vtk(vtk.read_bytes()), json.loads(rec.read_text(encoding="ascii")), result.output


def near_line(points: np.ndarray, a: tuple[float, float], b: tuple[float, float], d: float) -> int:
    """How many of ``points`` (``(n, >=2)``) lie within ``d`` of segment ``a b``."""
    p = np.asarray(points, dtype=np.float64)[:, :2]
    a_, u = np.asarray(a), np.asarray(b) - np.asarray(a)
    t = np.clip(((p - a_) @ u) / (u @ u), 0.0, 1.0)
    return int((np.hypot(*(p - (a_ + t[:, None] * u)).T) <= d).sum())


# ---------------------------------------------------------------- binding


def start(n: int = 25) -> tuple[Any, ...]:
    """Phase 1's input as ``Phase1`` builds it: the DEM (kept alive), its view,
    the start mesh from the grid's perimeter ring, and its constraint arrays."""
    r, c = np.indices((n, n), dtype=np.float64)
    dem = surface(c, r).astype(np.float32)
    ring = (
        [(0, j) for j in range(0, n - 1, 4)]
        + [(i, n - 1) for i in range(0, n - 1, 4)]
        + [(n - 1, j) for j in range(n - 1, 0, -4)]
        + [(i, 0) for i in range(n - 1, 0, -4)]
    )[::-1]
    xy = world(np.array([c for _, c in ring], float), np.array([r for r, _ in ring], float))
    built = _engine(xy, [([*range(len(ring))], ChainRole.Outer, 0)], True, DEFAULT_SNAP_SPACING)
    assert built.mesh is not None and built.noded is not None, built.message
    edges, masks = _constraint_arrays(built.mesh, built.noded)
    view = _core.raster_view(dem, x_min=X_MIN, y_max=Y_MAX, delta_x=float(H), delta_y=float(H))
    return dem, view, built.mesh, edges, masks


#: Across the 25 x 25 grid (x 500 000 to 500 720, y 6 999 280 to 7 000 000).
CROSSING = np.array([[X_MIN - 10.0, Y_MAX - 300.0, X_MIN + 730.0, Y_MAX - 500.0]])


def same_mesh(a: Any, b: Any) -> None:
    for name in ("vertices", "triangles", "z", "valid", "edges", "masks"):
        assert_array_equal(np.asarray(getattr(a, name)), np.asarray(getattr(b, name)), name)


class TestBinding:
    def test_a_field_with_n_equal_to_f_is_todays_refine(self) -> None:
        _dem, view, start_mesh, edges, masks = start()
        f = _core.LineTolerance(view, CROSSING, 4.0, 4.0, 0.0, 300.0, 1.0)
        plain = _core.refine(view, start_mesh, edges, masks, tolerance=4.0)
        with_field = _core.refine(view, start_mesh, edges, masks, tolerance=4.0, field=f)
        assert with_field.ok(), with_field.message
        same_mesh(plain, with_field)

    def test_the_ramp_holds_the_triangles_on_the_line_to_n(self) -> None:
        _dem, view, start_mesh, edges, masks = start()
        f = _core.LineTolerance(view, CROSSING, 0.5, 4.0, 0.0, 300.0, 1.0)
        out = _core.refine(view, start_mesh, edges, masks, tolerance=4.0, field=f)
        assert out.ok(), out.message
        assert isinstance(out.max_error_near, float)
        assert out.max_error_near <= 0.5
        assert out.max_error <= 4.0

    def test_refine_strip_and_refine_points_take_the_field(self) -> None:
        phase1 = Phase1(n=25, tolerance=4.0)
        args = phase1.args()
        dem = surface(*np.indices((25, 25), dtype=np.float64)[::-1]).astype(np.float32)
        view = _core.raster_view(dem, x_min=X_MIN, y_max=Y_MAX, delta_x=float(H), delta_y=float(H))
        f = _core.LineTolerance(view, CROSSING, 4.0, 4.0, 0.0, 300.0, 1.0)
        strip = _core.constraint_check_points(view, args[0], args[4])
        same_mesh(
            _core.refine_strip(view, strip, *args, tolerance=4.0),
            _core.refine_strip(view, strip, *args, tolerance=4.0, field=f),
        )
        xy, z = scattered(25, 2, seed=33)
        cp = phase1.store(_core.CheckPoints, xy, z)
        same_mesh(
            _core.refine_points(cp, *args, tolerance=4.0),
            _core.refine_points(cp, *args, tolerance=4.0, field=f),
        )

    @pytest.mark.parametrize(
        ("ramp", "word"),
        [
            ((5.0, 4.0, 0.0, 300.0, 1.0), "near"),
            ((-1.0, 4.0, 0.0, 300.0, 1.0), "near"),
            ((1.0, 4.0, 400.0, 300.0, 1.0), "start"),
            ((1.0, 4.0, 0.0, 300.0, -1.0), "margin"),
            ((1.0, float("nan"), 0.0, 300.0, 1.0), "far"),
        ],
    )
    def test_a_bad_ramp_is_a_value_error_naming_the_bound(
        self, ramp: tuple[float, ...], word: str
    ) -> None:
        _dem, view, *_ = start()
        with pytest.raises(ValueError, match=word):
            _core.LineTolerance(view, CROSSING, *ramp)

    def test_a_non_finite_coordinate_is_a_value_error(self) -> None:
        _dem, view, *_ = start()
        bad = CROSSING.copy()
        bad[0, 2] = np.inf
        with pytest.raises(ValueError, match="coordinate"):
            _core.LineTolerance(view, bad, 1.0, 4.0, 0.0, 300.0, 1.0)


class TestOtherGeometry:
    """Section 9.1, B: a field measures in the lattice frame it was built on,
    so a call on another raster geometry is refused, not silently wrong."""

    WORDS = "the tolerance field was built on another raster geometry"

    def test_refine_on_another_view(self) -> None:
        dem, view, start_mesh, edges, masks = start()
        f = _core.LineTolerance(view, CROSSING, 0.5, 4.0, 0.0, 300.0, 1.0)
        # One cell further west and one column wider, so the start mesh is
        # still inside the other view's node rectangle.
        wider = np.ascontiguousarray(np.hstack([dem[:, :1], dem]))
        other = _core.raster_view(
            wider, x_min=X_MIN - H, y_max=Y_MAX, delta_x=float(H), delta_y=float(H)
        )
        with pytest.raises(ValueError, match=self.WORDS):
            _core.refine(other, start_mesh, edges, masks, tolerance=4.0, field=f)

    def test_refine_points_with_a_store_on_another_geometry(self) -> None:
        phase1 = Phase1(n=25, tolerance=4.0)
        _dem, view, *_ = start()
        f = _core.LineTolerance(view, CROSSING, 0.5, 4.0, 0.0, 300.0, 1.0)
        xy, z = scattered(25, 2, seed=33)
        cp = _core.CheckPoints(x_min=X_MIN - H, y_max=Y_MAX, spacing=H, rows=25, cols=26)
        cp.add(xy, z)
        cp.freeze()
        with pytest.raises(ValueError, match=self.WORDS):
            _core.refine_points(cp, *phase1.args(), tolerance=4.0, field=f)


class TestStub:
    """``_core.pyi`` declares what the binding adds (mypy reads the stub, not
    the module)."""

    @staticmethod
    def declared() -> dict[str, ast.stmt]:
        tree = ast.parse(STUB.read_text(encoding="utf-8"))
        return {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef | ast.ClassDef)}

    def test_line_tolerance_is_declared(self) -> None:
        assert isinstance(self.declared().get("LineTolerance"), ast.ClassDef)

    def test_refine_takes_the_field_last_by_keyword(self) -> None:
        fn = self.declared().get("refine")
        assert isinstance(fn, ast.FunctionDef)
        assert [a.arg for a in fn.args.kwonlyargs][-1] == "field"

    def test_the_outcome_carries_max_error_near(self) -> None:
        cls = self.declared().get("RefineOutcome")
        assert isinstance(cls, ast.ClassDef)
        assert "max_error_near" in {n.name for n in cls.body if isinstance(n, ast.FunctionDef)}


# ---------------------------------------------------------------- test 12: refusals


def ramp_args(bumpy: Path, box: Path, line: Path, *flags: str) -> tuple[str, ...]:
    return ("--dem", str(bumpy), "--domain", str(box), *flags)


class TestRefusals:
    def test_without_tolerance(self, tmp_path: Path, bumpy: Path, line: Path) -> None:
        """No --domain, so no other flag asks for --tolerance first."""
        out = refused(
            tmp_path,
            "mesh",
            *("--dem", str(bumpy), "--tolerance-near", str(line), "1"),
            *("--tolerance-ramp", "0", "40"),
            says=("--tolerance-near needs --tolerance",),
        )
        assert "--tolerance-near needs --tolerance-ramp" not in out

    def test_without_dem(self, tmp_path: Path, line: Path) -> None:
        refused(
            tmp_path,
            "mesh",
            *("catchment", "--flat", "--tolerance-near", str(line), "1"),
            *("--tolerance-ramp", "0", "40"),
            says=("applies only with --dem", "--tolerance-near"),
        )

    def test_a_crs_without_lines(self, tmp_path: Path, bumpy: Path, box: Path, line: Path) -> None:
        """Section 9.1: as --domain-crs without --domain."""
        args = ramp_args(bumpy, box, line, "--tolerance", "20", "--tolerance-near-crs", "EPSG:4326")
        refused(
            tmp_path,
            "mesh",
            *args,
            says=("applies only with --tolerance-near", "--tolerance-near-crs"),
        )

    def test_without_the_ramp(self, tmp_path: Path, bumpy: Path, box: Path, line: Path) -> None:
        args = ramp_args(bumpy, box, line, "--tolerance", "20", "--tolerance-near", str(line), "1")
        refused(tmp_path, "mesh", *args, says=("--tolerance-near needs --tolerance-ramp",))

    def test_a_ramp_without_lines(self, tmp_path: Path, bumpy: Path, box: Path, line: Path) -> None:
        args = ramp_args(bumpy, box, line, "--tolerance", "20", "--tolerance-ramp", "0", "40")
        refused(tmp_path, "mesh", *args, says=("--tolerance-ramp needs --tolerance-near",))

    def test_n_above_f(self, tmp_path: Path, bumpy: Path, box: Path, line: Path) -> None:
        args = ramp_args(bumpy, box, line, "--tolerance", "20", "--tolerance-near", str(line), "25")
        refused(
            tmp_path,
            "mesh",
            *args,
            "--tolerance-ramp",
            "0",
            "40",
            says=("--tolerance-near", "above", "25", "20"),
        )

    def test_start_above_end(self, tmp_path: Path, bumpy: Path, box: Path, line: Path) -> None:
        args = ramp_args(bumpy, box, line, "--tolerance", "20", "--tolerance-near", str(line), "1")
        refused(
            tmp_path,
            "mesh",
            *args,
            "--tolerance-ramp",
            "500",
            "100",
            says=("--tolerance-ramp", "above", "500", "100"),
        )

    @pytest.mark.parametrize(
        ("near", "ramp", "flag"),
        [
            ("-1", ("0", "40"), "--tolerance-near"),
            ("nan", ("0", "40"), "--tolerance-near"),
            ("inf", ("0", "40"), "--tolerance-near"),
            ("1", ("-5", "40"), "--tolerance-ramp"),
            ("1", ("0", "inf"), "--tolerance-ramp"),
            ("1", ("nan", "40"), "--tolerance-ramp"),
        ],
    )
    def test_a_negative_or_non_finite_number(
        self,
        tmp_path: Path,
        bumpy: Path,
        box: Path,
        line: Path,
        near: str,
        ramp: tuple[str, str],
        flag: str,
    ) -> None:
        args = ramp_args(bumpy, box, line, "--tolerance", "20", "--tolerance-near", str(line), near)
        refused(
            tmp_path,
            "mesh",
            *args,
            "--tolerance-ramp",
            *ramp,
            says=(flag, "must be finite and >= 0"),
        )

    @pytest.mark.parametrize(
        "geometry",
        [
            Polygon([(TIE_X, TIE_Y), (TIE_X + 10, TIE_Y), (TIE_X + 10, TIE_Y - 10)]),
            Point(TIE_X, TIE_Y),
        ],
        ids=["polygons", "points"],
    )
    def test_a_file_with_no_lines(
        self, tmp_path: Path, bumpy: Path, box: Path, geometry: BaseGeometry
    ) -> None:
        bad = lines_file(tmp_path / "bergen_line.geojson", [geometry])
        args = ramp_args(bumpy, box, bad, "--tolerance", "20", "--tolerance-near", str(bad), "1")
        refused(
            tmp_path,
            "mesh",
            *args,
            "--tolerance-ramp",
            "0",
            "40",
            says=("bergen_line.geojson has no lines; polygons and points are not used here",),
            squash=True,
        )


# ---------------------------------------------------------------- test 12: runs and record


class TestRecord:
    FLAGS = ("--tolerance", "20", "--tolerance-ramp", "0", "40")

    def test_the_four_fields_and_the_bound_near_the_lines(
        self, tmp_path: Path, bumpy: Path, box: Path, line: Path
    ) -> None:
        _, record, _ = run(
            tmp_path,
            *ramp_args(bumpy, box, line, *self.FLAGS, "--tolerance-near", str(line), "1"),
        )
        assert record["tolerance_m"] == 20.0  # still the far value
        assert type(record["tolerance_near_m"]) is float
        assert record["tolerance_near_m"] == 1.0
        assert record["tolerance_ramp_m"] == "0 to 40"
        lines = record["tolerance_lines"]
        assert lines.startswith("line.geojson: "), lines
        assert lines == "line.geojson: 1 segment", lines  # section 9.1, pin 14
        assert 0.0 <= record["max_error_near_lines_m"] <= 1.0

    def test_none_of_them_without_the_flag(self, tmp_path: Path, bumpy: Path, box: Path) -> None:
        _, record, _ = run(tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "20")
        for name in FIELDS:
            assert name not in record, name

    @pytest.mark.parametrize("suffix", ["geojson", "gpkg"])
    def test_lines_out_of_reach_say_so_and_give_the_far_mesh(
        self,
        tmp_path: Path,
        bumpy: Path,
        box: Path,
        monkeypatch: pytest.MonkeyPatch,
        suffix: str,
    ) -> None:
        """G2: no segment within reach is the ``--tolerance F`` mesh. In a
        GeoPackage the R-tree query finds no row at all, which is still zero
        segments, not a file with no lines (code review round 1, fix 1)."""
        far_line = LineString([(TIE_X + 10_000, TIE_Y), (TIE_X + 10_000, TIE_Y - 100)])
        if suffix == "gpkg":
            layer = Layer("lines", 25833, [Row(1, far_line)])
            far = write_gpkg(tmp_path / "far.gpkg", [layer])
        else:
            far = lines_file(tmp_path / "far.geojson", [far_line])
        args = ("--dem", str(bumpy), "--domain", str(box), "--tolerance", "20")
        ramp = ("--tolerance-near", str(far), "1", "--tolerance-ramp", "0", "3000")
        phase = Spy(monkeypatch, cli, "refine")
        strip = Spy(monkeypatch, edge_strip, "refine_strip")
        vtk, record, output = run(tmp_path, *args, *ramp, out="ramp")
        # No field is built for zero segments (section 9.1).
        assert [k.get("field") for k in phase.kwargs + strip.kwargs] == [None, None]
        plain_vtk, _, _ = run(tmp_path, *args, out="plain")
        # Section 9.1, "Zero segments kept": three of the four fields, no
        # max_error_near_lines_m, since no triangle is near a line.
        assert record["tolerance_lines"] == f"far.{suffix}: 0 segments"
        assert record["tolerance_near_m"] == 1.0
        assert record["tolerance_ramp_m"] == "0 to 3000"
        assert "max_error_near_lines_m" not in record
        text = " ".join(output.split())
        assert (
            "no tolerance lines within 3000 m of the domain; every triangle is held to 20 m" in text
        ), text
        assert_array_equal(vtk.points, plain_vtk.points)
        assert len(vtk.polygons) == len(plain_vtk.polygons)
        for a, b in zip(vtk.polygons, plain_vtk.polygons, strict=True):
            assert_array_equal(a, b)

    def test_the_field_reaches_refine_and_the_edge_strip(
        self, tmp_path: Path, bumpy: Path, box: Path, line: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Spy(monkeypatch, cli, "refine")
        strip = Spy(monkeypatch, edge_strip, "refine_strip")
        run(tmp_path, *ramp_args(bumpy, box, line, *self.FLAGS, "--tolerance-near", str(line), "1"))
        assert [k.get("field") is not None for k in phase.kwargs] == [True]
        assert [k.get("field") is not None for k in strip.kwargs] == [True]

    def test_no_field_without_the_flag(
        self, tmp_path: Path, bumpy: Path, box: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        phase = Spy(monkeypatch, cli, "refine")
        run(tmp_path, "--dem", str(bumpy), "--domain", str(box), "--tolerance", "20")
        assert [k.get("field") for k in phase.kwargs] == [None]


# ---------------------------------------------------------------- test 13: resampled path


class TestResampled:
    """A geographic DEM meshed onto ``TARGET`` (EPSG:31983); the line is given
    in the DEM's CRS (EPSG:4326, no ``crs`` member), across the domain."""

    def test_the_line_is_where_the_mesh_coordinates_put_it(
        self,
        tmp_path: Path,
        geographic_dem: Path,
        domain_4674: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        lonlat = lonlat_ring([(30.0, 0.0), (30.0, 59.0)])
        path = lines_file(tmp_path / "across.geojson", [LineString(lonlat)], crs=None)
        to_target = Transformer.from_crs("EPSG:4326", TARGET, always_xy=True)
        (ax, bx), (ay, by) = to_target.transform([p[0] for p in lonlat], [p[1] for p in lonlat])
        args = ("--dem", str(geographic_dem), "--domain", str(domain_4674), "--out-crs", TARGET)
        final = Spy(monkeypatch, final_check, "refine_points")
        vtk, record, output = run(
            tmp_path,
            *args,
            *("--tolerance", "5", "--tolerance-near", str(path), "0.5"),
            *("--tolerance-ramp", "0", "200"),
            out="ramp",
        )
        assert "no tolerance lines" not in output
        assert not record["tolerance_lines"].endswith(": 0 segments"), record["tolerance_lines"]
        assert [k.get("field") is not None for k in final.kwargs] == [True]
        plain, _, _ = run(tmp_path, *args, "--tolerance", "5", out="plain")
        near = near_line(vtk.points, (ax, ay), (bx, by), 100.0)
        assert near > near_line(plain.points, (ax, ay), (bx, by), 100.0)


# ---------------------------------------------------------------- test 14: the effect


def test_more_vertices_near_the_line_than_the_far_mesh(
    tmp_path: Path, bumpy: Path, box: Path, line: Path
) -> None:
    args = ("--dem", str(bumpy), "--domain", str(box), "--tolerance", "20")
    ramp = ("--tolerance-near", str(line), "1", "--tolerance-ramp", "0", "40")
    with_lines, _, _ = run(tmp_path, *args, *ramp, out="ramp")
    without, _, _ = run(tmp_path, *args, out="plain")
    a, b = (LINE_X, TIE_Y + 10.0), (LINE_X, TIE_Y - 90.0)
    assert near_line(with_lines.points, a, b, 100.0) > near_line(without.points, a, b, 100.0)


def test_the_flags_are_documented_in_plain_words() -> None:
    code, output = invoke("mesh", "--help")
    assert code == 0, output
    for flag in ("--tolerance-near", "--tolerance-ramp", "--tolerance-near-crs"):
        assert flag in output, flag
