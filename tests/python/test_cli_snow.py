"""`rasputin snow`: the job from flags, the run, and an EnSight case ParaView opens.

`docs/increments/46-hbv-snow.md`, *4.6 The animation format* (`snow.case`,
`snow.geo`, one file per value per frame, `slope_deg`, `aspect_deg`,
`radiation_anomaly_mean` and `melt_out_day` without a time set, `snow.json`
with the job and the origin), *4.8 The Python surface* (the flags; `--out` a
new or empty directory), *4.9 Degeneracy and refusal policy* (a geographic
mesh refused; a zero-area triangle counted on stderr) and section 7's E1
(marked slow): the synthetic spring on the 2 m mesh against section 6's
table, rows 1-3.

`TestEndToEnd` is the demonstration: `synthetic-terrain`, `mesh`, `snow`, and
the case read back by the suite's own reader, by VTK's EnSight reader when the
`vtk` extra is installed, and by ParaView's `pvpython` when ParaView is
installed (`RASPUTIN_PVPYTHON`, else ParaView 6.1.1's default macOS path);
each of the last two skips when absent.

Every run here needs the sun (`_core.sun_positions`, PR 2).

Interface:

- The output folder holds exactly `snow.case`, `snow.geo`, `snow.json`,
  `swe<NNNN>.var` and `sca<NNNN>.var` for frames 0.., and `slope_deg.var`,
  `aspect_deg.var`, `radiation_anomaly_mean.var`, `melt_out_day.var`.
- `snow.json` is `{"job": <SnowJob dumped>, "origin": [x0, y0]}`, the
  origin the mesh's minimum corner; the geometry stores x - x0, y - y0.
- `melt_out_day` is the time in days of the first frame (frame 0 included)
  whose SCA is below 0.5, -1 if none.
- `--sun-correction` takes `none`, `melt-factor`, `threshold`; `--start` an
  ISO 8601 time with an offset. Refusals exit 2 and say why.
- The degenerate count is reported as "1 degenerate triangle". A triangle
  is degenerate by the receiver frame's own zero-area test (its ground edges
  span less than 1e-9 of the product of their lengths), so a near-zero-area
  sliver is counted and given no sun correction (dI exactly 0), never
  refused. An error naming a triangle gives its index in the mesh.
- `slope_deg.var` and `aspect_deg.var` are the receiver frame's per-triangle
  true east-north-up normals (`receivers_from_mesh`), never a grid normal
  turned by one convergence for the whole mesh (section 10.5).
- E1 runs only with `RASPUTIN_SLOW=1` (several minutes).

Cast shadows (section 10.5, 10.8's L5 and E2, 10.9's Bygdin run):

- `--dem PATH` (repeatable, as `mesh --dem`) and `--shadows on|off`, default
  `on` with `--dem` and `off` without. `snow.json`'s job records `shadows`: the
  default `RadiationRequest` dumped, or null.
- Refused before any file is written, in plain words: `--shadows on` without
  `--dem` ("shadows need the DEM the mesh was made from: give --dem"); a DEM in
  another CRS than the mesh (naming both, as `EPSG:25832` and `EPSG:25833`); a
  DEM that does not cover the mesh ("does not cover the mesh").
- E2 (`RASPUTIN_SLOW=1`): on the 2 m mesh, 2027-04-15, the plan-area-weighted
  domain mean of the day's dI falls by 0.5 to 2.5 W/m^2 with shadows on.
- The Bygdin run (`RASPUTIN_SLOW=1`, and DTM10 present): looks for
  `DTM10_UTM33_20260925` under `$RASPUTIN_DATA`, else under a `rasputin_data`
  directory beside any parent of this file (the main checkout's
  `../rasputin_data`, also from a worktree); skips when absent.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from cli_helpers import USAGE, invoke, ran, squashed
from ensightread import Animation
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.io.geotiff_write import encode_geotiff
from tin_engine.io.vtk_legacy import write_vtk
from tin_engine.radiation.request import RadiationRequest
from vtkread import read_vtk

X0, Y0, Z = 500000.0, 6800000.0, 300.0
TAN30 = math.tan(math.radians(30.0))
SIDE = 20.0
STATIC = ("slope_deg", "aspect_deg", "radiation_anomaly_mean", "melt_out_day")
SPRING_10 = ("--start", "2027-04-10T00:00:00+00:00")
PVPYTHON = Path(
    os.environ.get("RASPUTIN_PVPYTHON", "/Applications/ParaView-6.1.1.app/Contents/bin/pvpython")
)


def facets(
    extra: list[list[float]] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """A 30-degree south facet, then a 30-degree north one, centroid z 300 m.

    `extra`, three corners, is added as a third triangle.
    """
    rows: list[list[float]] = []
    for dx, sign in ((0.0, 1.0), (100.0, -1.0)):
        for x, y in ((0.0, 0.0), (SIDE, 0.0), (0.0, SIDE)):
            rows.append([X0 + dx + x, Y0 + y, Z + sign * TAN30 * (y - SIDE / 3)])
    tris = [[0, 1, 2], [3, 4, 5]]
    if extra:
        rows += extra
        tris.append([6, 7, 8])
    return np.array(rows), np.array(tris)


def mesh_file(
    path: Path,
    *,
    crs: str = "EPSG:25833",
    binary: bool = True,
    extra: list[list[float]] | None = None,
    arrays: tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]] | None = None,
) -> Path:
    """`facets(extra)` as a binary mesh in `crs`, or `arrays` (vertices,
    triangles) when given."""
    vertices, triangles = facets(extra) if arrays is None else arrays
    path.write_bytes(
        write_vtk(
            vertices,
            triangles=triangles,
            edges=np.empty((0, 2), dtype=np.int64),
            edge_masks=np.empty(0, dtype=np.uint32),
            vocabulary=DEFAULT_VOCABULARY,
            fields=(("crs", crs),),
            binary=binary,
        )
    )
    return path


def snow(mesh: Path, out: Path, *args: str) -> str:
    return ran("snow", "--mesh", str(mesh), "--out", str(out), *args)


BINS = "bins must be between 1 and 64, not 65"
SIDEWAYS = "'sideways' is not one of 'none', 'melt-factor', 'threshold'"


def refused(mesh: Path, out: Path, *args: str, says: str) -> str:
    before = sorted(out.iterdir()) if out.is_dir() else None
    code, output = invoke("snow", "--mesh", str(mesh), "--out", str(out), *args)
    assert code == USAGE, output
    assert "Traceback" not in output and "No such option" not in output, output
    # Squashed: Rich wraps a long path or a hyphenated word anywhere in its panel.
    assert squashed(says) in squashed(output), f"{says!r} not in {output!r}"
    if before is None:
        assert not out.exists() or not any(out.iterdir())
    else:
        assert sorted(out.iterdir()) == before
    return output


@pytest.fixture
def mesh(tmp_path: Path) -> Path:
    return mesh_file(tmp_path / "facets.vtk")


@pytest.fixture
def two_days(mesh: Path, tmp_path: Path) -> Path:
    out = tmp_path / "anim"
    snow(mesh, out, "--sun-correction", "melt-factor", *SPRING_10, "--days", "2")
    return out


class TestOutput:
    def test_the_folder_holds_the_case_and_nothing_else(self, two_days: Path) -> None:
        frames = [f"{v}{q:04d}.var" for v in ("swe", "sca") for q in range(3)]
        expected = {"snow.case", "snow.geo", "snow.json", *frames, *(f"{s}.var" for s in STATIC)}
        assert {p.name for p in two_days.iterdir()} == expected

    def test_the_case_reads_back_with_three_daily_frames(self, two_days: Path) -> None:
        anim = Animation.read(two_days / "snow.case")
        assert anim.case.times == (0.0, 24.0, 48.0)
        assert set(anim.case.timed) == {"swe", "sca"}
        assert set(anim.case.static) == set(STATIC)
        assert_array_equal(anim.geometry.triangles, [[0, 1, 2], [3, 4, 5]])

    def test_the_geometry_is_stored_from_the_minimum_corner(self, two_days: Path) -> None:
        vertices, _ = facets()
        origin = (vertices[:, 0].min(), vertices[:, 1].min())
        anim = Animation.read(two_days / "snow.case")
        assert_array_equal(
            anim.geometry.points[:, 0], (vertices[:, 0] - origin[0]).astype(np.float32)
        )
        assert_array_equal(
            anim.geometry.points[:, 1], (vertices[:, 1] - origin[1]).astype(np.float32)
        )
        record = json.loads((two_days / "snow.json").read_text(encoding="utf-8"))
        assert tuple(record["origin"]) == origin

    def test_snow_json_records_the_job_from_the_flags(self, two_days: Path) -> None:
        from tin_engine.snow.params import SnowJob

        record = json.loads((two_days / "snow.json").read_text(encoding="utf-8"))
        job = SnowJob.model_validate(record["job"])
        assert job.correction.mode == "melt_factor"
        assert job.spring.days == 2
        assert job.spring.start.isoformat() == "2027-04-10T00:00:00+00:00"
        assert job.initial_swe_mm == 600.0
        assert len(job.distribution.factors) == 10

    def test_the_frames_start_white_and_melt(self, two_days: Path) -> None:
        anim = Animation.read(two_days / "snow.case")
        assert_array_equal(anim.frame(two_days, "swe", 0), [600.0, 600.0])
        assert_array_equal(anim.frame(two_days, "sca", 0), [1.0, 1.0])
        last = anim.frame(two_days, "swe", 2)
        assert (last <= 600.0).all() and (last > 0.0).all()

    def test_slope_aspect_and_the_sun_per_facet(self, two_days: Path) -> None:
        anim = Animation.read(two_days / "snow.case")
        # True slope and aspect: UTM's scale 0.9996 makes a grid 30-degree
        # facet 29.99 degrees on the ground, and the north facet, 110 m east
        # of the central meridian, faces 0.002 degrees east of true north.
        assert_allclose(anim.static(two_days, "slope_deg"), [30.0, 30.0], atol=0.02)
        assert_allclose(anim.static(two_days, "aspect_deg"), [180.0, 0.0], atol=0.01)
        sun = anim.static(two_days, "radiation_anomaly_mean")
        assert sun[0] > 20.0 and sun[1] < -20.0
        assert_array_equal(anim.static(two_days, "melt_out_day"), [-1.0, -1.0])


class TestMeltOutDay:
    def test_it_is_the_first_frame_below_half_cover(self, mesh: Path, tmp_path: Path) -> None:
        out = tmp_path / "thin"
        snow(
            mesh,
            out,
            "--sun-correction",
            "none",
            *SPRING_10,
            "--days",
            "10",
            "--initial-swe",
            "5",
            "--snow-cv",
            "0",
        )
        anim = Animation.read(out / "snow.case")
        sca = np.array([anim.frame(out, "sca", q) for q in range(len(anim.case.times))])
        below = sca < 0.5
        assert below[-1].all(), "5 mm should melt out within ten days"
        expected = below.argmax(axis=0).astype(np.float32)  # daily frames: index = day
        assert_array_equal(anim.static(out, "melt_out_day"), expected)

    def test_no_snow_melts_out_on_day_0(self, mesh: Path, tmp_path: Path) -> None:
        out = tmp_path / "bare"
        snow(mesh, out, "--sun-correction", "none", *SPRING_10, "--days", "1", "--initial-swe", "0")
        assert_array_equal(
            Animation.read(out / "snow.case").static(out, "melt_out_day"), [0.0, 0.0]
        )


class TestFlags:
    def test_frame_hours_and_step_minutes_set_the_time_axis(
        self, mesh: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "six"
        snow(
            mesh,
            out,
            "--sun-correction",
            "threshold",
            *SPRING_10,
            "--days",
            "1",
            "--step-minutes",
            "30",
            "--frame-hours",
            "6",
            "--threads",
            "2",
        )
        assert Animation.read(out / "snow.case").case.times == (0.0, 6.0, 12.0, 18.0, 24.0)

    def test_bins_and_cv_set_the_distribution(self, mesh: Path, tmp_path: Path) -> None:
        from tin_engine.snow.params import SnowDistribution, SnowJob

        out = tmp_path / "bins"
        snow(
            mesh,
            out,
            "--sun-correction",
            "none",
            *SPRING_10,
            "--days",
            "1",
            "--bins",
            "4",
            "--snow-cv",
            "0.2",
        )
        job = SnowJob.model_validate(json.loads((out / "snow.json").read_text("utf-8"))["job"])
        assert job.distribution == SnowDistribution.lognormal(0.2, 4)

    def test_an_empty_existing_folder_is_accepted(self, mesh: Path, tmp_path: Path) -> None:
        out = tmp_path / "empty"
        out.mkdir()
        snow(mesh, out, "--sun-correction", "none", *SPRING_10, "--days", "1")
        assert (out / "snow.case").exists()

    def test_a_degenerate_triangle_is_counted_and_written(self, tmp_path: Path) -> None:
        line = [[X0 + 300.0, Y0, Z], [X0 + 310.0, Y0 + 10.0, Z], [X0 + 320.0, Y0 + 20.0, Z]]
        mesh = mesh_file(tmp_path / "degenerate.vtk", extra=line)
        out = tmp_path / "anim"
        output = snow(mesh, out, "--sun-correction", "melt-factor", *SPRING_10, "--days", "1")
        assert "1 degenerate triangle" in output
        anim = Animation.read(out / "snow.case")
        assert anim.static(out, "radiation_anomaly_mean")[2] == 0.0
        assert np.isfinite(anim.frame(out, "swe", 1)).all()


def true_facet(
    x: float, y: float, slope_deg: float, aspect_deg: float, *, side: float = SIDE
) -> npt.NDArray[np.float64]:
    """Three corners at (x, y) in EPSG:25833 on the plane of that true slope and
    true aspect, built through J = (1/m) R(gamma) at the centroid (pyproj)."""
    from pyproj import Proj

    proj = Proj("EPSG:25833")
    lon, lat = proj(x + side / 3, y + side / 3, inverse=True)
    f = proj.get_factors(lon, lat)
    g, m = math.radians(f.meridian_convergence), f.meridional_scale
    j = np.array([[math.cos(g), math.sin(g)], [-math.sin(g), math.cos(g)]]) / m
    s, a = math.radians(slope_deg), math.radians(aspect_deg)
    grad = -math.tan(s) * j.T @ np.array([math.sin(a), math.cos(a)])
    plan = ((0.0, 0.0), (side, 0.0), (0.0, side))
    return np.array(
        [[x + dx, y + dy, Z + grad @ (dx - side / 3, dy - side / 3)] for dx, dy in plan]
    )


#: Bygdin's centre in EPSG:25833 (convergence -5.724 degrees), and a point on
#: the central meridian (convergence 0): no one turn serves both.
BYGDIN_XY = (152147.55, 6821319.44)
#: A sliver near (X0, Y0): edges (10, 10) and (20, 20 + 1e-8) m, a cross
#: product of 1e-7 m^2 (not 0) and a sine of 2.5e-10, under the frame's 1e-9.
SLIVER = [[X0 + 300.0, Y0, Z], [X0 + 310.0, Y0 + 10.0, Z], [X0 + 320.0, Y0 + 20.0 + 1e-8, Z]]
#: Float32 in the .var files: a step of 1.5e-5 degrees near 180.
VAR_DEG = 1e-4


class TestOneSourceOfNormals:
    """Section 10.5: slope and aspect from the receiver frame's normals."""

    def test_slope_and_aspect_files_are_the_frames_true_normals(self, tmp_path: Path) -> None:
        from tin_engine.radiation.frame import receivers_from_mesh

        vertices = np.concatenate(
            [true_facet(*BYGDIN_XY, 30.0, 180.0), true_facet(X0, Y0, 20.0, 90.0)]
        )
        triangles = np.array([[0, 1, 2], [3, 4, 5]])
        mesh = mesh_file(tmp_path / "wide.vtk", arrays=(vertices, triangles))
        out = tmp_path / "anim"
        snow(mesh, out, "--sun-correction", "none", *SPRING_10, "--days", "1")
        anim = Animation.read(out / "snow.case")
        slope, aspect = anim.static(out, "slope_deg"), anim.static(out, "aspect_deg")
        assert_allclose(slope, [30.0, 20.0], rtol=0, atol=VAR_DEG)
        assert_allclose(aspect, [180.0, 90.0], rtol=0, atol=VAR_DEG)
        e, n, u = receivers_from_mesh(vertices, triangles, "EPSG:25833").normal.T
        assert_array_equal(slope, np.degrees(np.arctan2(np.hypot(e, n), u)).astype(np.float32))
        assert_array_equal(aspect, (np.degrees(np.arctan2(e, n)) % 360.0).astype(np.float32))


class TestSlivers:
    """Section 4.9: a near-zero-area triangle is counted, not refused."""

    @pytest.mark.parametrize("shadows", ["off", "on"])
    def test_a_sliver_runs_through_and_gets_no_sun_correction(
        self, tmp_path: Path, shadows: str
    ) -> None:
        mesh = mesh_file(tmp_path / "sliver.vtk", extra=SLIVER)
        dem = ("--dem", str(dem_file(tmp_path / "dem.tif"))) if shadows == "on" else ()
        out = tmp_path / "anim"
        output = snow(mesh, out, *dem, "--shadows", shadows, "--sun-correction", "melt-factor",
                      *SPRING_10, "--days", "1")  # fmt: skip
        assert "1 degenerate triangle" in output
        anim = Animation.read(out / "snow.case")
        mean = anim.static(out, "radiation_anomaly_mean")
        assert mean[2] == 0.0
        assert mean[0] > 20.0 and mean[1] < -20.0
        assert anim.static(out, "slope_deg")[2] == 0.0
        assert np.isfinite(anim.frame(out, "swe", 1)).all()

    def test_an_error_names_the_triangle_by_its_index_in_the_mesh(self, tmp_path: Path) -> None:
        # Triangle 0 has no area and is left out of the frame; triangle 2 lies
        # at 89.97 N (northing 9 995 000), past the frame's 89.9: it is
        # triangle 2 of the mesh, not 1 of the triangles kept.
        line = [[X0, Y0, Z], [X0 + 10.0, Y0 + 10.0, Z], [X0 + 20.0, Y0 + 20.0, Z]]
        fine = [[X0 + 100.0, Y0, Z], [X0 + 120.0, Y0, Z], [X0 + 100.0, Y0 + 20.0, Z]]
        pole = [[X0, 9_995_000.0, Z], [X0 + 20.0, 9_995_000.0, Z], [X0, 9_995_020.0, Z]]
        arrays = (np.array(line + fine + pole), np.arange(9).reshape(3, 3))
        mesh = mesh_file(tmp_path / "pole.vtk", arrays=arrays)
        refused(mesh, tmp_path / "out", "--sun-correction", "none", *SPRING_10,
                "--days", "1", says="triangle 2 has its centroid beyond 89.9")  # fmt: skip


class TestRefusals:
    def test_a_folder_with_files_in_it(self, mesh: Path, tmp_path: Path) -> None:
        out = tmp_path / "used"
        out.mkdir()
        (out / "keep.txt").write_text("x", encoding="ascii")
        refused(mesh, out, "--sun-correction", "none", says="used must be a new or empty directory")

    def test_a_mesh_in_degrees(self, tmp_path: Path) -> None:
        mesh = mesh_file(tmp_path / "geo.vtk", crs="EPSG:4326")
        refused(
            mesh,
            tmp_path / "out",
            "--sun-correction",
            "none",
            says="snow needs a mesh in a projected CRS in metres; this one is in degrees",
        )

    @pytest.mark.parametrize(
        ("crs", "unit"), [("EPSG:2263", "US survey foot"), ("EPSG:3359", "foot")]
    )
    def test_a_mesh_in_feet(self, tmp_path: Path, crs: str, unit: str) -> None:
        mesh = mesh_file(tmp_path / "feet.vtk", crs=crs)
        output = refused(
            mesh,
            tmp_path / "out",
            "--sun-correction",
            "none",
            "--days",
            "1",
            says="snow needs a mesh in a projected CRS in metres",
        )
        assert squashed(unit) in squashed(output), output

    def test_an_ascii_mesh(self, tmp_path: Path) -> None:
        mesh = mesh_file(tmp_path / "ascii.vtk", binary=False)
        refused(
            mesh,
            tmp_path / "out",
            "--sun-correction",
            "none",
            says="write the mesh without --ascii",
        )

    def test_a_start_without_a_time_zone(self, mesh: Path, tmp_path: Path) -> None:
        refused(
            mesh,
            tmp_path / "out",
            "--sun-correction",
            "none",
            "--start",
            "2027-04-10T00:00:00",
            says="2027-04-10T00:00:00 has no time zone; give it an offset or Z",
        )

    def test_65_bins(self, mesh: Path, tmp_path: Path) -> None:
        refused(mesh, tmp_path / "out", "--sun-correction", "none", "--bins", "65", says=BINS)

    def test_an_unknown_correction(self, mesh: Path, tmp_path: Path) -> None:
        refused(mesh, tmp_path / "out", "--sun-correction", "sideways", says=SIDEWAYS)

    def test_frames_between_steps(self, mesh: Path, tmp_path: Path) -> None:
        refused(
            mesh,
            tmp_path / "out",
            "--sun-correction",
            "none",
            "--step-minutes",
            "90",
            "--frame-hours",
            "2",
            says="a frame every 2 h is not a whole number of 90-minute steps",
        )

    def test_a_spring_that_is_not_whole_frames(self, mesh: Path, tmp_path: Path) -> None:
        refused(
            mesh,
            tmp_path / "out",
            "--sun-correction",
            "none",
            "--days",
            "91",
            "--frame-hours",
            "48",
            says="91 days is not a whole number of 48-hour frames",
        )

    def test_more_frames_than_four_digits_can_number(self, mesh: Path, tmp_path: Path) -> None:
        refused(
            mesh,
            tmp_path / "out",
            "--sun-correction",
            "none",
            "--days",
            "417",
            "--frame-hours",
            "1",
            says="417 days with a frame every 1 h is 10,009 frames, "
            "more than the 10,000 the animation can number",
        )


def reference_figures(vtk_path: Path, out: Path) -> dict[str, Any]:
    """Section 6's figures from a run's output: inclined-area weights, classes S and N."""
    mesh = read_vtk(vtk_path.read_bytes())
    pts = mesh.points
    tri = np.array(mesh.polygons, dtype=np.int64)
    a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    nrm = np.cross(b - a, c - a)
    nrm[nrm[:, 2] < 0] *= -1
    area3 = 0.5 * np.linalg.norm(nrm, axis=1)
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    zc = (a[:, 2] + b[:, 2] + c[:, 2]) / 3
    slope = np.degrees(np.arccos(np.clip(nrm[:, 2], -1, 1)))
    aspect = np.degrees(np.arctan2(nrm[:, 0], nrm[:, 1])) % 360
    band = (zc >= 900) & (zc <= 1100) & (slope >= 25) & (slope <= 35)
    south = band & (np.abs(aspect - 180) <= 22.5)
    north = band & ((aspect <= 22.5) | (aspect >= 337.5))
    anim = Animation.read(out / "snow.case")
    w = area3 / area3.sum()
    sca = {d: anim.frame(out, "sca", d).astype(np.float64) for d in (30, 60, 90)}
    swe = {d: anim.frame(out, "swe", d).astype(np.float64) for d in (30, 60, 90)}
    mo = anim.static(out, "melt_out_day")
    return {
        "sca": [float(w @ sca[d]) for d in (30, 60, 90)],
        "swe": [float(w @ swe[d]) for d in (30, 60, 90)],
        "south": float(np.median(mo[south])),
        "north": float(np.median(mo[north])),
    }


#: Section 6's table, rows 1-3 (2 m mesh, N = 10, NOAA's sun in the probe).
TABLE = {
    "none": ([1.0, 0.6298, 0.0829], [565.29, 203.58, 23.56], 83.0, 82.0),
    "melt-factor": ([1.0, 0.6559, 0.1152], [565.14, 226.07, 36.05], 75.5, -1.0),
    "threshold": ([1.0, 0.6337, 0.0932], [558.13, 206.93, 28.22], 75.0, -1.0),
}


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory: pytest.TempPathFactory) -> Path:
    folder = tmp_path_factory.mktemp("synthetic")
    ran("synthetic-terrain", "--out", str(folder / "synth.tif"))
    return folder


@pytest.fixture(scope="module")
def mesh2(synthetic: Path) -> Path:
    """The synthetic terrain meshed at --tolerance 2 (294 770 triangles)."""
    target = synthetic / "mesh2.vtk"
    ran("mesh", "--dem", str(synthetic / "synth.tif"), "--tolerance", "2", "--out", str(target))
    return target


class TestEndToEnd:
    """The demonstration case: terrain, mesh, two days of snow, opened by the readers."""

    @pytest.fixture(scope="class")
    def case(self, synthetic: Path) -> Path:
        ran(
            "mesh",
            "--dem",
            str(synthetic / "synth.tif"),
            "--tolerance",
            "20",
            "--out",
            str(synthetic / "mesh20.vtk"),
        )
        out = synthetic / "anim20"
        ran(
            "snow",
            "--mesh",
            str(synthetic / "mesh20.vtk"),
            "--out",
            str(out),
            "--sun-correction",
            "melt-factor",
            "--days",
            "2",
        )
        return out / "snow.case"

    def test_our_reader_opens_it(self, case: Path, synthetic: Path) -> None:
        mesh = read_vtk((synthetic / "mesh20.vtk").read_bytes())
        anim = Animation.read(case)
        assert len(anim.geometry.triangles) == len(mesh.polygons) > 1000
        assert anim.case.times == (0.0, 24.0, 48.0)
        swe = anim.frame(case.parent, "swe", 2)
        assert (swe <= 600.0).all() and (swe >= 0.0).all() and (swe < 600.0).any()
        slope = anim.static(case.parent, "slope_deg")
        assert slope.min() >= 0.0 and slope.max() < 90.0

    def test_vtk_opens_it(self, case: Path) -> None:
        vtk = pytest.importorskip("vtk")
        reader = vtk.vtkGenericEnSightReader()
        errors: list[int] = []
        reader.AddObserver("ErrorEvent", lambda *_: errors.append(1))
        reader.SetCaseFileName(str(case))
        reader.UpdateInformation()
        reader.SetTimeValue(48.0)
        reader.Update()
        block = reader.GetOutput().GetBlock(0)
        anim = Animation.read(case)
        assert block.GetNumberOfCells() == len(anim.geometry.triangles)
        swe = block.GetCellData().GetArray("swe")
        assert swe.GetValue(0) == anim.frame(case.parent, "swe", 2)[0]
        assert errors == []

    def test_paraview_opens_it(self, case: Path) -> None:
        if not PVPYTHON.exists():
            pytest.skip(f"no pvpython at {PVPYTHON} (set RASPUTIN_PVPYTHON)")
        script = (
            "from paraview.simple import OpenDataFile\n"
            "from paraview import servermanager as sm\n"
            f"r = OpenDataFile({str(case)!r})\nr.UpdatePipeline()\n"
            "ts = list(r.TimestepValues)\nr.UpdatePipeline(ts[-1])\nd = sm.Fetch(r)\n"
            "b = d.GetBlock(0) if d.IsA('vtkMultiBlockDataSet') else d\n"
            "print('READ', ts, b.GetNumberOfCells())\n"
        )
        done = subprocess.run(
            [str(PVPYTHON), "-c", script], capture_output=True, text=True, timeout=300
        )
        line = next((ln for ln in done.stdout.splitlines() if ln.startswith("READ")), "")
        assert line, done.stdout + done.stderr
        cells = len(Animation.read(case).geometry.triangles)
        assert line == f"READ [0.0, 24.0, 48.0] {cells}"


@pytest.mark.skipif(
    os.environ.get("RASPUTIN_SLOW") != "1", reason="slow (minutes): RASPUTIN_SLOW=1"
)
class TestE1:
    """Section 6's table, rows 1-3, from the CLI on the 2 m mesh. A miss is first
    checked against the design probe run on the same mesh."""

    def test_the_2_m_mesh_has_294_770_triangles(self, mesh2: Path) -> None:
        assert len(read_vtk(mesh2.read_bytes()).polygons) == 294_770

    @pytest.mark.parametrize("mode", list(TABLE))
    def test_the_spring_matches_section_6(self, mesh2: Path, synthetic: Path, mode: str) -> None:
        out = synthetic / f"e1-{mode}"
        ran("snow", "--mesh", str(mesh2), "--out", str(out), "--sun-correction", mode)
        got = reference_figures(mesh2, out)
        sca, swe, south, north = TABLE[mode]
        assert_allclose(got["sca"], sca, atol=0.01)
        assert_allclose(got["swe"], swe, rtol=0.02)
        assert got["south"] == pytest.approx(south, abs=1.0)
        assert got["north"] == pytest.approx(north, abs=1.0)


def dem_file(path: Path, *, x0: float = 499950.0, y0: float = 6800150.0, epsg: int = 25833) -> Path:
    """A flat 300 m DEM, 41 by 41 nodes 10 m apart, node (0, 0) at (x0, y0);
    by default it covers `facets()`."""
    path.write_bytes(encode_geotiff(np.full((41, 41), 300.0, dtype=np.float32), x0, y0, 10.0, epsg))
    return path


class TestShadowFlags:
    def test_dem_turns_shadows_on_by_default(self, mesh: Path, tmp_path: Path) -> None:
        out = tmp_path / "on"
        snow(mesh, out, "--dem", str(dem_file(tmp_path / "dem.tif")), *SPRING_10, "--days", "1")
        job = json.loads((out / "snow.json").read_text())["job"]
        assert job["shadows"] == RadiationRequest().model_dump(mode="json")

    def test_shadows_off_with_a_dem(self, mesh: Path, tmp_path: Path) -> None:
        out = tmp_path / "off"
        dem = str(dem_file(tmp_path / "dem.tif"))
        snow(mesh, out, "--dem", dem, "--shadows", "off", *SPRING_10, "--days", "1")
        assert json.loads((out / "snow.json").read_text())["job"]["shadows"] is None

    def test_no_dem_is_no_shadows(self, two_days: Path) -> None:
        assert json.loads((two_days / "snow.json").read_text())["job"]["shadows"] is None

    def test_a_flat_dem_around_the_facets_changes_little(
        self, mesh: Path, tmp_path: Path, two_days: Path
    ) -> None:
        # The flat DEM at the facets' own height hides only a sun below about
        # half a degree; the two days' mean dI moves by well under 1 W/m^2.
        out = tmp_path / "on"
        dem = str(dem_file(tmp_path / "dem.tif"))
        snow(mesh, out, "--dem", dem, "--sun-correction", "melt-factor", *SPRING_10, "--days", "2")
        on = Animation.read(out / "snow.case").static(out, "radiation_anomaly_mean")
        off = Animation.read(two_days / "snow.case").static(two_days, "radiation_anomaly_mean")
        assert (on <= off).all()
        assert_allclose(on, off, atol=1.0)


class TestShadowRefusals:
    """L5: each in plain words, before any file is written."""

    def test_shadows_on_without_a_dem(self, mesh: Path, tmp_path: Path) -> None:
        refused(
            mesh, tmp_path / "out", "--shadows", "on",
            says="shadows need the DEM the mesh was made from: give --dem",
        )  # fmt: skip

    def test_a_dem_in_another_crs(self, mesh: Path, tmp_path: Path) -> None:
        dem = dem_file(tmp_path / "dem.tif", x0=290000.0, y0=6800150.0, epsg=25832)
        output = refused(mesh, tmp_path / "out", "--dem", str(dem), says="EPSG:25832")
        assert squashed("EPSG:25833") in squashed(output)

    def test_a_dem_that_does_not_cover_the_mesh(self, mesh: Path, tmp_path: Path) -> None:
        dem = dem_file(tmp_path / "dem.tif", x0=600000.0)
        refused(mesh, tmp_path / "out", "--dem", str(dem), says="does not cover the mesh")

    def test_a_dem_that_covers_only_part_of_the_mesh(self, mesh: Path, tmp_path: Path) -> None:
        # The north facet starts at x 500100: a DEM ending at 500050 misses it.
        dem = dem_file(tmp_path / "dem.tif", x0=499650.0)
        refused(mesh, tmp_path / "out", "--dem", str(dem), says="does not cover the mesh")

    def test_an_unknown_shadows_value(self, mesh: Path, tmp_path: Path) -> None:
        refused(
            mesh, tmp_path / "out", "--shadows", "maybe", says="'maybe' is not one of 'on', 'off'"
        )


def _data_dir() -> Path | None:
    """DTM10's directory: `$RASPUTIN_DATA`, else `rasputin_data` beside a parent."""
    name = "DTM10_UTM33_20260925"
    if root := os.environ.get("RASPUTIN_DATA"):
        return Path(root) / name if (Path(root) / name).is_dir() else None
    for parent in Path(__file__).resolve().parents:
        if (parent / "rasputin_data" / name).is_dir():
            return parent / "rasputin_data" / name
    return None


def _plan_area(vtk_path: Path) -> npt.NDArray[np.float64]:
    mesh = read_vtk(vtk_path.read_bytes())
    pts, tri = mesh.points, np.array(mesh.polygons, dtype=np.int64)
    a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    return 0.5 * np.abs(np.cross(b - a, c - a)[:, 2])


@pytest.mark.skipif(
    os.environ.get("RASPUTIN_SLOW") != "1", reason="slow (minutes): RASPUTIN_SLOW=1"
)
class TestE2:
    """Section 10.8's E2: shadows on the 2 m mesh, one day, 2027-04-15."""

    def test_shadows_lower_the_days_mean_by_half_to_two_and_a_half_w(
        self, mesh2: Path, synthetic: Path
    ) -> None:
        means = {}
        for flag in ("on", "off"):
            out = synthetic / f"e2-{flag}"
            ran(
                "snow", "--mesh", str(mesh2), "--out", str(out), "--dem",
                str(synthetic / "synth.tif"), "--shadows", flag,
                "--start", "2027-04-15T00:00:00+00:00", "--days", "1",
            )  # fmt: skip
            means[flag] = Animation.read(out / "snow.case").static(out, "radiation_anomaly_mean")
        on, off = (means[f].astype(np.float64) for f in ("on", "off"))
        assert (on <= off).all()
        w = _plan_area(mesh2)
        change = float(w @ (on - off) / w.sum())
        # Probe S1, on 3000 DEM nodes: -1.1 W/m^2. Narrowed at the acceptance run.
        assert -2.5 <= change <= -0.5, change


@pytest.mark.skipif(
    os.environ.get("RASPUTIN_SLOW") != "1", reason="slow (minutes): RASPUTIN_SLOW=1"
)
class TestBygdin:
    """Section 10.9's demonstration: Bygdin on DTM10, 122 days, shadows on and off."""

    DOMAIN = Path(__file__).resolve().parents[2] / (
        "docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson"
    )

    @pytest.fixture(scope="class")
    def runs(self, tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, Path]:
        dtm10 = _data_dir()
        if dtm10 is None:
            pytest.skip("DTM10_UTM33_20260925 not found ($RASPUTIN_DATA or ../rasputin_data)")
        folder = tmp_path_factory.mktemp("bygdin")
        mesh = folder / "bygdin.vtk"
        ran(
            "mesh", "--dem", str(dtm10), "--domain", str(self.DOMAIN), "--tolerance", "10",
            "--out", str(mesh),
        )  # fmt: skip
        for flag in ("on", "off"):
            ran(
                "snow", "--mesh", str(mesh), "--dem", str(dtm10), "--shadows", flag,
                "--days", "122", "--out", str(folder / f"snow-{flag}"),
            )  # fmt: skip
        return mesh, folder / "snow-on", folder / "snow-off"

    def test_the_mesh_is_bygdin_at_tolerance_10(self, runs: tuple[Path, Path, Path]) -> None:
        # Tens of thousands of triangles, not an exact count: the count belongs
        # to the mesh's own tests and has moved between about 51,500 and 56,000
        # across increments; this only catches the wrong domain or tolerance.
        n = len(read_vtk(runs[0].read_bytes()).polygons)
        assert 20_000 <= n <= 120_000, n

    def test_122_days_of_daily_frames(self, runs: tuple[Path, Path, Path]) -> None:
        anim = Animation.read(runs[1] / "snow.case")
        assert anim.case.times == tuple(24.0 * d for d in range(123))

    def test_shadows_only_take_away(self, runs: tuple[Path, Path, Path]) -> None:
        mesh, on_dir, off_dir = runs
        on = Animation.read(on_dir / "snow.case").static(on_dir, "radiation_anomaly_mean")
        off = Animation.read(off_dir / "snow.case").static(off_dir, "radiation_anomaly_mean")
        assert (on <= off).all()
        w = _plan_area(mesh)
        # Probe B2 on DTM10 nodes: -2.7 to -3.3 W/m^2 on the three dates; the
        # sign and a loose band only, the figure is the acceptance run's.
        change = float(w @ (on.astype(np.float64) - off) / w.sum())
        assert -10.0 < change < -0.5, change
