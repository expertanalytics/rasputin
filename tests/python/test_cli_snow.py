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
- The degenerate count is reported as "1 degenerate triangle".
- E1 runs only with `RASPUTIN_SLOW=1` (several minutes).
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
from tin_engine.io.vtk_legacy import write_vtk
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
) -> Path:
    vertices, triangles = facets(extra)
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
        assert_allclose(anim.static(two_days, "slope_deg"), [30.0, 30.0], atol=1e-3)
        assert_allclose(anim.static(two_days, "aspect_deg"), [180.0, 0.0], atol=1e-3)
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

    @pytest.fixture(scope="class")
    def mesh2(self, synthetic: Path) -> Path:
        target = synthetic / "mesh2.vtk"
        ran("mesh", "--dem", str(synthetic / "synth.tif"), "--tolerance", "2", "--out", str(target))
        return target

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
