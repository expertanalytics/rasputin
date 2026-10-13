"""The EnSight Gold encoders: `tin_engine.io.ensight`, test R6.

`docs/increments/46-hbv-snow.md`, *4.6 The animation format: EnSight Gold,
geometry once*, and section 7's R6: exact byte sizes (244 + 4 T bytes per
value per frame; 2 587 724 for T = 646 870), a case file VTK's reader
accepts (time values one per line), coordinates stored as float32 offsets
from the mesh's minimum corner, and, with the `vtk` extra, a read back
through `vtkGenericEnSightReader`: the time values, a value planted at step 50
and the cell count. That last test skips without `vtk`.

Every file is also read by `ensightread`, the suite's own strict reader, which
takes every count from the bytes. Its layout is the probe writer's
(`46-probes/probe_sizes.py`), the one ParaView 6.1.1 read back with 92 time
steps; `TestReader` checks the reader on bytes laid out as that probe lays
them out, so the reader is known to accept what ParaView accepted.

Interface assumed (names chosen here and listed in the handback):

- `geometry_bytes(vertices, triangles, *, origin) -> bytes`: vertices
  `(n, 3)` float64, triangles `(T, 3)` 0-based; `origin` the `(x0, y0)`
  subtracted before the float32 cast, written in the second description line
  as `origin <x0> <y0>` (Python's `repr` of each float, so it parses back
  exactly). One part, number 1; `node id off`, `element id off`; one `tria3`
  block of 1-based corners. Size 732 + 12 n + 12 T bytes.
- `variable_bytes(name, values) -> bytes`: one float32 per triangle; the
  description line is `name`. Size 244 + 4 T bytes.
- `case_text(*, geometry, timed, static, times_hours) -> str`: `timed`
  variables on time set 1 with the file pattern `<name>****.var` (frames
  numbered from 0 by 1), `static` variables with no time set in `<name>.var`,
  one time value per line.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_array_equal

from ensightread import (
    VTK_CASE_LINE,
    Animation,
    EnsightError,
    read_case,
    read_geometry,
    read_variable,
)

# Two triangles at UTM 33 northings: float32 alone would lose 0.5 m here.
VERTICES = np.array(
    [
        [494000.25, 6806000.75, 1200.5],
        [494010.5, 6806000.75, 1201.25],
        [494000.25, 6806010.0, 1199.75],
        [494012.125, 6806011.5, 1202.0],
    ]
)
TRIANGLES = np.array([[0, 1, 2], [1, 3, 2]])
ORIGIN = (494000.25, 6806000.75)


@pytest.fixture(scope="module")
def ensight() -> ModuleType:
    return importlib.import_module("tin_engine.io.ensight")


def s80(text: str) -> bytes:
    return text.encode("ascii").ljust(80, b"\0")


def probe_geometry(points: npt.NDArray[np.float64], tri: npt.NDArray[np.int64]) -> bytes:
    """The probe writer's geometry layout, the one ParaView read back."""
    ox, oy = points[:, 0].min(), points[:, 1].min()
    return b"".join(
        [
            s80("C Binary"),
            s80("rasputin snow"),
            s80(f"origin {ox} {oy}"),
            s80("node id off"),
            s80("element id off"),
            s80("part"),
            np.int32(1).tobytes(),
            s80("mesh"),
            s80("coordinates"),
            np.int32(len(points)).tobytes(),
            (points[:, 0] - ox).astype("<f4").tobytes(),
            (points[:, 1] - oy).astype("<f4").tobytes(),
            points[:, 2].astype("<f4").tobytes(),
            s80("tria3"),
            np.int32(len(tri)).tobytes(),
            (tri + 1).astype("<i4").tobytes(),
        ]
    )


class TestReader:
    """The suite's reader accepts the probe's layout and refuses damage to it."""

    def test_it_reads_the_probe_geometry(self) -> None:
        g = read_geometry(probe_geometry(VERTICES, TRIANGLES))
        assert g.descriptions[1] == f"origin {ORIGIN[0]} {ORIGIN[1]}"
        assert_array_equal(g.triangles, TRIANGLES)
        assert g.points[1, 0] == np.float32(10.25)

    def test_it_refuses_a_truncated_or_padded_file(self) -> None:
        data = probe_geometry(VERTICES, TRIANGLES)
        with pytest.raises(EnsightError, match="ends at byte"):
            read_geometry(data[:-1])
        with pytest.raises(EnsightError, match="trailing bytes"):
            read_geometry(data + b"\0")

    def test_it_refuses_a_corner_past_the_coordinates(self) -> None:
        with pytest.raises(EnsightError, match="corner outside"):
            read_geometry(probe_geometry(VERTICES, np.array([[0, 1, 4]])))

    def test_it_refuses_a_case_line_vtk_cannot_read(self) -> None:
        values = " ".join(str(24 * q) for q in range(92))
        text = (
            "FORMAT\ntype: ensight gold\nGEOMETRY\nmodel: snow.geo\n"
            "TIME\ntime set: 1\nnumber of steps: 92\nfilename start number: 0\n"
            f"filename increment: 1\ntime values: {values}\n"
        )
        with pytest.raises(EnsightError, match="VTK reads"):
            read_case(text)


class TestGeometry:
    def test_its_size_is_732_plus_12_per_node_and_per_triangle(self, ensight: ModuleType) -> None:
        data = ensight.geometry_bytes(VERTICES, TRIANGLES, origin=ORIGIN)
        assert len(data) == 732 + 12 * len(VERTICES) + 12 * len(TRIANGLES)

    def test_coordinates_are_float32_offsets_from_the_origin(self, ensight: ModuleType) -> None:
        g = read_geometry(ensight.geometry_bytes(VERTICES, TRIANGLES, origin=ORIGIN))
        assert_array_equal(g.points[:, 0], (VERTICES[:, 0] - ORIGIN[0]).astype(np.float32))
        assert_array_equal(g.points[:, 1], (VERTICES[:, 1] - ORIGIN[1]).astype(np.float32))
        assert_array_equal(g.points[:, 2], VERTICES[:, 2].astype(np.float32))

    def test_the_offsets_keep_millimetres_at_utm_northings(self, ensight: ModuleType) -> None:
        # 0.75 m and 11.5 m above the origin: exact in float32 only as offsets.
        g = read_geometry(ensight.geometry_bytes(VERTICES, TRIANGLES, origin=ORIGIN))
        assert g.points[3, 1] == np.float32(10.75)
        assert np.float32(6806011.5) - np.float32(6806000.75) != np.float32(10.75)

    def test_the_origin_is_in_the_description_and_parses_back_exactly(
        self, ensight: ModuleType
    ) -> None:
        origin = (494000.1, 6806000.3)
        g = read_geometry(ensight.geometry_bytes(VERTICES, TRIANGLES, origin=origin))
        words = g.descriptions[1].split()
        assert words[0] == "origin"
        assert (float(words[1]), float(words[2])) == origin

    def test_corners_are_one_based_and_the_part_is_1(self, ensight: ModuleType) -> None:
        g = read_geometry(ensight.geometry_bytes(VERTICES, TRIANGLES, origin=ORIGIN))
        assert g.part == 1
        assert_array_equal(g.triangles, TRIANGLES)

    def test_an_empty_mesh_is_a_valid_file(self, ensight: ModuleType) -> None:
        empty = np.empty((0, 3))
        data = ensight.geometry_bytes(empty, np.empty((0, 3), dtype=np.int64), origin=(0.0, 0.0))
        g = read_geometry(data)
        assert len(g.points) == 0 and len(g.triangles) == 0


class TestVariable:
    def test_its_size_is_244_plus_4_per_triangle(self, ensight: ModuleType) -> None:
        values = np.arange(7, dtype=np.float64)
        assert len(ensight.variable_bytes("swe", values)) == 244 + 4 * 7

    def test_the_probe_meshs_frame_is_2_587_724_bytes(self, ensight: ModuleType) -> None:
        # Section 4.6: T = 646 870. Encoded for real, not computed.
        data = ensight.variable_bytes("swe", np.zeros(646_870))
        assert len(data) == 2_587_724

    def test_values_are_float32_on_part_1(self, ensight: ModuleType) -> None:
        values = np.array([0.1, 600.0, -1.0])
        v = read_variable(ensight.variable_bytes("sca", values), 3)
        assert v.description == "sca"
        assert v.part == 1
        assert_array_equal(v.values, values.astype(np.float32))


class TestCase:
    def case(self, ensight: ModuleType, frames: int) -> str:
        text: str = ensight.case_text(
            geometry="snow.geo",
            timed=("swe", "sca"),
            static=("slope_deg", "aspect_deg", "radiation_anomaly_mean", "melt_out_day"),
            times_hours=[24.0 * q for q in range(frames)],
        )
        return text

    def test_it_names_the_geometry_the_variables_and_the_times(self, ensight: ModuleType) -> None:
        case = read_case(self.case(ensight, 92))
        assert case.geometry == "snow.geo"
        assert case.timed == {"swe": "swe****.var", "sca": "sca****.var"}
        assert case.static == {
            "slope_deg": "slope_deg.var",
            "aspect_deg": "aspect_deg.var",
            "radiation_anomaly_mean": "radiation_anomaly_mean.var",
            "melt_out_day": "melt_out_day.var",
        }
        assert case.times == tuple(24.0 * q for q in range(92))
        assert (case.start, case.increment) == (0, 1)
        assert case.frame_file("swe", 50) == "swe0050.var"

    def test_the_times_are_one_per_line(self, ensight: ModuleType) -> None:
        # With the 92 values on one line VTK's reader fails ("'time values' not found").
        lines = self.case(ensight, 92).splitlines()
        at = next(i for i, line in enumerate(lines) if line.startswith("time values"))
        assert [float(line) for line in lines[at + 1 : at + 93]] == [24.0 * q for q in range(92)]

    def test_a_year_of_hourly_frames_keeps_every_line_short(self, ensight: ModuleType) -> None:
        text = self.case(ensight, 8761)
        assert max(len(line) for line in text.splitlines()) < VTK_CASE_LINE
        assert len(read_case(text).times) == 8761


class TestVtkReadsItBack:
    """R6 through VTK's EnSight reader, the one ParaView uses. Needs the `vtk` extra."""

    def test_time_values_a_planted_value_and_the_cell_count(
        self, ensight: ModuleType, tmp_path: Path
    ) -> None:
        vtk = pytest.importorskip("vtk")
        frames = 60
        (tmp_path / "snow.geo").write_bytes(
            ensight.geometry_bytes(VERTICES, TRIANGLES, origin=ORIGIN)
        )
        for q in range(frames):
            values = np.full(len(TRIANGLES), float(q))
            (tmp_path / f"swe{q:04d}.var").write_bytes(ensight.variable_bytes("swe", values))
        case = tmp_path / "snow.case"
        case.write_text(
            ensight.case_text(
                geometry="snow.geo",
                timed=("swe",),
                static=(),
                times_hours=[24.0 * q for q in range(frames)],
            ),
            encoding="ascii",
        )
        assert Animation.read(case).case.times[50] == 1200.0  # our reader agrees first

        reader = vtk.vtkGenericEnSightReader()
        errors: list[str] = []
        reader.AddObserver("ErrorEvent", lambda *_: errors.append("error"))
        reader.SetCaseFileName(str(case))
        reader.UpdateInformation()
        times = reader.GetTimeSets().GetItem(0)
        assert [times.GetValue(i) for i in range(times.GetNumberOfTuples())] == [
            24.0 * q for q in range(frames)
        ]
        reader.SetTimeValue(1200.0)
        reader.Update()
        output = reader.GetOutput()
        block = output.GetBlock(0)
        assert block.GetNumberOfCells() == len(TRIANGLES)
        swe = block.GetCellData().GetArray("swe")
        assert [swe.GetValue(i) for i in range(len(TRIANGLES))] == [50.0, 50.0]
        assert errors == []
