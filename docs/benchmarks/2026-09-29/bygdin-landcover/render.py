"""Render the 10 m Bygdin mesh in natural land-cover colours, offscreen
(@perf, 2026-09-29, increment 16c acceptance item 5). Usage, from the
repository root, .venv active (the `viewer` extra's vtk):

    python docs/benchmarks/2026-09-29/bygdin-landcover/render.py <mesh.vtk> <out_dir>

Writes bygdin_landcover_oblique.png (from the south-east, 2x vertical
exaggeration) and bygdin_landcover_top.png (map view). The colours come from
tin_engine.palettes.CORINE_NATURAL, the table `rasputin palette corine` writes.
Only the triangles are drawn; constraint lines (code 0) are left out, since on
a shaded 3D surface they fight the triangles for depth.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
from vtkmodules.vtkCommonCore import vtkLookupTable, vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
from vtkmodules.vtkIOImage import vtkPNGWriter
from vtkmodules.vtkIOLegacy import vtkPolyDataReader
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor
from vtkmodules.vtkRenderingCore import (
    vtkActor, vtkPolyDataMapper, vtkRenderer, vtkRenderWindow, vtkWindowToImageFilter,
)  # fmt: skip
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401  (registers the OpenGL backend)
import vtkmodules.vtkRenderingFreeType  # noqa: F401  (text for the legend)

from tin_engine.palettes import CORINE_NATURAL

EXAGGERATION = 2.0


def surface(path: Path) -> tuple[vtkPolyData, list[int]]:
    """The triangles alone, centred on the origin, z exaggerated, with the
    triangle part of `land_cover_code`; and the codes present."""
    r = vtkPolyDataReader()
    r.SetFileName(str(path))
    r.Update()
    pd = r.GetOutput()
    n_lines = pd.GetNumberOfLines()
    codes = vtk_to_numpy(pd.GetCellData().GetArray("land_cover_code"))[n_lines:].astype(np.int32)
    pts = vtk_to_numpy(pd.GetPoints().GetData()).astype(np.float64)
    pts[:, :2] -= pts[:, :2].mean(axis=0)  # UTM offsets would cost float32 precision
    pts[:, 2] = (pts[:, 2] - pts[:, 2].min()) * EXAGGERATION
    tri = vtk_to_numpy(pd.GetPolys().GetConnectivityArray()).astype(np.int64).reshape(-1, 3)
    out = vtkPolyData()
    p = vtkPoints()
    p.SetData(numpy_to_vtk(pts, deep=True))
    out.SetPoints(p)
    cells = vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, 3 * len(tri) + 1, 3), deep=True),
                  numpy_to_vtkIdTypeArray(tri.ravel(), deep=True))  # fmt: skip
    out.SetPolys(cells)
    arr = numpy_to_vtk(codes, deep=True)
    arr.SetName("land_cover_code")
    out.GetCellData().SetScalars(arr)
    return out, sorted(set(codes.tolist()))


def lookup(present: list[int]) -> vtkLookupTable:
    """Indexed mode: one annotated colour per code present."""
    lut = vtkLookupTable()
    lut.IndexedLookupOn()
    lut.SetNumberOfTableValues(len(present))
    for i, code in enumerate(present):
        label, hexa = CORINE_NATURAL[code]
        rgb = [int(hexa[j : j + 2], 16) / 255 for j in (1, 3, 5)]
        lut.SetTableValue(i, *rgb, 1.0)
        lut.SetAnnotation(code, f"{code} {label}")
    lut.SetNanColor(1.0, 0.0, 1.0, 1.0)
    return lut


def render(poly: vtkPolyData, lut: vtkLookupTable, out: Path, oblique: bool) -> None:
    normals = vtkPolyDataNormals()
    normals.SetInputData(poly)
    normals.SplittingOff()
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.SetLookupTable(lut)
    mapper.SetScalarModeToUseCellData()
    mapper.UseLookupTableScalarRangeOn()
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetAmbient(0.25 if oblique else 0.45)
    actor.GetProperty().SetDiffuse(0.8 if oblique else 0.6)
    bar = vtkScalarBarActor()
    bar.SetLookupTable(lut)
    bar.SetNumberOfLabels(0)
    bar.SetMaximumWidthInPixels(60)
    # Annotations are drawn left of the swatches; leave them room there.
    bar.SetPosition(0.19, 0.52 if oblique else 0.04)
    bar.SetWidth(0.21)
    bar.SetHeight(0.44 if oblique else 0.40)
    bar.GetAnnotationTextProperty().SetColor(0.1, 0.1, 0.1)
    bar.GetAnnotationTextProperty().SetFontSize(14)
    bar.GetAnnotationTextProperty().ShadowOff()
    bar.SetTitle("")
    ren = vtkRenderer()
    ren.AddActor(actor)
    ren.AddViewProp(bar)
    if oblique:
        ren.GradientBackgroundOn()
        ren.SetBackground(0.93, 0.95, 0.97)
        ren.SetBackground2(0.62, 0.74, 0.88)
    else:
        ren.SetBackground(1.0, 1.0, 1.0)
    win = vtkRenderWindow()
    win.SetOffScreenRendering(1)
    win.SetSize(1400, 900)
    win.SetMultiSamples(8)
    win.AddRenderer(ren)
    x0, x1, y0, y1, z0, z1 = poly.GetBounds()
    cx, cy, cz = (x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2
    cam = ren.GetActiveCamera()
    cam.SetFocalPoint(cx + (0.05 * (x1 - x0) if oblique else 0), cy - (0.08 * (y1 - y0) if oblique else 0), cz)
    if oblique:
        # From the south-east, looking north-west, about 30 degrees down.
        d = 1.45 * max(x1 - x0, y1 - y0)
        cam.SetPosition(cx + 0.55 * d, cy - 0.75 * d, cz + 0.55 * d)
        cam.SetViewUp(0, 0, 1)
        cam.SetViewAngle(26)
    else:
        cam.SetFocalPoint(cx - 0.04 * (x1 - x0), cy, cz)
        cam.SetPosition(cx - 0.04 * (x1 - x0), cy, cz + 10 * (y1 - y0))
        cam.SetViewUp(0, 1, 0)
        cam.ParallelProjectionOn()
        cam.SetParallelScale(0.55 * (y1 - y0))
    ren.ResetCameraClippingRange()
    win.Render()
    grab = vtkWindowToImageFilter()
    grab.SetInput(win)
    grab.Update()
    png = vtkPNGWriter()
    png.SetCompressionLevel(9)
    png.SetFileName(str(out))
    png.SetInputConnection(grab.GetOutputPort())
    png.Write()


def main(mesh: Path, out_dir: Path) -> None:
    poly, present = surface(mesh)
    lut = lookup(present)
    render(poly, lut, out_dir / "bygdin_landcover_oblique.png", oblique=True)
    render(poly, lut, out_dir / "bygdin_landcover_top.png", oblique=False)
    print(f"codes drawn: {present}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
