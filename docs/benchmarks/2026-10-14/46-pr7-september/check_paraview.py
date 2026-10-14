"""Increment 46 PR 7 (copied from PR 6's check_paraview.py; the label takes a date): open the Bygdin runs with shadows on and off
in ParaView, read every frame of both, and save one screenshot of three
panels, seen from above with north up: SWE on a given day with shadows on,
the same day with shadows off, and their difference (on - off, mm, 0 to 250,
darker red where the shadows keep more snow; it is never negative).

Run with ParaView's pvpython:
    pvpython check_paraview.py <on dir> <off dir> <screenshot.png> <frame index> <date label>
"""

import sys

from paraview.simple import (
    AppendAttributes, AssignViewToLayout, Calculator, ColorBy, CreateLayout, CreateView,
    GetColorTransferFunction, GetScalarBar, OpenDataFile, SaveScreenshot, Show, UpdatePipeline,
)

on_dir, off_dir, png, mid, label = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
readers = {}
for name, d in (("on", on_dir), ("off", off_dir)):
    r = OpenDataFile(f"{d}/snow.case")
    times = list(r.TimestepValues)
    print(f"{name}: {len(times)} time steps, {times[0]} to {times[-1]} hours; "
          f"cell arrays {list(r.CellData.keys())}")
    cells = set()
    for t in times:
        UpdatePipeline(time=t, proxy=r)
        cells.add(r.GetDataInformation().GetNumberOfCells())
    UpdatePipeline(time=times[mid], proxy=r)
    print(f"{name}: frames read {len(times)}, cell counts {sorted(cells)}, "
          f"frame {mid} swe range {r.CellData['swe'].GetRange()}, "
          f"melt_out_day range {r.CellData['melt_out_day'].GetRange()}")
    readers[name] = (r, times)

both = AppendAttributes(Input=[readers["on"][0], readers["off"][0]])
UpdatePipeline(time=readers["on"][1][mid], proxy=both)
print(f"appended arrays: {list(both.CellData.keys())}")
diff = Calculator(Input=both, AttributeType="Cell Data", ResultArrayName="swe_shadow_minus_open",
                  Function="swe - swe_input_1")
t_mid = readers["on"][1][mid]
UpdatePipeline(time=t_mid, proxy=diff)
print(f"day {t_mid / 24:.0f}: swe on - off range {diff.CellData['swe_shadow_minus_open'].GetRange()}")

layout = CreateLayout("bygdin")
layout.SplitHorizontal(0, 1.0 / 3.0)
layout.SplitHorizontal(2, 0.5)
panels = ((readers["on"][0], "swe", "Blues", (0.0, 600.0), "SWE, mm, shadows on"),
          (readers["off"][0], "swe", "Blues", (0.0, 600.0), "SWE, mm, shadows off"),
          (diff, "swe_shadow_minus_open", "Reds", (0.0, 250.0), "SWE on - off, mm"))
views = []
for i, (src, arr, preset, rng, title) in enumerate(panels):
    view = CreateView("RenderView")
    AssignViewToLayout(view=view, layout=layout, hint=(1, 5, 6)[i])
    view.ViewTime = t_mid
    view.OrientationAxesVisibility = 0
    display = Show(src, view)
    ColorBy(display, ("CELLS", arr))
    lut = GetColorTransferFunction(arr)
    lut.ApplyPreset(preset, True)
    if preset == "Reds":
        lut.InvertTransferFunction()  # white at 0, dark red where the shadows keep most snow
    lut.RescaleTransferFunction(*rng)
    lut.AutomaticRescaleRangeMode = "Never"
    bar = GetScalarBar(lut, view)
    bar.Title = f"{title}, day {t_mid / 24:.0f} ({label})"
    bar.ComponentTitle = ""
    bar.WindowLocation = "Lower Center"
    bar.Orientation = "Horizontal"
    bar.ScalarBarLength = 0.6
    bar.TitleFontSize = 12
    bar.LabelFontSize = 11
    bar.RangeLabelFormat = "{:.0f}"
    display.SetScalarBarVisibility(view, True)
    view.ResetCamera()
    cam = view.GetActiveCamera()
    fp = cam.GetFocalPoint()
    cam.SetPosition(fp[0], fp[1], fp[2] + 1.0)
    cam.SetViewUp(0.0, 1.0, 0.0)
    view.ResetCamera()
    cam.Zoom(1.15)
    views.append(view)
SaveScreenshot(png, layout, ImageResolution=[2400, 900])
print(f"screenshot: {png}")
