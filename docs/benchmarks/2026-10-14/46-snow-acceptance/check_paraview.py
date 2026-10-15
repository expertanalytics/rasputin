"""Increment 46 PR 3 acceptance: open snow.case in ParaView, read every frame,
and save one screenshot of a mid-spring frame coloured by snow water equivalent.

Run with ParaView's pvpython:
    pvpython check_paraview.py <dir>/snow.case <screenshot.png> <frame index>
"""

import sys

from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetScalarBar,
    OpenDataFile,
    SaveScreenshot,
    Show,
    UpdatePipeline,
)

case, png, mid = sys.argv[1], sys.argv[2], int(sys.argv[3])
reader = OpenDataFile(case)
times = list(reader.TimestepValues)
print(f"time steps: {len(times)}, first {times[0]}, last {times[-1]} hours")
print(f"cell arrays: {list(reader.CellData.keys())}")
bad = []
for i, t in enumerate(times):
    UpdatePipeline(time=t, proxy=reader)
    info = reader.GetDataInformation()
    rng = reader.CellData["swe"].GetRange()
    if info.GetNumberOfCells() != 646870:
        bad.append(i)
    if i in (0, mid, len(times) - 1):
        print(f"frame {i} t={t} h cells={info.GetNumberOfCells()} swe range={rng}")
print(f"frames read: {len(times)}, frames with a wrong cell count: {bad}")

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1400, 1100]
view.ViewTime = times[mid]
display = Show(reader, view)
ColorBy(display, ("CELLS", "swe"))
lut = GetColorTransferFunction("swe")
lut.ApplyPreset("Blues", True)
lut.RescaleTransferFunction(0.0, 600.0)
bar = GetScalarBar(lut, view)
bar.Title = f"SWE, mm (day {times[mid] / 24:.0f} of the spring)"
bar.ComponentTitle = ""
display.SetScalarBarVisibility(view, True)
view.OrientationAxesVisibility = 1
view.ResetCamera()
cam = view.GetActiveCamera()
cam.Elevation(-50)
cam.Azimuth(20)
view.ResetCamera()
cam.Zoom(1.5)
GetAnimationScene().AnimationTime = times[mid]
SaveScreenshot(png, view, ImageResolution=[1400, 1100])
print(f"screenshot: {png}")
