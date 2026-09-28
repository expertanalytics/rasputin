"""Read a rasputin .vtk back with VTK's vtkPolyDataReader (as the readback tests do)
and print what the 16b-1/2 acceptance checks: counts, the file's fields (R10), and
how many constraint edges carry land_cover (bit 7) and water (bit 8)."""
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

r = vtk.vtkPolyDataReader()
r.SetFileName(sys.argv[1])
r.ReadAllFieldsOn(); r.ReadAllScalarsOn()
r.Update()
pd = r.GetOutput()
print(f"VTK read: {pd.GetNumberOfPoints()} points, {pd.GetNumberOfPolys()} triangles, "
      f"{pd.GetNumberOfLines()} lines, z bounds {pd.GetBounds()[4]:.2f}..{pd.GetBounds()[5]:.2f}")
fd = pd.GetFieldData()
for i in range(fd.GetNumberOfArrays()):
    a = fd.GetAbstractArray(i)
    vals = [a.GetVariantValue(j).ToString() for j in range(min(a.GetNumberOfValues(), 8))]
    print(f"field {a.GetName()}: {' | '.join(vals)[:300]}")
cd = pd.GetCellData()
names = [cd.GetAbstractArray(i).GetName() for i in range(cd.GetNumberOfArrays())]
print("cell arrays:", ", ".join(names))
nl = pd.GetNumberOfLines()
# VTK orders cells lines first.
m = vtk_to_numpy(cd.GetArray("feature_mask")).astype(np.int64)[:nl]
lc, wa = (m >> 7) & 1, (m >> 8) & 1
print(f"constraint edges {nl}: land_cover {int(lc.sum())}, water {int(wa.sum())}, "
      f"water without land_cover {int((wa & (1 - lc)).sum())}, mask 0 {int((m == 0).sum())}")
for n in ("land_cover", "water"):
    arr = cd.GetArray(n)
    print(f"cell array {n}: " + ("absent" if arr is None else f"sum over lines {int(vtk_to_numpy(arr)[:nl].sum())}"))
