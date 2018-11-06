from typing import Optional
from pathlib import Path
from rasputin import triangulate_dem


def write_ascii_off(pts: triangulate_dem.PointVector,
                    faces: triangulate_dem.FaceVector,
                    filepath: Path) -> None:
    with filepath.open("w") as tf:
        tf.write("OFF\n")
        tf.write(f"{len(pts)} {len(faces)} 0\n")
        for pt in pts:
            tf.write(f"{pt[0]} {pt[1]} {pt[2]}\n")
        for fc in faces:
            tf.write(f"3 {fc[0]} {fc[1]} {fc[2]}\n")


def write_vtk(pts: triangulate_dem.PointVector,
              faces: triangulate_dem.FaceVector,
              shade: Optional[triangulate_dem.IntVector],
              filepath: Path) -> None:
    with filepath.open("w") as tf:
        tf.write("# vtk DataFile Version 4.2\n")
        tf.write("vtk output\n")
        tf.write("ASCII\n")
        tf.write("DATASET UNSTRUCTURED_GRID\n")
        tf.write(f"POINTS {len(pts)} double\n")
        tf.write(" ".join([f"{pt[0]} {pt[1]} {pt[2]}" for pt in pts]))
        tf.write("\n\n")
        tf.write(f"CELLS {len(faces)} {4*len(faces)}\n")
        for fc in faces:
            tf.write(f"3 {fc[0]} {fc[1]} {fc[2]}\n")
        tf.write(f"CELL_TYPES {len(faces)}\n")
        for i in range(len(faces)):
            tf.write("5\n")
        if shade is not None:
            tf.write("\n")
            tf.write(f"CELL_DATA {len(shade)}\n")
            tf.write("SCALARS Shade double\n")
            tf.write("LOOKUP_TABLE default\n")
            tf.write(" ".join([f"{s}" for s in shade]))
            tf.write("\n")

