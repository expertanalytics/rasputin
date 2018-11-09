from typing import List, Tuple
from pathlib import Path
import numpy as np
from meshio import XdmfTimeSeriesWriter, write_points_cells
from rasputin import triangulate_dem


def write_mesh(*,
               pts: triangulate_dem.PointVector,
               faces: triangulate_dem.FaceVector,
               shades: List[Tuple[int, np.ndarray]],
               filepath: Path) -> None:
    pts = np.asarray(pts)
    cells = {"triangle": np.asarray(faces)}
    if not shades:
        write_points_cells(str(filepath), pts, cells)
    else:
        writer = XdmfTimeSeriesWriter(str(filepath))
        writer.write_points_cells(pts, cells)
        for time, shade in shades:
            writer.write_data(time, cell_data={"triangle": {"shade": shade}})
