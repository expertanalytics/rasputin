from typing import List, Tuple, Union, Dict, Any
from pathlib import Path
import numpy as np
from meshio import XdmfTimeSeriesWriter, write_points_cells
from rasputin import triangulate_dem


def write_mesh(*,
               pts: triangulate_dem.point3_vector,
               faces: triangulate_dem.face_vector,
               shades: List[Tuple[int, np.ndarray]],
               filepath: Path) -> None:
    """Write mesh to file...

    :param pts:
    :param faces:
    :param shades:
    :param filepath:
    :return:
    """
    pts = np.asarray(pts)
    cells = {"triangle": np.asarray(faces)}
    if not shades:
        write_points_cells(str(filepath), pts, cells)
    else:
        writer = XdmfTimeSeriesWriter(str(filepath))
        writer.write_points_cells(pts, cells)
        for time, shade in shades:
            writer.write_data(time, cell_data={"triangle": {"shade": shade}})


class Writer(object):
    def __init__(self, *, filepath: Path):
        self._filepath = filepath
        self._point_data = dict()
        self._cell_data = dict()
        self._xdmfwriter = None

    def __enter__(self):
        self._xdmfwriter = XdmfTimeSeriesWriter(str(self._filepath))
        return self

    def __exit__(self, *args):
        self._xdmfwriter = None

    @property
    def xdmfwriter(self):
        if self._xdmfwriter is not None:
            return self._xdmfwriter
        else:
            raise RuntimeError("No file open!")

    def set_tin(self, pts, faces):
        self._pts = np.asarray(pts)
        self._faces = np.asarray(faces)

    def add_point_data(self, data, name):
        self._point_data[name] = np.asarray(data)

    def add_cell_data(self, data, name):
        self._cell_data[name] = np.asarray(data)

    def write_tin(self):
        self.xdmfwriter.write_points_cells(self._pts,
                                           {"triangle": self._faces})

    def write_all(self, t: float = 0.0):
        if not self.xdmfwriter.has_mesh:
            self.write_tin()
        self.xdmfwriter.write_data(t,
                                   cell_data={"triangle": self._cell_data},
                                   point_data=self._point_data)


def write(*,
          filepath: Path,
          pts: triangulate_dem.point3_vector,
          faces: triangulate_dem.face_vector,
          t: float = 0.0,
          fields: Dict["str", Any] = {}):
    """
    Convenience function for writing to file.
    """
    with Writer(filepath=filepath) as w:
        w.set_tin(pts, faces)
        for name, data in fields.items():
            # Try to determine if field has point or cell data
            # by counting number of values.
            # Can fail for very small meshes
            if len(data) == len(faces):
                w.add_cell_data(data, name)
            else:
                w.add_point_data(data, name)

        w.write_all(t=t)
