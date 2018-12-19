from typing import List, Tuple
import numpy as np
from rasputin import triangulate_dem


def compute_shade(*,
                  pts: triangulate_dem.PointVector,
                  faces: triangulate_dem.FaceVector,
                  sun_ray: List[float],
                  time: int) -> List[Tuple[int, np.ndarray]]:
    assert len(sun_ray) == 3
    shades = triangulate_dem.compute_shadow(pts, faces, sun_ray)
    has_shades = np.zeros(len(faces), dtype='int')
    has_shades[shades] = 1
    return [(time, has_shades)]
