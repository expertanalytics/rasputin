import numpy as np
from rasputin import triangulate_dem


def compute_shade(pts, faces, sun_ray) -> np.ndarray:
    shades = triangulate_dem.compute_shadow(pts, faces, sun_ray)
    has_shades = np.zeros(len(faces), dtype='int')
    has_shades[shades] = 1
    return has_shades
