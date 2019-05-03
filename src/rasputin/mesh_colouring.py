import numpy as np
import seaborn as sns
from typing import Optional

from rasputin import triangulate_dem
from rasputin.avalanche import varsom_angles


def color_field_by_height(*, points: triangulate_dem.PointVector) -> np.ndarray:
    vertex_heights = np.asarray(points)[:, 2].copy()
    factor = 2
    vertex_heights -= min(vertex_heights)
    vertex_heights /= factor*max(vertex_heights)
    vertex_heights += 1 - 1/factor
    colors = np.empty((len(vertex_heights), 3))
    colors[:, 0] = vertex_heights
    colors[:, 1] = vertex_heights
    colors[:, 2] = vertex_heights
    return colors


def add_slope_colors(*,
                     normals: triangulate_dem.PointVector,
                     colors: np.ndarray) -> np.ndarray:
    slopes = np.asanyarray(triangulate_dem.compute_slopes(normals))
    colors[slopes >= 55/180*np.pi] = [0.2, 0.2, 0.2]
    return colors


def color_field_by_slope(*,
                         normals: triangulate_dem.PointVector) -> np.ndarray:
    slopes = np.asanyarray(triangulate_dem.compute_slopes(normals))
    colors = np.empty((len(slopes), 3))
    colors[slopes < 5.0e-2] = [0, 0, 1]
    colors[slopes > 5.0e-2] = [1, 1, 1]
    colors[slopes > 30/180*np.pi] = [1, 0, 0]
    return colors


def color_field_by_avalanche_danger(*,
                                    normals: triangulate_dem.PointVector,
                                    points: triangulate_dem.PointVector,
                                    faces: triangulate_dem.FaceVector,
                                    avalanche_problems: list,
                                    colors: Optional[np.ndarray]) -> np.ndarray:
    aspects = np.asarray(triangulate_dem.compute_aspects(normals))
    slopes = np.asarray(triangulate_dem.compute_slopes(normals))
    if colors is None:
        colors = np.ones((len(aspects), 3))
    else:
        colors = colors.copy()
    angles = np.asarray(varsom_angles)/180*np.pi
    for problem in avalanche_problems:
        expositions = [bool(int(s)) for s in problem["expositions"]]
        exposed_angles = angles[expositions]
        level = problem["level"]
        heights = problem["heights"]
        for i, face in enumerate(faces):
            avg_h = np.mean([points[f][2] for f in face])
            risk = False
            for hmin, hmax in heights:
                if hmin < avg_h < hmax:
                    aspect = aspects[i]
                    for a_min, a_max in exposed_angles:
                        if a_min < aspect < a_max:
                            if slopes[i] > 30/180*np.pi:
                                risk = True
                                break
                        elif a_min > a_max:
                            if a_min < aspect or aspect < a_max:
                                if slopes[i] > 30/180*np.pi:
                                    risk = True
                                    break
                if risk:
                    break
            if risk:
                if level == 1:
                    colors[i] = [0, 1, 0]
                elif level == 2:
                    colors[i] = [1, 1, 0]
                elif level == 3:
                    colors[i] = [1, 0, 0]
                elif level == 4:
                    colors[i] = [0, 0, 0]
    return colors


def color_field_by_aspect(*, normals: triangulate_dem.PointVector)->np.ndarray:

    palette = sns.color_palette("pastel", n_colors=len(varsom_angles))
    aspects = np.asanyarray(triangulate_dem.compute_aspects(normals))
    colors = np.empty((len(aspects), 3))
    for angle, col in zip(varsom_angles, palette):
        range = np.logical_and(aspects >= angle[0]/180*np.pi, aspects < angle[1]/180*np.pi)
        colors[range] = col
    return colors
