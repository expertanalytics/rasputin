import os
from pathlib import Path
import numpy as np
import PIL
import pyproj

from rasputin import triangulate_dem
from rasputin.html_writer import (write_mesh, add_slope_colors, color_field_by_avalanche_danger)
from rasputin.reader import read_raster_file, extract_geo_keys, GeoKeysInterpreter
from rasputin.reader import RasterRepository
from rasputin.triangulate_dem import lindstrom_turk_by_ratio
from rasputin import avalanche


"""
Avalance Forecast Visualization Example.

Items to develop:
 * Accept an area of interest as start_x, start_y, stop_x, stop_y, with coordinate system
 * Construct several color maps for each avalanche danger, and make the web app offer selections
 * Alternative approach is to partition the whole area into disjoint toplogies and switch between these
 * Use better texturing for different terrain types (under development)
 * Use more of the information from varsom (and perhaps alpha blending) to better display avalanche dangers
 
"""
if "RASPUTIN_DATA_DIR" in os.environ:
    data_dir = Path(os.environ["RASPUTIN_DATA_DIR"])
else:
    data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"


# Define area of interest
# Orig data set size
x0 = 8.530918
y0 = 60.898468
input_coordinate_system = pyproj.Proj(init="EPSG:4326").definition_string()
target_coordinate_system = pyproj.Proj(init="EPSG:32633").definition_string()

# Distances in meters
dx = 5000.0
dy = 5000.0

dangers = []
expositions = "00000000"

tmp = avalanche.get_forecasts(x=x0, y=y0, proj=input_coordinate_system)[0]
avalanche_problems = []
if tmp:
    if "AvalancheProblems" in tmp:
        for p in tmp["AvalancheProblems"]:
            expositions = p["ValidExpositions"]
            level = p["AvalProbabilityId"]
            hf = p["ExposedHeightFill"]
            if hf == 1:
                min_h = p["ExposedHeight1"]
                max_h = 20000
                dangers = [(min_h, max_h)]
            elif hf == 2:
                min_h = 0
                max_h = p["ExposedHeight1"]
                dangers = [(min_h, max_h)]
            elif hf == 3:
                min_h1 = 0
                max_h1 = p["ExposedHeight1"]
                min_h2 = p["ExposedHeight2"]
                max_h2 = 20000
                dangers = [(min_h1, max_h1), (min_h2, max_h2)]
            elif hf == 4:
                min_h = p["ExposedHeight1"]
                max_h = p["ExposedHeight2"]
                dangers = [(min_h, max_h)]
            avalanche_problems.append({"expositions": expositions,
                                       "level": level,
                                       "heights": dangers})

raster_coords = RasterRepository(directory=data_dir).read(x=x0,
                                                          y=y0,
                                                          dx=dx,
                                                          dy=dy,
                                                          input_coordinate_system=input_coordinate_system,
                                                          target_coordinate_system=target_coordinate_system)

points, faces = lindstrom_turk_by_ratio(raster_coords, 0.5)
lakes, terrain = triangulate_dem.extract_lakes(points, faces)

normals = triangulate_dem.surface_normals(points, terrain)
point_normals = triangulate_dem.point_normals(points, terrain)

colors = np.ones(np.asarray(terrain).shape)
colors = add_slope_colors(normals=normals, colors=colors)
colors = color_field_by_avalanche_danger(normals=normals,
                                              points=points,
                                              faces=terrain,
                                              avalanche_problems=avalanche_problems,
                                              colors=colors)

outputpath = Path.cwd() / "general_web_small"
write_mesh(pts=points,
           faces=terrain,
           normals=point_normals,
           features=[lakes],
           output_dir=outputpath,
           face_field=colors)
