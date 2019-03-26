from pathlib import Path
import numpy as np
import PIL

from rasputin import triangulate_dem
from rasputin.html_writer import (write_mesh, add_slope_colors, color_field_by_avalanche_danger)
from rasputin.reader import read_raster_file, extract_geo_keys, GeoKeysInterpreter
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

# Define area of interest
# Orig data set size
x0 = 432000.0
x1 = 491990.0
y0 = 6790010.0
y1 = 6850000.0

# Select a smaller part of the raster
y1 = y0 + (y1 - y0)/10
x1 = x0 + (x1 - x0)/10

# TODO: Use the lat-lon from area of interest to read raster files from senorge files.
raster_file = Path("../../rasputin_data/jotunheimen.tif")
assert raster_file.exists()

interp = GeoKeysInterpreter(extract_geo_keys(image=PIL.Image.open(raster_file)))

dangers = []
expositions = "00000000"

tmp = avalanche.get_forecasts(x=x0, y=y0, proj=interp.to_proj4())[0]
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

raster_coords, info = read_raster_file(filepath=raster_file,
                                       x0=x0,
                                       y0=y0,
                                       x1=x1,
                                       y1=y1)
points, faces = lindstrom_turk_by_ratio(raster_coords, 0.5)

normals = triangulate_dem.surface_normals(points, faces)
point_normals = triangulate_dem.point_normals(points, faces)

# Commenting these out for now, as we need to partition the mesh into disjoint topologies for it to work.
#meshes = mesh_with_avalanche_danger(points=points, faces=faces, normals=normals,
#                                    avalanche_problems=avalanche_problems)

colors = np.ones(np.asarray(faces).shape)
colors = add_slope_colors(normals=normals, colors=colors)
colors = color_field_by_avalanche_danger(normals=normals,
                                              points=points,
                                              faces=faces,
                                              avalanche_problems=avalanche_problems,
                                              colors=colors)

outputpath = Path.cwd() / "jotunheimen_web"
write_mesh(pts=points,
           faces=faces,
           normals=point_normals,
           features=[],
           #features=meshes,
           output_dir=outputpath,
           #vertex_field=color_field_by_height(points=points),
           #vertex_field=color_field_by_height(points=points),
           face_field=colors)
