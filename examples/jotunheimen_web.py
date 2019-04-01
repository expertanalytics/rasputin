from pathlib import Path
import numpy as np
import PIL

from rasputin import triangulate_dem
from rasputin.html_writer import (write_mesh, color_field_by_height, color_field_by_slope,
                                  color_field_by_aspect, add_slope_colors, color_field_by_avalanche_danger,
                                  mesh_with_avalanche_danger)
from rasputin.reader import read_raster_file, extract_geo_keys, GeoKeysInterpreter
from rasputin.triangulate_dem import lindstrom_turk_by_ratio
from rasputin import avalanche


"""
Most minimal example
"""

raster_file = Path("../../rasputin_data/jotunheimen.tif")
assert raster_file.exists()

interp = GeoKeysInterpreter(extract_geo_keys(image=PIL.Image.open(raster_file)))

# Orig data set size
x0 = 432000.0
x1 = 491990.0 
y0 = 6790010.0 
y1 = 6850000.0

y1 = y0 + (y1 - y0)/10
x1 = x0 + (x1 - x0)/10
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
lakes, terrain = triangulate_dem.extract_lakes(points, faces)

normals = triangulate_dem.surface_normals(points, terrain)
point_normals = triangulate_dem.point_normals(points, terrain)
face_field = color_field_by_slope(normals=normals)

#face_field2 = color_field_by_aspect(normals=normals)
#face_field2 = add_slope_colors(normals=normals, colors=face_field2)


#meshes = mesh_with_avalanche_danger(points=points, faces=faces, normals=normals,
#                                    avalanche_problems=avalanche_problems)

colors = np.ones(np.asarray(terrain).shape)

colors = add_slope_colors(normals=normals, colors=colors)
colors = color_field_by_avalanche_danger(normals=normals,
                                              points=points,
                                              faces=terrain,
                                              avalanche_problems=avalanche_problems,
                                              colors=colors)

outputpath = Path.cwd() / "jotunheimen_web"
write_mesh(pts=points,
           faces=terrain,
           normals=point_normals,
           features=[lakes],
           #features=meshes,
           output_dir=outputpath,
           #vertex_field=color_field_by_height(points=points),
           #vertex_field=color_field_by_height(points=points),
           face_field=colors)
