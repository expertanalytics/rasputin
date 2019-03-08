from pathlib import Path

from rasputin import triangulate_dem
from rasputin.html_writer import write_mesh, color_field_by_height
from rasputin.reader import read_raster_file
from rasputin.triangulate_dem import lindstrom_turk_by_ratio


"""
Most minimal example
"""

raster_file = Path("../../rasputin_data/jotunheimen.tif")
assert raster_file.exists()

# Orig data set size
x0 = 432000.0
x1 = 491990.0 
y0 = 6790010.0 
y1 = 6850000.0

y1 = y0 + (y1 - y0)/40
x1 = x0 + (x1 - x0)/40

raster_coords, info = read_raster_file(filepath=raster_file,
                                       x0=x0,
                                       y0=y0,
                                       x1=x1,
                                       y1=y1)
points, faces = lindstrom_turk_by_ratio(raster_coords, 0.5)

normals = triangulate_dem.surface_normals(points, faces)

outputpath = Path.cwd() / "jotunheimen_web"
write_mesh(pts=points,
           faces=faces,
           normals=normals,
           output_dir=outputpath,
           vertex_field=color_field_by_height(points=points))
