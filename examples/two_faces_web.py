from pathlib import Path

from rasputin import triangulate_dem
from rasputin.html_writer import write_mesh


"""
Most minimal example
"""
vs = triangulate_dem.PointVector([(1,0,0), (0,1,0), (0,0,0), (0.25,0.25,1)])

points, faces = triangulate_dem.lindstrom_turk_by_ratio(vs, 2.0)
normals = triangulate_dem.surface_normals(points, faces)


output_dir = Path.cwd() / "two_faces_web"
write_mesh(pts=points, faces=faces, normals=normals, output_dir=output_dir)
