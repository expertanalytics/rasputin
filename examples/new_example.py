import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

from shapely.geometry import Polygon
from rasputin.geometry import GeoPolygon, GeoPoints
from rasputin.reader import RasterRepository
from rasputin.gml_repository import GMLRepository, LandCoverType
from rasputin.mesh import Mesh
from pathlib import Path
from rasputin.writer import write_mesh
import os
import pyproj
import numpy as np



path = Path(os.environ["RASPUTIN_DATA_DIR"]) / "corine"
input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
target_crs = pyproj.CRS.from_epsg(32633)

x = np.array([8.5, 8.52, 8.52, 8.5])
y = np.array([60.55, 60.55, 60.50, 60.50])

x, y = pyproj.Transformer.from_crs(input_crs, target_crs).transform(x, y)
domain = GeoPolygon(polygon=Polygon(shell=list(zip(x, y))), crs=target_crs)

repos = GMLRepository(path=path)

dem_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "dem_archive"
rr = RasterRepository(directory=dem_archive)
raster_domain = domain.transform(target_crs=pyproj.CRS.from_proj4(rr.coordinate_system(domain=domain)))
raster_data_list = rr.read(domain=raster_domain)
mesh = Mesh.from_raster(data=raster_data_list, domain=raster_domain)

cmesh = mesh.simplify(ratio=0.1)
geo_cell_centers = GeoPoints(xy=cmesh.cell_centers[:, :2], crs=target_crs)
terrain_cover = repos.land_cover(land_types=None, geo_points=geo_cell_centers, domain=domain)


print(terrain_cover)


output_dir = Path.cwd() / "new_example_web"
write_mesh(pts=points, faces=faces, normals=normals, output_dir=output_dir)

