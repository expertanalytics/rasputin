from rasputin import reader                                                                                                                
from rasputin import geometry                                                                                                              
from shapely.geometry import Polygon
from rasputin.mesh import Mesh
from pathlib import Path

poly = Polygon([(59.887961, 10.859963), (59.886346, 10.861894), (59.886551, 10.867806), (59.889162, 10.86949), (59.889097, 10.863804)])

p = Path.home() / "rasputin_data/dem_archive/66m1_4_10m_z33.tif"

raster = reader.read_raster_file(filepath=p, polygon=poly)

Geopoly = geometry.Geopolygon()
