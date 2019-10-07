from pathlib import Path
import os
import numpy as np
import pyproj
from shapely.geometry import Polygon
from rasputin.mesh import Mesh
from rasputin.geometry import GeoPoints
from rasputin.gml_repository import GMLRepository, LandCoverType
from rasputin.geometry import GeoPolygon
from rasputin.reader import RasterRepository
from descartes import PolygonPatch
import matplotlib.pyplot as plt


def test_gml_repository():
    if "RASPUTIN_DATA_DIR" not in os.environ:
        raise RuntimeError("Please set RASPUTIN_DATA_DIR")
    path = Path(os.environ["RASPUTIN_DATA_DIR"]) / "corine"
    input_crs = pyproj.CRS.from_string("+init=EPSG:4326")
    target_crs = pyproj.CRS.from_epsg(32633)

    x = np.array([8.5, 8.52, 8.52, 8.5])
    y = np.array([60.55, 60.55, 60.50, 60.50])

    x, y = pyproj.Transformer.from_crs(input_crs, target_crs).transform(x, y)
    domain = GeoPolygon(polygon=Polygon(shell=list(zip(x, y))), crs=target_crs)

    repos = GMLRepository(path=path)
    plot = False
    if plot:
        response = repos.read(domain=domain)
        fig = plt.figure()
        ax = fig.gca()
        for key in response:
            color = [c/255 for c in LandCoverType.color(land_cover_type=key)]
            for p in response[key]:
                ax.add_patch(PolygonPatch(p, fc=color, alpha=0.5))
        ax.set_xbound(min(x), max(x))
        ax.set_ybound(min(y), max(y))
        plt.show()

    dem_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "dem_archive"
    rr = RasterRepository(directory=dem_archive)
    raster_domain = domain.transform(target_crs=pyproj.CRS.from_proj4(rr.coordinate_system(domain=domain)))
    raster_data_list = rr.read(domain=raster_domain)
    mesh = Mesh.from_raster(data=raster_data_list, domain=raster_domain)

    cmesh = mesh.simplify(ratio=0.1)
    geo_cell_centers = GeoPoints(xy=cmesh.cell_centers[:, :2], crs=target_crs)
    terrain_cover = repos.land_cover(land_types=None, geo_points=geo_cell_centers, domain=domain)
    assert terrain_cover is not None
