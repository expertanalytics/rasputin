import os
import sys
from pathlib import Path
import numpy as np
import pprint
import pyproj
import argparse
from rasputin.reader import RasterRepository
from rasputin.tin_repository import TinRepository
from rasputin.triangulate_dem import lindstrom_turk_by_ratio, extract_lakes, cell_centers, face_vector
from rasputin.geometry import Geometry, write_scene, lake_material, terrain_material
from rasputin.globcov_repository import GlobCovRepository, GeoPoints, LandCoverType


def store_tin():
    """
    Avalance Forecast Visualization Example.

    Items to develop:
     * Construct several color maps for each avalanche danger, and make the web app offer selections
     * Alternative approach is to partition the whole area into disjoint topologies and switch between these
     * Use different textures for different terrain types (under development)
     * Use more of the information from varsom.no (and perhaps alpha blending) to better display avalanche dangers

    """
    if "RASPUTIN_DATA_DIR" in os.environ:
        dem_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "dem_archive"
        tin_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "tin_archive"
        lt_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "globcov"
    else:
        #  data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"
        dem_archive = Path(".") / "dem_archive"
        tin_archive = Path(".") / "tin_archive"
        lt_archive = Path(".") / "globcov"
        print(f"WARNING: No data directory specified, assuming dem_archive {dem_archive.absolute()}")
        print(f"WARNING: No data directory specified, assuming tin_archive {tin_archive.absolute()}")
        print(f"WARNING: No data directory specified, assuming land_type_archive {lt_archive.absolute()}")
    try:
        next(dem_archive.glob("*.tif"))
    except StopIteration as si:
        raise RuntimeError(f"No GeoTIFF files found in {dem_archive.absolute()}, giving up.")
    if not tin_archive.exists():
        tin_archive.mkdir(parents=True)
    elif tin_archive.exists() and not tin_archive.is_dir():
        raise RuntimeError(f"{tin_archive} exists and is not a directory, giving up.")

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-lat", type=float, default=60.898468, help="Latitude of center coordinate")
    arg_parser.add_argument("-lon", type=float, default=8.530918, help="Longitude of center coordinate")
    arg_parser.add_argument("-dx", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-dy", type=float, default=5000, help="Distance in meters")
    arg_parser.add_argument("-ratio", type=float, default=0.4, help="Mesh coarsening factor in [0, 1]")
    arg_parser.add_argument("-override", action="store_true", help="Replace existing archive entry")
    arg_parser.add_argument("-land-type-partition", action="store_true", help="Partition mesh my land type")
    arg_parser.add_argument("uid", type=str, help="Unique ID for the result TIN")
    res = arg_parser.parse_args(sys.argv[1:])

    tr = TinRepository(path=tin_archive)
    if res.uid in tr.content:
        if not res.override:
            raise RuntimeError(f"Tin archive {tin_archive.absolute()} already contains uid {res.uid}.")
        else:
            tr.delete(res.uid)
    # Define area of interest
    x0 = res.lon
    y0 = res.lat
    dx = res.dx
    dy = res.dy

    input_coordinate_system = pyproj.Proj(init="EPSG:4326").definition_string()
    target_coordinate_system = pyproj.Proj(init="EPSG:32633").definition_string()

    raster_coords = RasterRepository(directory=dem_archive).read(x=x0,
                                                                 y=y0,
                                                                 dx=dx,
                                                                 dy=dy,
                                                                 input_coordinate_system=input_coordinate_system,
                                                                 target_coordinate_system=target_coordinate_system)
    points, faces = lindstrom_turk_by_ratio(raster_coords, res.ratio)
    if not res.land_type_partition:
        tr.save(uid=res.uid, geometries={"terrain": Geometry(points=points,
                                                             faces=faces, projection=target_coordinate_system,
                                                             base_color=(1.0, 1.0, 1.0),
                                                             material=None)})
    else:
        geometries = {}
        lakes, terrain = extract_lakes(points, faces)
        geometries["lake"] = Geometry(points=points,
                                      faces=lakes,
                                      base_color=(0, 0, 1),
                                      material=lake_material,
                                      projection=target_coordinate_system).consolidate()
        lc_repo = GlobCovRepository(path=lt_archive)
        tin_cell_centers = cell_centers(points, terrain)
        geo_cell_centers = GeoPoints(xy=np.asarray(tin_cell_centers)[:, :2],
                                     projection=pyproj.Proj(target_coordinate_system))
        terrain_cover = lc_repo.read_types(land_types=None, geo_points=geo_cell_centers)
        terrains = {lt.value: face_vector() for lt in LandCoverType}
        for i, cell in enumerate(terrain_cover):
            terrains[cell].append(terrain[i])
        for t in terrains:
            if not terrains[t]:
                continue
            geometries[LandCoverType(t).name] = Geometry(points=points,
                                                         faces=terrains[t],
                                                         base_color=[c/255 for c in LandCoverType.color(land_cover_type=LandCoverType(t))],
                                                         material=terrain_material,
                                                         projection=target_coordinate_system).consolidate()
        tr.save(uid=res.uid, geometries=geometries)
    meta = tr.content[res.uid]
    print(f"Successfully added uid='{res.uid}' to the tin archive {tin_archive.absolute()}, with meta info:")
    pprint.PrettyPrinter(indent=4).pprint(meta)
