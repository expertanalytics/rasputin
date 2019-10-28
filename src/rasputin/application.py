import os
import sys
import logging
from pathlib import Path
import numpy as np
import pprint
import pyproj
from shapely.geometry import Polygon
import argparse
from datetime import datetime, timedelta

from rasputin.reader import RasterRepository
from rasputin.tin_repository import TinRepository, ShadeRepository
from rasputin.geometry import Geometry, GeoPoints, GeoPolygon
from rasputin.mesh import Mesh

from rasputin import gml_repository
from rasputin import globcov_repository


def store_tin():
    """
    Avalance Forecast Visualization Example.

    Items to develop:
     * Construct several color maps for each avalanche danger, and make the web app offer selections
     * Alternative approach is to partition the whole area into disjoint topologies and switch between these
     * Use different textures for different terrain types (under development)
     * Use more of the information from varsom.no (and perhaps alpha blending) to better display avalanche dangers

    """
    logging.basicConfig(level=logging.CRITICAL, format='Rasputin[%(levelname)s]: %(message)s')
    logger = logging.getLogger()


    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-x", nargs="+", type=float, help="x-coordinates of polygon", default=None)
    arg_parser.add_argument("-y", nargs="+", type=float, help="y-coordinates of polygon", default=None)
    arg_parser.add_argument("-polyfile", type=str, help="Polygon definition in WKT or WKB format", default="")
    arg_parser.add_argument("-target-coordinate-system", type=str, default="EPSG:32633", help="Target coordinate system")
    arg_parser.add_argument("-ratio", type=float, default=0.4, help="Mesh coarsening factor in [0, 1]")

    arg_parser.add_argument("-override", action="store_true", help="Replace existing archive entry")
    arg_parser.add_argument("-land-type-partition",
                            type=str,
                            default="",
                            choices=["corine", "globcov"],
                            help="Partition mesh by land type")
    arg_parser.add_argument("uid", type=str, help="Unique ID for the result TIN")
    arg_parser.add_argument("-silent", action="store_true", help="Run in silent mode")
    res = arg_parser.parse_args(sys.argv[1:])
    if not res.silent:
        logger.setLevel(logging.INFO)

    if "RASPUTIN_DATA_DIR" in os.environ:
        dem_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "dem_archive"
        tin_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "tin_archive"
        gc_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "globcov"
        corine_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "corine"
    else:
        #  data_dir = Path(os.environ["HOME"]) /"projects" / "rasputin_data" / "dem_archive"
        dem_archive = Path(".") / "dem_archive"
        tin_archive = Path(".") / "tin_archive"
        gc_archive = Path(".") / "globcov"
        corine_archive = Path(".") / "corine"
        logger.critical(f"WARNING: No data directory specified, assuming dem_archive {dem_archive.absolute()}")
        logger.critical(f"WARNING: No data directory specified, assuming tin_archive {tin_archive.absolute()}")
        logger.critical(f"WARNING: No data directory specified, assuming globcov_archive {gc_archive.absolute()}")
        logger.critical(f"WARNING: No data directory specified, assuming corine_archive {corine_archive.absolute()}")
    # Some sanity checks
    try:
        next(dem_archive.glob("*.tif"))
    except StopIteration as si:
        raise RuntimeError(f"No GeoTIFF files found in {dem_archive.absolute()}, giving up.")
    if not tin_archive.exists():
        tin_archive.mkdir(parents=True)
    elif tin_archive.exists() and not tin_archive.is_dir():
        raise RuntimeError(f"{tin_archive} exists and is not a directory, giving up.")

    tr = TinRepository(path=tin_archive)
    if res.uid in tr.content:
        if not res.override:
            raise RuntimeError(f"Tin archive {tin_archive.absolute()} already contains uid {res.uid}.")
        else:
            tr.delete(res.uid)

    # Determine region of interest
    if res.polyfile:
        input_domain = GeoPolygon.from_polygon_file(filepath=Path(res.polyfile),
                                                    crs=pyproj.CRS.from_string("+init=EPSG:4326"))

    elif (res.x and res.y):
        assert 3 <= len(res.x) == len(res.y), "x and y coordinates must have equal length greater or equal to 3"
        source_polygon = Polygon((x, y) for (x, y) in zip(res.x, res.y))
        input_domain = GeoPolygon(polygon=source_polygon,
                                  crs=pyproj.CRS.from_string("+init=EPSG:4326"))

    else:
        raise RuntimeError("A constraining polygon is needed")

    target_crs = pyproj.CRS.from_string(f"+init={res.target_coordinate_system}")
    target_domain = input_domain.transform(target_crs=target_crs)

    raster_repo = RasterRepository(directory=dem_archive)
    raster_crs = pyproj.CRS.from_string(raster_repo.coordinate_system(domain=target_domain))
    raster_domain = input_domain.transform(target_crs=raster_crs)

    raster_data_list = raster_repo.read(domain=raster_domain)

    mesh = (Mesh.from_raster(data=raster_data_list,
                             domain=raster_domain)
            .simplify(ratio=res.ratio))


    assert len(mesh.points), "No tin extracted, something went wrong..."
    if raster_crs.to_authority() != target_crs.to_authority():
        points, faces = mesh.points, mesh.faces
        proj = pyproj.Transformer.from_crs(raster_crs, target_crs)
        x, y, z = proj.transform(points[:, 0],
                                 points[:, 1],
                                 points[:, 2])
        points = np.dstack([x, y, z]).reshape(-1, 3)
        mesh = Mesh.from_points_and_faces(points=points, faces=faces, proj4_str=target_crs.to_proj4())

    if res.land_type_partition:
        if res.land_type_partition == "corine":
            lt_repo = gml_repository.GMLRepository(path=corine_archive)
        else:
            lt_repo = globcov_repository.GlobCovRepository(path=gc_archive)
        geo_cell_centers = GeoPoints(xy=mesh.cell_centers[:, :2],
                                     crs=target_crs)
        terrain_cover = lt_repo.land_cover(land_types=None,
                                           geo_points=geo_cell_centers,
                                           domain=target_domain)
        terrain_colors = np.empty((terrain_cover.shape[0], 3), dtype='d')
        extracted_terrain_types = set()
        for i, cell in enumerate(terrain_cover):
            if cell not in extracted_terrain_types:
                extracted_terrain_types.add(cell)
        meta_info = lt_repo.land_cover_meta_info_type
        for tt in extracted_terrain_types:
            cover = lt_repo.land_cover_type(tt)
            color = [c/255 for c in meta_info.color(land_cover_type=cover)]
            terrain_colors[terrain_cover == tt] = color
        tr.save(uid=res.uid,
                geometry=Geometry(mesh=mesh, crs=target_crs),
                land_cover_repository=lt_repo,
                face_fields={"cover_type": terrain_cover, "cover_color": terrain_colors})

    meta = tr.content[res.uid]
    logger.info(f"Successfully added uid='{res.uid}' to the tin archive {tin_archive.absolute()}, with meta info:")
    pprint.PrettyPrinter(indent=4).pprint(meta)


def compute_shades():
    """Compute shades for each cell in a given interval and with given frequency"""

    logging.basicConfig(level=logging.CRITICAL, format='Rasputin[%(levelname)s]: %(message)s')
    logger = logging.getLogger()
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-uid", type=str, help="uid of mesh")
    arg_parser.add_argument("-start-year", type=int, help="Start year")
    arg_parser.add_argument("-start-month", type=int, help="Start month")
    arg_parser.add_argument("-start-day", type=int, help="Start day")
    arg_parser.add_argument("-end-year", type=int, help="Start year")
    arg_parser.add_argument("-end-month", type=int, help="Start month")
    arg_parser.add_argument("-end-day", type=int, help="Start day")
    arg_parser.add_argument("-frequency", type=int, help="Frequency in seconds")
    arg_parser.add_argument("-silent", action="store_true", help="Run in silent mode")
    arg_parser.add_argument("-overwrite", action="store_true", help="Overwrite shade_uid")
    arg_parser.add_argument("shade_uid", type=str, help="Uid of saved shade")
    res = arg_parser.parse_args(sys.argv[1:])
    if not res.silent:
        logger.setLevel(logging.INFO)
    if "RASPUTIN_DATA_DIR" in os.environ:
        tin_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "tin_archive"
        shade_repo_archive = Path(os.environ["RASPUTIN_DATA_DIR"]) / "shade_archive"
    else:
        tin_archive = Path(".") / "tin_archive"
        shade_repo_archive = Path(".") / "shade_archive"
        logger.critical(f"WARNING: No data directory specified, assuming tin_archive {tin_archive.absolute()}")
    tin_repo = TinRepository(path=tin_archive)
    mesh = tin_repo.read(uid=res.uid).mesh
    start_time = datetime(res.start_year, res.start_month, res.start_day)
    end_time = datetime(res.end_year, res.end_month, res.end_day)
    dt = timedelta(seconds=res.frequency)
    t = start_time
    with ShadeRepository(path=shade_repo_archive,
                         tin_repo=tin_repo,
                         tin_uid=res.uid,
                         shade_uid=res.shade_uid,
                         overwrite=res.overwrite) as shade_archive:
        while t <= end_time:
            shade = mesh.shade(t.timestamp())
            shade_archive.save(t.timestamp(), shade)
            t += dt

