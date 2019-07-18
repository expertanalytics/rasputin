import os
import sys
import logging
from pathlib import Path
import numpy as np
import pprint
import pyproj
from shapely.geometry import Polygon
from shapely import wkt, wkb
import argparse

from rasputin.reader import RasterRepository
from rasputin.tin_repository import TinRepository
from rasputin.geometry import Geometry, GeoPoints, GeoPolygon
from rasputin.mesh import Mesh

from rasputin import gml_repository
from rasputin import globcov_repository


def read_poly_file(*, path: Path) -> Polygon:
    if path.suffix.lower() == ".wkb":
        with path.open("rb") as pfile:
            polygon = wkb.loads(pfile.read())
    elif path.suffix.lower() == ".wkt":
        with path.open("r") as pfile:
            polygon = wkt.loads(pfile.read())
    return polygon


def store_tin():
    """
    Avalance Forecast Visualization Example.

    Items to develop:
     * Construct several color maps for each avalanche danger, and make the web app offer selections
     * Alternative approach is to partition the whole area into disjoint topologies and switch between these
     * Use different textures for different terrain types (under development)
     * Use more of the information from varsom.no (and perhaps alpha blending) to better display avalanche dangers

    """

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-x", nargs="+", type=float,
                            help="x-coordinates of constraining polygon", default=None)
    arg_parser.add_argument("-y", nargs="+", type=float,
                            help="y-coordinates of constraining polygon", default=None)
    arg_parser.add_argument("-polyfile", type=str,
                            help="Constraining polygon definition in WKT or WKB format", default="")
    arg_parser.add_argument("-input-coordinate-system", type=str, default="EPSG:4326",
                            help="Coordinate system for constraining polygon (default: %(default)s)")
    arg_parser.add_argument("-target-coordinate-system", type=str, default="EPSG:32633",
                            help="Target coordinate system for generated tin (default: %(default)s)")
    arg_parser.add_argument("-ratio", type=float, default=0.4,
                            help="Mesh coarsening factor in (0, 1] (default: %(default)s)")
    arg_parser.add_argument("-override", action="store_true",
                            help="Replace existing archive entry")
    arg_parser.add_argument("-d", action="store_true",
                            help="Run in DEBUG mode")
    arg_parser.add_argument("-land-type-partition", type=str, default="", choices=["corine", "globcov"],
                            help="Partition mesh by land type. Requires additional downloaded data.")
    arg_parser.add_argument("-transpose", action="store_true",
                            help="Use the transpose of the raster image")
    arg_parser.add_argument("uid", type=str,
                            help="Unique ID for the result TIN")

    res = arg_parser.parse_args(sys.argv[1:])
    ll = logging.INFO
    if res.d:
        ll = logging.DEBUG
    logging.basicConfig(level=ll)
    logger = logging.getLogger("rasputin_store")
    if res.d:
        logger.setLevel(logging.DEBUG)

    if not (res.x and res.y or res.polyfile):
        logger.critical("\nEither both -x and -y, or -polyfile must be given to define constraining boundary.\n\n")
        arg_parser.print_help()
        sys.exit()

    # Set the data paths
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
        logger.warning(f"WARNING: No data directory specified, assuming dem_archive {dem_archive.absolute()}")
        logger.warning(f"WARNING: No data directory specified, assuming tin_archive {tin_archive.absolute()}")
        logger.warning(f"WARNING: No data directory specified, assuming globcov_archive {gc_archive.absolute()}")
        logger.warning(f"WARNING: No data directory specified, assuming corine_archive {corine_archive.absolute()}")

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
    if not res.x or not res.y:
        if not res.polyfile:
            raise RuntimeError("A constraining polygon is needed")
        else:
            source_polygon = read_poly_file(path=Path(res.polyfile))
    else:
        assert 3 <= len(res.x) == len(res.y), "x and y coordinates must have equal length greater or equal to 3"
        source_polygon = Polygon((x, y) for (x, y) in zip(res.x, res.y))

    target_coordinate_system = pyproj.Proj(init=res.target_coordinate_system)
    input_domain = GeoPolygon(polygon=source_polygon, projection=pyproj.Proj(init=res.input_coordinate_system))
    target_domain = input_domain.transform(target_projection=target_coordinate_system)

    # The transpose is probably flipped in the logic, so negating it should give correct behaviour.
    raster_repo = RasterRepository(directory=dem_archive, transpose=not res.transpose)
    raster_coordinate_system = pyproj.Proj(raster_repo.coordinate_system(domain=target_domain))
    raster_domain = input_domain.transform(target_projection=raster_coordinate_system)

    constraints = []
    lt_repo = None
    if res.land_type_partition:
        if res.land_type_partition == "corine":
            lt_repo = gml_repository.GMLRepository(path=corine_archive)
        else:
            lt_repo = globcov_repository.GlobCovRepository(path=gc_archive)
        constraints = lt_repo.constraints(domain=target_domain)
        if res.d:
            # Check stuff out.
            with open(f"domain.wkt", "w") as tf:
                tf.write(wkt.dumps(target_domain.polygon))
            for i, c in enumerate(constraints):
                with open(f"constr_{i}.wkt", "w") as tf:
                    tf.write(wkt.dumps(c.polygon))

    raster_data_list = raster_repo.read(domain=raster_domain)

    mesh = (Mesh.from_raster(data=raster_data_list,
                             domain=raster_domain)
            .simplify(ratio=res.ratio))

    assert len(mesh.points), "No tin extracted, something went wrong..."
    if raster_coordinate_system.definition_string() != target_coordinate_system.definition_string():
        points, faces = mesh.points, mesh.faces
        x, y, z = pyproj.transform(raster_coordinate_system,
                                   target_coordinate_system,
                                   points[:, 0],
                                   points[:, 1],
                                   points[:, 2])
        points = np.dstack([x, y, z])[0]
        mesh = Mesh.from_points_and_faces(points=points, faces=faces)

    if lt_repo is None:
        tr.save(uid=res.uid, geometries={"terrain": Geometry(mesh=mesh,
                                                             projection=target_coordinate_system,
                                                             base_color=(1.0, 1.0, 1.0),
                                                             material=None)})
    else:
        geometries = {}
        geo_cell_centers = GeoPoints(xy=mesh.cell_centers[:, :2],
                                     projection=target_coordinate_system)
        terrain_cover = lt_repo.land_cover(land_types=None,
                                           geo_points=geo_cell_centers,
                                           domain=target_domain)
        terrains = {lt.value: [] for lt in lt_repo.land_cover_type}

        for i, cell in enumerate(terrain_cover):
            terrains[cell].append(i)
        for t in terrains:
            if not terrains[t]:
                continue
            cover = lt_repo.land_cover_type(t)
            colors = [c/255 for c in lt_repo.land_cover_meta_info_type.color(land_cover_type=cover)]
            material = lt_repo.land_cover_meta_info_type.material(land_cover_type=cover)
            sub_mesh = mesh.extract_sub_mesh(np.array(terrains[t]))
            geometries[cover.name] = Geometry(mesh=sub_mesh,
                                              base_color=colors,
                                              material=material,
                                              projection=target_coordinate_system)
        tr.save(uid=res.uid, geometries=geometries)

    meta = tr.content[res.uid]
    logging.info(f"Successfully added uid='{res.uid}' to the tin archive {tin_archive.absolute()}, with meta info:")
    logging.info(pprint.pformat(meta, indent=4))
    #pprint.PrettyPrinter(indent=4).pprint(meta)
