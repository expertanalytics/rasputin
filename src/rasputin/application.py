import os
import sys
from pathlib import Path
import numpy as np
import pprint
import pyproj
from shapely.geometry import Polygon
from shapely import wkt, wkb
import argparse
from rasputin.reader import RasterRepository
from rasputin.tin_repository import TinRepository
from rasputin.triangulate_dem import lindstrom_turk_by_ratio, cell_centers, face_vector, point3_vector
from rasputin.geometry import Geometry, GeoPoints, GeoPolygon

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
        print(f"WARNING: No data directory specified, assuming dem_archive {dem_archive.absolute()}")
        print(f"WARNING: No data directory specified, assuming tin_archive {tin_archive.absolute()}")
        print(f"WARNING: No data directory specified, assuming globcov_archive {gc_archive.absolute()}")
        print(f"WARNING: No data directory specified, assuming corine_archive {corine_archive.absolute()}")

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
    res = arg_parser.parse_args(sys.argv[1:])

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

    raster_repo = RasterRepository(directory=dem_archive)
    raster_coordinate_system = pyproj.Proj(raster_repo.coordinate_system())
    target_coordinate_system = pyproj.Proj(init=res.target_coordinate_system)
    input_domain = GeoPolygon(polygon=source_polygon, projection=pyproj.Proj(init="EPSG:4326"))
    target_domain = input_domain.transform(target_projection=target_coordinate_system)
    raster_domain = input_domain.transform(target_projection=raster_coordinate_system)
    if res.land_type_partition:
        if res.land_type_partition == "corine":
            lt_repo = gml_repository.GMLRepository(path=corine_archive)
        else:
            lt_repo = globcov_repository.GlobCovRepository(path=gc_archive)
        constraints = lt_repo.constraints(domain=target_domain)
    else:
        constraints = []
    raster_data_list, cpp_polygon = raster_repo.read(domain=raster_domain)
    points, faces = lindstrom_turk_by_ratio(raster_data_list,
                                            cpp_polygon,
                                            res.ratio)
    assert len(points), "No tin extracted, something went wrong..."
    p = np.asarray(points)
    x, y, z = pyproj.transform(raster_coordinate_system,
                               target_coordinate_system,
                               p[:, 0],
                               p[:, 1],
                               p[:, 2])
    points = point3_vector.from_numpy(np.dstack([x, y, z])[0])

    if not res.land_type_partition:
        tr.save(uid=res.uid, geometries={"terrain": Geometry(points=points,
                                                             faces=faces,
                                                             projection=target_coordinate_system,
                                                             base_color=(1.0, 1.0, 1.0),
                                                             material=None)})
    else:

        geometries = {}
        constraints = lt_repo.constraints(domain=target_domain)
        tin_cell_centers = cell_centers(points, faces)
        geo_cell_centers = GeoPoints(xy=np.asarray(tin_cell_centers)[:, :2],
                                     projection=target_coordinate_system)
        terrain_cover = lt_repo.land_cover(land_types=None,
                                           geo_points=geo_cell_centers,
                                           domain=target_domain)
        terrains = {lt.value: face_vector() for lt in lt_repo.land_cover_type}

        for i, cell in enumerate(terrain_cover):
            terrains[cell].append(faces[i])
        for t in terrains:
            if not terrains[t]:
                continue
            cover = lt_repo.land_cover_type(t)
            colors = [c/255 for c in lt_repo.land_cover_meta_info_type.color(land_cover_type=cover)]
            material = lt_repo.land_cover_meta_info_type.material(land_cover_type=cover)
            geometries[cover.name] = Geometry(points=points,
                                              faces=terrains[t],
                                              base_color=colors,
                                              material=material,
                                              projection=target_coordinate_system).consolidate()
        tr.save(uid=res.uid, geometries=geometries)
    meta = tr.content[res.uid]
    print(f"Successfully added uid='{res.uid}' to the tin archive {tin_archive.absolute()}, with meta info:")
    pprint.PrettyPrinter(indent=4).pprint(meta)
