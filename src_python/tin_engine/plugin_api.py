"""What a plugin may import from ``tin_engine``: re-exports only.

A plugin package imports the core through this module and nothing else; every
other module of ``tin_engine`` may change in any commit. A name here is the
core's own object, never a copy or a wrapper, and it is promised: removed or
changed only in a release that says so. The list is what plugins use today and
grows on request.
"""

from tin_engine.catchment_core import accumulate
from tin_engine.cli_kit import THREADS, cli_resources, destination, words
from tin_engine.crs import frame_factors, parse_crs, reprojector, same_crs, site_of, srs_entry
from tin_engine.dem_input import DemRequest, open_dem, repository_for
from tin_engine.feature_input import FeatureError
from tin_engine.io.feature_table import FeatureTable, read_feature_table
from tin_engine.io.geopackage_write import create_features_table, create_geopackage, insert_features
from tin_engine.io.geotiff_write import encode_geotiff
from tin_engine.io.models import Bounds, DemTile, RasterMeta
from tin_engine.io.vtk_read import read_mesh
from tin_engine.raster import core_args
from tin_engine.resources import Resources
from tin_engine.run_record import plural
from tin_engine.sources import (
    RemoteSource,
    cache_root,
    catalogue,
    features_file,
    kept_named,
    plugin_source_files,
    user_source_file,
)

__all__ = [
    "THREADS",
    "Bounds",
    "DemRequest",
    "DemTile",
    "FeatureError",
    "FeatureTable",
    "RasterMeta",
    "RemoteSource",
    "Resources",
    "accumulate",
    "cache_root",
    "catalogue",
    "cli_resources",
    "core_args",
    "create_features_table",
    "create_geopackage",
    "destination",
    "encode_geotiff",
    "features_file",
    "frame_factors",
    "insert_features",
    "kept_named",
    "open_dem",
    "parse_crs",
    "plugin_source_files",
    "plural",
    "read_feature_table",
    "read_mesh",
    "repository_for",
    "reprojector",
    "same_crs",
    "site_of",
    "srs_entry",
    "user_source_file",
    "words",
]
