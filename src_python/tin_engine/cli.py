"""Command-line interface for the rasputin terrain engine.

This module is the single composition root: joining the file system to the
engine is its whole job. Paths from the command line also reach the domain,
DEM-input and repository modules, which read the files they name, while the
core never sees a file, a path or a CRS. Everything that has to know both sides lives here.
"""

from __future__ import annotations

import asyncio
import csv
import datetime
import functools
import hashlib
import json
import math
import os
import shlex
import sys
import time
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import closing, contextmanager
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, Final, Literal, TextIO, cast, get_args

import numpy as np
import numpy.typing as npt
import shapely
import typer
from pydantic import ValidationError
from shapely.geometry import LineString, Point, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient

from tin_engine import cli_kit, edge_strip, final_check, installed_version, sub_catchments
from tin_engine._core import (
    ChainRole,
    CheckPoints,
    IndexedMesh2,
    LineTolerance,
    NodedPslg,
    RefineOutcome,
    SlopeTolerance,
    build_pslg,
    describe,
    hardening,
    node,
    refine,
    sample,
    triangulate,
)
from tin_engine.budget import FLOOR_TOLERANCE, Trial, UnderFloor, Window, search
from tin_engine.catchment import (
    Catchment,
    CatchmentRequest,
    GaugeResult,
    LakeError,
    _lake,
    check_reach_crs,
    delineate,
)
from tin_engine.catchment_batch import (
    BatchRequest,
    NoRiverLine,
    StationResult,
    cell_at,
    run_batch,
    seed_for,
)
from tin_engine.cells import categories_of, cell_table, overlay
from tin_engine.chains import start_chains
from tin_engine.cli_kit import THREADS
from tin_engine.crs import (
    crs_label,
    parse_crs,
    reprojector,
    same_crs,
    srs_entry,
    transform_description,
    transform_label,
)
from tin_engine.decompose import DEFAULT_PIECE_BUDGET, Partition, bytes_per_node, partition
from tin_engine.decompose import lines as cut_lines
from tin_engine.dem_input import (
    CachedSource,
    DemPlan,
    DemRequest,
    NoAreaError,
    OutCrsRequiredError,
    PieceDem,
    TargetRegionError,
    assemble_dem,
    load_piece,
    piece_plan,
    plan_dem,
    repository_for,
)
from tin_engine.domain import BoxError, CrsBox, DomainError, DomainPolygon
from tin_engine.elevation import Trimmed, trim
from tin_engine.feature_input import (
    SUFFIXES,
    Category,
    ClassPolygons,
    EdgeChoice,
    EdgeCleanup,
    FeatureError,
    FeatureSet,
    LandCoverClasses,
    SourceLines,
    edges_of,
    read_class_polygons,
    read_lake_polygons,
    read_linestrings,
    read_polygons,
    read_source,
    source_region,
)
from tin_engine.features import (
    DEFAULT_VOCABULARY,
    PIECE_VOCABULARY,
    EdgeProperty,
    EdgeVocabulary,
)
from tin_engine.gauge import Gauge, Placement
from tin_engine.geodesy import area_m2
from tin_engine.gml_import import gml_to_geopackage
from tin_engine.grid_domain import default_stride, refine_start_stride, subsample
from tin_engine.hydrography import RiverSegment, Station
from tin_engine.io.cog import NotCached
from tin_engine.io.domain_file import read_domain
from tin_engine.io.geojson import catchment_geojson
from tin_engine.io.geopackage_write import (
    Palette,
    create_geopackage,
    write_cells,
    write_classes,
    write_fractions,
    write_land_cover,
    write_lines,
    write_run,
    write_styles,
    write_triangles,
)
from tin_engine.io.models import (
    BY_NAME,
    UNITS,
    ClassFractions,
    DemTile,
    IndexWindow,
    RasterMeta,
    ZSpelling,
    ZUnit,
)
from tin_engine.io.network_file import outlet_point, read_links, unit_geometries, units_at
from tin_engine.io.ply import field_comments, write_ply
from tin_engine.io.repository import (
    CacheRepository,
    CacheWriter,
    DemRepository,
    TiffDemRepository,
    create_sqlite,
)
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_nve_lakes, read_references, read_stations
from tin_engine.io.vtk_legacy import write_vtk
from tin_engine.land_cover_classes import derived_cleanup, known_codes, shipped_names
from tin_engine.land_cover_classes import load as load_classes
from tin_engine.land_cover_fractions import ClassRaster
from tin_engine.land_cover_fractions import assign as assign_fractions
from tin_engine.land_cover_map import LandCoverMap, build
from tin_engine.landcover import label_triangles
from tin_engine.legends import load_legend
from tin_engine.legends import shipped_names as legend_names
from tin_engine.mosaic import Seam, check_canvas
from tin_engine.network import (
    Network,
    check_area,
    compare,
    home_unit,
    join,
    outlet_distance_m,
)
from tin_engine.palettes import paraview_preset
from tin_engine.polygon_fractions import assign as assign_polygon_fractions
from tin_engine.preferred import areas_inside
from tin_engine.raster import to_core
from tin_engine.refine_job import How, Once, Refined, RunError, error_values
from tin_engine.resources import Resources, parse_bytes
from tin_engine.run_record import (
    SEAMS_AGREE,
    RunRecord,
    Value,
    as_json,
    escaped_ascii,
    file_fields,
    ordinal,
    plural,
    refined_record,
    stride_record,
    summary,
)
from tin_engine.sources import (
    SOURCES,
    STATION_SOURCES,
    VECTOR_SUFFIXES,
    RemoteSource,
    catalogue,
    check_values,
    data_root,
    kept_named,
    plugin_source_files,
    user_source_file,
    xdg_base,
)
from tin_engine.stats import PhaseClock, Report, Sizes, quality, render
from tin_engine.target_grid import Block, TargetGrid
from tin_engine.tolerance_field import ToleranceLines, ToleranceSlope, line_segments

if TYPE_CHECKING:  # the cut path's types; its modules load only for a cut run
    from tin_engine.basin_run import CutRun, Partitioned, Piece
    from tin_engine.io.mesh_index import PieceEntry

#: A cut run, given how to load a piece's DEM and write a finished piece.
type RunCut = Callable[[Callable[[IndexWindow], PieceDem], Callable[[Piece], PieceEntry]], CutRun]

#: What rasputin does, in one sentence: the root help and `pyproject.toml` say it.
DESCRIPTION: Final = (
    "Rasputin turns a digital elevation model (DEM) into a triangle mesh of the terrain "
    "that stays within a height tolerance you choose."
)
#: The panels of `rasputin --help`, in the order it shows them.
MESHING, CATCHMENTS, DATA, OTHER = "Meshing", "Catchments", "Data", "Other"

app = typer.Typer(name="rasputin", help=DESCRIPTION, no_args_is_help=True)

#: The one refusal when no cache directory can be found.
NO_CACHE: Final = (
    "no cache directory: pass --cache, or set RASPUTIN_DATA (the cache is then "
    "$RASPUTIN_DATA/cache) or XDG_CACHE_HOME (then $XDG_CACHE_HOME/rasputin)"
)
CACHE_DEFAULT: Final = (
    "Default: $RASPUTIN_DATA/cache, else $XDG_CACHE_HOME/rasputin, else ~/.cache/rasputin."
)
CACHE = Annotated[
    Path | None,
    typer.Option("--cache", help=f"The tile cache a catalogue --dem is read from. {CACHE_DEFAULT}"),
]
BBOX_HELP: Final = (
    "XMIN YMIN XMAX YMAX in --bbox-crs (required with it); x is easting or longitude "
    "whatever the CRS's axis order"
)
BBOX_CRS = Annotated[
    str | None,
    typer.Option(
        "--bbox-crs",
        help="The --bbox's CRS, anything pyproj reads; EPSG:4326 is longitude and latitude "
        "in degrees. Never assumed. Required with --bbox.",
    ),
]
NO_BOX_CRS: Final = (
    "--bbox needs --bbox-crs, the CRS its four numbers are in (EPSG:4326 for longitude and "
    "latitude in degrees); a box's CRS is never assumed"
)
SEED_CRS = Annotated[
    str | None,
    typer.Option(
        "--seed-crs",
        help="The --seed's CRS, anything pyproj reads; EPSG:4326 is longitude and latitude "
        "in degrees. Never assumed. Required with --seed.",
    ),
]
NO_SEED_CRS: Final = (
    "--seed needs --seed-crs, the CRS its two numbers are in (EPSG:4326 for longitude and "
    "latitude in degrees); a seed's CRS is never assumed"
)
DEM_HELP: Final = (
    "A GeoTIFF DEM, several (repeat --dem), one directory of tiles, or a source from "
    "rasputin sources"
)

MEMORY_BUDGET = Annotated[
    str | None,
    typer.Option(
        "--memory-budget",
        metavar="BYTES",
        help="The most one DEM window may take: a number of bytes, or with K, M, G or T "
        "(8G). Default: half the machine's memory, a container's limit counting as the "
        "machine's.",
    ),
]
Z_UNIT = Annotated[
    ZSpelling | None,
    typer.Option(
        "--z-unit",
        help="The DEM's height unit when the file does not say: m, ft or us-ft (US survey "
        "foot). Heights are converted to metres. Default: the file's, else metres.",
    ),
]

#: The chain roles' names, mapped onto the enum. A role outside it is a
#: ``KeyError`` here rather than a silently unconstrained chain downstream.
ROLES = {
    "outer": ChainRole.Outer,
    "hole": ChainRole.Hole,
    "breakline": ChainRole.Breakline,
}

#: Which roles name chains that are rings, and therefore carry a closing edge no
#: graph stores.
CORE_CLOSED_ROLES: tuple[object, ...] = (ChainRole.Outer, ChainRole.Hole)

#: Metres: input coordinates are UTM-shaped, so this is a millimetre on the
#: ground.
#:
#: A policy choice made here on the record, because there is no default spacing
#: anywhere in C++ and three of the nine noder statuses point at this number --
#: two of them in opposite directions. Coarser merges features that are genuinely
#: distinct (a road four millimetres from a river becomes one edge carrying both
#: bits); finer runs out of lattice, ``kMaxGridIndex`` being ``2**51``, so at
#: spacing ``s`` the representable coordinate is ``s * 2.25e15`` and a UTM
#: easting overflows well before ``1e-12``. A millimetre sits four orders of
#: magnitude clear of that.
DEFAULT_SNAP_SPACING = 1e-3
# The start mesh's minimum angle, degrees.
DEFAULT_START_MIN_ANGLE = 25.0
# The start mesh's minimum angle with land cover simplified, degrees.
LANDCOVER_START_MIN_ANGLE = 15.0
# The gain a start-quality point must bring, degrees; a negative gain adds every
# point.
DEFAULT_START_QUALITY_GAIN = 0.0
#: The four land-cover flags have their own help panel, so the long
#: --[no-]edges-merge-same-class does not narrow the main one's names.
CLEAN_UP = "Land-cover clean-up"
#: The three clean-up distances default to values derived from the class file.
DERIVED = (
    "Default: from the land-cover class file's positional accuracy and the DEM cell (CORINE: {})."
)
#: Their flags and the bound each given value stays under.
OPTION_TOPS = (
    ("--edges-repair", math.inf),
    ("--edges-max-shift", math.inf),
    ("--edges-outline-snap", 100.0),
)
#: The line kinds a --lines key may name: the edge names that are no area's.
LINE_KEYS = tuple(
    p.name for p in DEFAULT_VOCABULARY.properties if p.name not in ("land_cover", "water")
)


def main() -> None:
    """Keep Typer in sub-command mode.

    With a single registered command and no callback, Typer collapses that
    command into the root, and `rasputin version` becomes a usage error.
    """


def bounds_checks() -> str:
    """``off``, or ``on (<mode>)`` naming the checks the loaded ``_core`` has."""
    return "off" if hardening == "none" else f"on ({hardening})"


def version() -> None:
    """Print the installed rasputin version, then whether it is bounds-checked."""
    typer.echo(installed_version())
    typer.echo(f"bounds checks: {bounds_checks()}")
    from tin_engine.plugins import loaded  # then each plugin; the hook imports this module

    for line in loaded():
        typer.echo(line)


#: The install line for the ``mcp`` extra. Rich help reads ``[mcp]`` as a markup
#: tag and drops it, so the help text escapes the bracket exactly when the app
#: renders help through Rich (``rich_markup_mode`` is None without Rich, where
#: Click's plain help would print the backslash).
INSTALL_MCP = 'pip install "rasputin[mcp]"'
_INSTALL_MCP_HELP = INSTALL_MCP.replace("[", r"\[") if app.rich_markup_mode else INSTALL_MCP


MCP_HELP: Final = (
    "Serve the commands to an agent over the Model Context Protocol, on stdin and stdout."
    "\n\nOne tool per command, each call run in a process of its own. "
    "With --http, over Streamable HTTP instead: each user's bearer token from "
    "RASPUTIN_MCP_TOKENS (user=token pairs), calls run in <out-root>/<user>/work, and "
    "each output comes back as a signed download link valid 1 hour; files are kept 7 days. "
    f"Needs the mcp extra: {_INSTALL_MCP_HELP}"
)


def mcp(
    jobs: Annotated[
        int, typer.Option("--jobs", min=1, help="Tool calls run at once; one beyond this waits.")
    ] = 1,
    out_root: Annotated[
        Path | None,
        typer.Option(
            "--out-root",
            metavar="DIR",
            help="Tool calls may write only inside this folder, through .. and links too; "
            "missing folders inside it are made. Default: ~/rasputin-output.",
        ),
    ] = None,
    http: Annotated[
        bool, typer.Option("--http", help="Serve over Streamable HTTP, not stdin and stdout.")
    ] = False,
    host: Annotated[str, typer.Option("--host", help="With --http: the address to bind.")] = (
        "127.0.0.1"
    ),
    port: Annotated[int, typer.Option("--port", help="With --http: the port to bind.")] = 8765,
    public_url: Annotated[
        str | None,
        typer.Option(
            "--public-url",
            metavar="URL",
            help="With --http: the address clients use, which links are built on; https "
            "unless loopback. Default: http://<host>:<port>.",
        ),
    ] = None,
    quota: Annotated[
        str, typer.Option("--quota", help="With --http: each user's folder size, as 20G.")
    ] = "20G",
) -> None:
    """Serve the commands to an agent over the Model Context Protocol (``--help`` says how)."""
    try:
        from tin_engine.mcp import http as mcp_http
        from tin_engine.mcp import server
    except ImportError:
        typer.echo(f"Error: the mcp extra is not installed; run {INSTALL_MCP}", err=True)
        raise typer.Exit(2) from None
    from tin_engine import plugins  # the hook imports this module

    served = plugins.served or plugins.application()
    if not http:
        asyncio.run(server.serve(served, jobs, out_root))
        return
    from tin_engine.resources import parse_bytes

    try:
        settings = mcp_http.HttpSettings(
            host=host,
            port=port,
            public_url=public_url or f"http://{host}:{port}",
            out_root=out_root or Path.home() / server.DEFAULT_OUT_ROOT,
            quota=parse_bytes(quota),
        )
        mcp_http.serve_http(served, settings, mcp_http.tokens_from_env(os.environ), jobs=jobs)
    except ValueError as exc:
        typer.echo(f"Error: {cli_kit.words(exc)}", err=True)
        raise typer.Exit(2) from None


@dataclass(frozen=True, slots=True)
class _Run:
    """What ``build_pslg`` -> ``node`` -> ``triangulate`` made of some chains.

    ``noded`` is None when the validator or the noder refused; ``mesh`` is None
    when anything refused.
    """

    noded: NodedPslg | None
    mesh: IndexedMesh2 | None
    ok: bool
    status: str
    message: str


def _engine(
    vertices: npt.ArrayLike,
    chains: list[tuple[list[int], ChainRole, int]],
    delaunay: bool,
    spacing: float,
    clock: PhaseClock | None = None,
) -> _Run:
    """Validate, node and triangulate, reducing every refusal to words.

    ``clock`` records the three calls as ``start mesh: ...`` rows; without one
    a fresh clock is used and discarded."""
    clock = clock or PhaseClock()
    with clock.phase("start mesh: build"):
        result = build_pslg(np.asarray(vertices), chains)
    if not result.ok or result.pslg is None:
        return _Run(
            noded=None,
            mesh=None,
            ok=False,
            status=", ".join(d.error.name for d in result.diagnostics),
            message="; ".join(d.message for d in result.diagnostics),
        )
    with clock.phase("start mesh: node"):
        noded = node(result.pslg, spacing)
    if not noded.ok() or noded.pslg is None:
        return _Run(
            noded=None,
            mesh=None,
            ok=False,
            status=noded.status.name,
            message=_band(describe(noded.status), noded.message),
        )
    with clock.phase("start mesh: triangulate"):
        outcome = triangulate(noded.pslg, delaunay)
    return _Run(
        noded=noded.pslg,
        mesh=outcome.mesh if outcome.ok() else None,
        ok=outcome.ok(),
        status=outcome.status.name,
        message=_band(describe(outcome.status), outcome.message),
    )


def _band(sentence: str, detail: str) -> str:
    """The engine's prose and then the engine's detail, as one sentence."""
    return " ".join(part for part in (sentence, detail) if part)


def _snap_spacing(value: float) -> float:
    """A usage error, exited before the engine runs.

    NOT duplication of the engine's ``InvalidSnapSpacing``: the two answer
    different questions to different callers. This is a mis-typed option, and
    ``typer`` refuses it the same way a bad ``--out`` is refused. The C++ check
    is a library precondition with callers that are not this CLI, and neither
    can be deleted in favour of the other.
    """
    if not math.isfinite(value) or value <= 0.0:
        raise typer.BadParameter(f"snap-spacing must be finite and positive, got {value}")
    return value


#: What the suffix of ``--out`` selects.
MESH_SUFFIXES = (".vtk", ".ply", ".gpkg")


def _undirected(a: int, b: int) -> tuple[int, int]:
    """The canonical form of an edge, so the two directions are one key."""
    return (a, b) if a < b else (b, a)


def _masked_pairs(mesh: IndexedMesh2) -> set[tuple[int, int]]:
    """Every constrained mesh edge, deduplicated across the triangles sharing it.

    Bit ``e`` of a triangle's mask is the edge ``(v[e], v[(e + 1) % 3])``, per
    ``_core.pyi`` -- NOT CGAL's "edge ``e`` is opposite vertex ``e``", which is
    a rotation of it and would write a plausible file with every constraint on
    the wrong edge. An interior constraint is flagged by both its triangles,
    which is why this returns a set and not a list.
    """
    triangles = np.asarray(mesh.triangles)
    flags = np.asarray(mesh.constrained_edges)
    return {
        _undirected(int(triangle[e]), int(triangle[(e + 1) % 3]))
        for triangle, mask in zip(triangles, flags, strict=True)
        for e in range(3)
        if int(mask) >> e & 1
    }


def _chain_masks(pslg: NodedPslg) -> dict[tuple[int, int], int]:
    """Join each undirected node pair to the feature mask the noder gave it.

    ``edge_properties`` is dense and index-aligned with the flat edge
    enumeration: chain ``c``'s edge ``k`` sits at ``sum(edge_count(j) for j <
    c) + k``, and ``edge_count`` is the chain's vertex count for a ring and one
    less for a breakline -- a ring's closing edge is enumerated but not stored,
    so ``(k + 1) % n`` is unconditional. Getting that wrong shifts every mask
    after the first ring onto the wrong edge.

    A pair reached by two chains takes the UNION: the
    set means "every property of every chain that contributed geometry here",
    and union is the only merge whose answer does not depend on chain order.
    """
    properties = np.asarray(pslg.edge_properties)
    masks: dict[tuple[int, int], int] = {}
    at = 0
    for c, chain in enumerate(pslg.chains):
        walk = [int(i) for i in pslg.indices_of(c)]
        count = len(walk) if chain.role in CORE_CLOSED_ROLES else len(walk) - 1
        for k in range(count):
            pair = _undirected(walk[k], walk[(k + 1) % len(walk)])
            masks[pair] = masks.get(pair, 0) | int(properties[at + k])
        at += count
    return masks


def _constraint_arrays(
    mesh: IndexedMesh2, pslg: NodedPslg
) -> tuple[npt.NDArray[np.uint32], npt.NDArray[np.uint32]]:
    """The edge block: the mesh's own constrained edges, and their feature bits.

    A constrained mesh edge with no entry in the noded graph gets 0, which
    ``_core.pyi`` defines as *unclassified* rather than *wrong*.
    """
    masks = _chain_masks(pslg)
    pairs = sorted(_masked_pairs(mesh))
    return (
        np.array(pairs, dtype=np.uint32).reshape(-1, 2),
        np.array([masks.get(pair, 0) for pair in pairs], dtype=np.uint32),
    )


def mesh(
    out: Annotated[
        Path,
        typer.Option(
            "--out",
            writable=True,
            help="Where to write: .vtk for ParaView, .ply or .gpkg for QGIS. A .gpkg with "
            "--land-cover, its land-cover class file and its shares also holds the layer "
            "cells, what Shyft reads: one row per triangle with its land-type shares glacier, "
            "lake, reservoir, forest and unspecified (savanna, listed too, lies inside "
            "unspecified).",
        ),
    ],
    dem: Annotated[
        list[str],
        typer.Option(
            "--dem",
            help=f"{DEM_HELP}: mesh its extent with z sampled from it. A file named like a "
            "source is written as a path: ./glo30.",
        ),
    ],
    cache: CACHE = None,
    bbox: Annotated[
        tuple[float, float, float, float] | None,
        typer.Option(
            "--bbox",
            metavar="BOX",
            help=f"Mesh only {BBOX_HELP}; in another CRS than the mesh's, the box holding "
            "the image of its edges there is meshed, snapped outward to the DEM's nodes, or "
            "with --out-crs to its grid's.",
        ),
    ] = None,
    bbox_crs: BBOX_CRS = None,
    stride: Annotated[
        int | None,
        typer.Option("--stride", help="Every Nth DEM node. Default: <= 256 a side."),
    ] = None,
    out_edges: Annotated[
        Path | None,
        typer.Option(
            "--out-edges",
            writable=True,
            help="With .ply, write the mesh's edges to this PLY.",
        ),
    ] = None,
    binary: Annotated[
        bool,
        typer.Option(
            "--binary/--ascii",
            help="Binary files, small and fast; --ascii writes text you can read.",
        ),
    ] = True,
    delaunay: Annotated[
        bool, typer.Option("--delaunay/--no-delaunay", help="Triangulate with or without it.")
    ] = True,
    snap_spacing: Annotated[
        float,
        typer.Option(
            "--snap-spacing",
            callback=_snap_spacing,
            help="Snap grid spacing for the noder, in the input's units.",
        ),
    ] = DEFAULT_SNAP_SPACING,
    tolerance: Annotated[
        float | None,
        typer.Option(
            "--tolerance",
            help="Refine until every triangle is within this many metres "
            "of the DEM; --stride then sets the start grid. Default: no refinement.",
        ),
    ] = None,
    triangles: Annotated[
        int | None,
        typer.Option(
            "--triangles",
            min=1,
            metavar="N",
            help="Refine to about N triangles instead: the tolerance that gives them is "
            "found and reported, and the mesh is the one --tolerance writes at it.",
        ),
    ] = None,
    max_triangles: Annotated[
        int | None,
        typer.Option(
            "--max-triangles",
            metavar="M",
            help="With --triangles, never more than M triangles. Default: 5 % over N.",
        ),
    ] = None,
    rms: Annotated[
        bool,
        typer.Option(
            "--rms",
            help="Also report the root-mean-square height error over the DEM nodes "
            "inside the mesh; one more pass over the mesh. Always on with --triangles.",
        ),
    ] = False,
    tolerance_near: Annotated[
        tuple[Path, float] | None,
        typer.Option(
            "--tolerance-near",
            metavar="FILE N",
            help="Hold the mesh to N metres on these lines (.geojson, .gpkg or .gml). "
            "Needs --tolerance-ramp.",
        ),
    ] = None,
    tolerance_near_crs: Annotated[
        str | None,
        typer.Option(
            "--tolerance-near-crs", help="The --tolerance-near file's CRS when it does not say."
        ),
    ] = None,
    tolerance_ramp: Annotated[
        tuple[float, float] | None,
        typer.Option(
            "--tolerance-ramp",
            metavar="START END",
            help="From START metres from the lines, where N stops, to END metres, where "
            "--tolerance starts.",
        ),
    ] = None,
    tolerance_slope: Annotated[
        tuple[float, float, float] | None,
        typer.Option(
            "--tolerance-slope",
            metavar="N START END",
            help="Hold DEM nodes on steep ground to N metres: the general --tolerance up to "
            "START degrees of slope, N from END degrees, linear between (START = END: a step). "
            "The slope is the DEM's own, at each node.",
            rich_help_panel="Tolerance on steep ground",
        ),
    ] = None,
    domain: Annotated[
        Path | None,
        typer.Option(
            "--domain",
            help="With --tolerance, mesh only this polygon (.geojson, .json or "
            ".wkt), starting from its boundary alone.",
        ),
    ] = None,
    domain_crs: Annotated[
        str | None,
        typer.Option(
            "--domain-crs",
            help="The --domain file's CRS, anything pyproj reads; required for .wkt.",
        ),
    ] = None,
    out_crs: Annotated[
        str | None,
        typer.Option(
            "--out-crs",
            help="Mesh in this CRS: a DEM in another CRS (required for one in degrees, "
            "such as glo30) is resampled onto a square grid here, and with --tolerance "
            "checked against its own nodes. For example "
            "EPSG:3006 in Sweden, EPSG:25833 in Norway.",
        ),
    ] = None,
    start_min_angle: Annotated[
        float | None,
        typer.Option(
            "--start-min-angle",
            help="With --tolerance, add DEM nodes to the start mesh until its triangles "
            "have no angle under this many degrees, or a stated reason prevents it; "
            "0 is off, at most 35. Default: 25, or 15 with land cover simplified "
            "(--edges-max-shift above 0).",
        ),
    ] = None,
    start_quality_gain: Annotated[
        float | None,
        typer.Option(
            "--start-quality-gain",
            help="With --tolerance, add such a node only if it raises the smallest angle "
            "around it by at least this many degrees and, unless --no-constraint-feet, "
            "split a line it lies beyond; negative adds every node, with no line split "
            "even with constraint feet. "
            "At most 10. Default: 0.",
        ),
    ] = None,
    no_constraint_feet: Annotated[
        bool,
        typer.Option(
            "--no-constraint-feet",
            help="With --tolerance, insert a DEM node close to a constraint segment "
            "itself, instead of its foot on the segment.",
        ),
    ] = False,
    land_cover: Annotated[
        list[str] | None,
        typer.Option(
            "--land-cover",
            help="What is on the ground: a cached catalogue source of classes "
            f"({', '.join(k for k, s in SOURCES.items() if s.values == 'classes')}) or of "
            "class polygons "
            f"({', '.join(k for k, s in SOURCES.items() if s.values == 'features' and s.classes)}"
            "), a uint8 GeoTIFF or a file of class polygons. Without class borders in --edges "
            "it gives the exact share of each class in every triangle, in the mesh file. "
            "Repeatable for class polygons, read as one map: where two overlap, the first "
            "given keeps the area.",
        ),
    ] = None,
    land_cover_year: Annotated[
        int | None,
        typer.Option(
            "--land-cover-year",
            help="With --land-cover, the year of the classes. Default: the source's last year.",
        ),
    ] = None,
    land_cover_classes: Annotated[
        str | None,
        typer.Option(
            "--land-cover-classes",
            help="With --land-cover, the land-cover class file naming its classes: "
            f"{', '.join(shipped_names())}, or a land-cover class file (.toml). "
            "Default: the source's.",
        ),
    ] = None,
    land_cover_layer: Annotated[
        str | None,
        typer.Option(
            "--land-cover-layer", help="The --land-cover GeoPackage files' features table."
        ),
    ] = None,
    land_cover_crs: Annotated[
        str | None,
        typer.Option("--land-cover-crs", help="The --land-cover files' CRS when they do not say."),
    ] = None,
    prefer: Annotated[
        list[str] | None,
        typer.Option(
            "--prefer",
            metavar="KEY=SOURCE",
            help="With class polygons from --land-cover, take the class its class file names "
            "for KEY (lake, glacier) from SOURCE, a file or cached source of polygons, and the "
            "rest from --land-cover; the land cover's own polygons of that class go to their "
            "neighbours. Repeatable; where two sources overlap, the first given wins.",
        ),
    ] = None,
    lines: Annotated[
        list[str] | None,
        typer.Option(
            "--lines",
            metavar="KEY=SOURCE",
            help=f"Lines of one kind ({', '.join(LINE_KEYS)}) from SOURCE, a file or cached "
            "source of lines, for --edges. Repeatable, one source per kind.",
        ),
    ] = None,
    edges: Annotated[
        list[str] | None,
        typer.Option(
            "--edges",
            metavar="SELECTOR",
            help="With --domain, which borders and lines the mesh must follow: all, water, a "
            "category of the class file (lake, forest, ...), a class code, or a --lines kind. "
            "Repeatable. Default: none.",
        ),
    ] = None,
    shares: Annotated[
        bool,
        typer.Option(
            "--shares",
            help="With --land-cover and class borders in --edges, also the exact share of "
            "each class in every triangle, and in a .gpkg the cells layer.",
        ),
    ] = False,
    edges_repair: Annotated[
        float | None,
        typer.Option(
            "--edges-repair",
            rich_help_panel=CLEAN_UP,
            show_default=False,
            help="Metres: land-cover borders closer than this become one "
            "border, and gaps narrower than it are filled; 0 is off. " + DERIVED.format(1),
        ),
    ] = None,
    edges_merge_same_class: Annotated[
        bool,
        typer.Option(
            help="Drop borders between neighbouring land-cover polygons of one class. Default: on.",
            rich_help_panel=CLEAN_UP,
            show_default=False,
        ),
    ] = True,
    edges_max_shift: Annotated[
        float | None,
        typer.Option(
            "--edges-max-shift",
            rich_help_panel=CLEAN_UP,
            show_default=False,
            help="Metres: the farthest a land-cover border may move from where it was (after "
            "the repair and the outline rule) when simplified, each class keeping its area; "
            "simplifying can bring a border back within the outline-snap distance of the "
            "outline, but not closer than the repair distance. --lines lines are kept as given. "
            "0: no simplification, borders as read. " + DERIVED.format(50),
        ),
    ] = None,
    edges_outline_snap: Annotated[
        float | None,
        typer.Option(
            "--edges-outline-snap",
            rich_help_panel=CLEAN_UP,
            show_default=False,
            help="Metres, under 100: move land-cover borders this close to the domain "
            "outline onto it; 0 is off. " + DERIVED.format(5),
        ),
    ] = None,
    reservoirs: Annotated[
        str | None,
        typer.Option(
            "--reservoirs",
            help="Reservoir polygons, a file or cached source, splitting the cells' lake share.",
        ),
    ] = None,
    catchments: Annotated[
        str | None,
        typer.Option(
            "--catchments",
            metavar="UNITS",
            help="Sub-catchment polygons in a projected CRS, a file or cached source: those "
            "whose inner point lies in --domain are joined into the domain, their borders, "
            "simplified within --catchments-max-shift, become constraints, and each cell's "
            "catchment_id is its unit's id.",
        ),
    ] = None,
    catchments_max_shift: Annotated[
        float | None,
        typer.Option(
            "--catchments-max-shift",
            show_default=False,
            help="Metres: with --catchments, the farthest a unit border may move from where it "
            "was when simplified, each unit keeping its area. 0: no simplification, borders as "
            "read. Default: an eighth of a triangle's edge at --triangles.",
        ),
    ] = None,
    simplify_domain: Annotated[
        bool,
        typer.Option(
            "--simplify-domain",
            help="Simplify the --domain outline, each ring moving at most --domain-max-shift "
            "and the area kept. --catchments always simplifies its outline, within "
            "--catchments-max-shift.",
        ),
    ] = False,
    domain_max_shift: Annotated[
        float | None,
        typer.Option(
            "--domain-max-shift",
            show_default=False,
            help="Metres: with --simplify-domain, the farthest the outline may move from where "
            "it was when simplified, the area kept. 0: no simplification, outline as read. "
            "Default: an eighth of a triangle's edge at --triangles.",
        ),
    ] = None,
    threads: THREADS = None,
    memory_budget: MEMORY_BUDGET = None,
    z_unit: Z_UNIT = None,
    pieces: Annotated[
        int | None,
        typer.Option(
            "--pieces",
            min=1,
            metavar="N",
            help="With --domain, cut the run on DEM lattice lines into about this many pieces, "
            "written to <out>.pieces/ with index.json, and nothing at --out. Default: 1.",
        ),
    ] = None,
    piece_budget: Annotated[
        str | None,
        typer.Option(
            "--piece-budget",
            metavar="BYTES",
            help="The estimated memory one piece of a cut mesh may take; a domain estimated "
            "over it is cut into pieces. A number of bytes, or with K, M, G or T (16G). "
            "Default 16G, the same on every machine.",
        ),
    ] = None,
    stats: Annotated[
        Path | None,
        typer.Option(
            "--stats",
            writable=True,
            allow_dash=True,
            help="Also write sizes, quality and timings as Markdown to this path (.md "
            "recommended); - prints it on stdout after the path line(s).",
        ),
    ] = None,
    record_path: Annotated[
        Path | None,
        typer.Option(
            "--record",
            writable=True,
            allow_dash=True,
            help="Also write the run's record as JSON to this path, for reproducibility.",
        ),
    ] = None,
) -> None:
    """Mesh a DEM into a TIN and write it as .vtk, .ply or .gpkg, chosen by the suffix of --out.

    The DEM is one GeoTIFF, several, a directory of tiles, or a source: the
    package's or one in your sources file (see rasputin sources). --bbox or
    --domain cuts it. Without --tolerance every --stride-th DEM node is
    triangulated; with it the mesh is refined until every triangle is within
    that many metres of the DEM. --land-cover says what is on the ground,
    --prefer takes a class from another source and --lines adds rivers and
    roads; --edges chooses the borders and lines the mesh follows. Every
    triangle gets a land_cover_code, coloured by rasputin palette, and
    without class borders (or with --shares) each class's share. --pieces and
    --piece-budget cut a --domain run into pieces, each its own file under
    <out>.pieces/ with an index.json, which agree bit for bit along every seam.

    .vtk is one file for ParaView. .ply is the surface for QGIS, with the
    constraint lines in a second file given by --out-edges. .gpkg is a GeoPackage
    with triangle and line layers. Each carries what the run did; --stats and
    --record write it in full. A refused run exits non-zero and writes no file.
    """
    clock = PhaseClock()
    box = _crs_box(bbox, bbox_crs)
    res = cli_kit.cli_resources(threads, memory_budget)
    codes: npt.NDArray[np.int32] | None = None
    codes_text = ""
    table: ClassFractions | None = None
    palette: Palette | None = None
    cover: list[tuple[int, BaseGeometry]] = []
    seams: tuple[Seam, ...] = ()
    known = _known_sources(os.environ)
    root = cache_root(cache, os.environ)
    classes = _class_source(
        land_cover or [], land_cover_year, land_cover_classes, land_cover_layer, land_cover_crs,
        root, known, lambda value: _fetch_run(value, domain, box, out_crs, cache),
    )  # fmt: skip
    line_sources = _line_sources(
        lines or [], known, root, lambda value: _fetch_run(value, domain, box, out_crs, cache)
    )
    choice = _edge_choice(edges or [], classes, [s.key for s in line_sources])
    if edges and domain is None:
        raise typer.BadParameter("needs --dem, --domain and --tolerance", param_hint="--edges")
    for s in line_sources:
        if s.key not in choice.lines:
            raise typer.BadParameter(
                f"--lines {s.key} is followed by no --edges; add --edges {s.key}",
                param_hint="--lines",
            )
    if shares and classes is None:
        raise typer.BadParameter("needs --land-cover", param_hint="--shares")
    borders = bool(choice.every_class or choice.codes)  # the run follows class borders
    computing = classes is not None and (shares or not borders)  # the run gives the shares
    vocabulary = DEFAULT_VOCABULARY
    if catchments:
        catchment = EdgeProperty(name="catchment", bit=10)  # bit 9 is the cut's seam
        vocabulary = EdgeVocabulary(properties=(*DEFAULT_VOCABULARY.properties, catchment))
    categories = None  # the cells layer's, checked here before any DEM is read
    named = classes.classes if classes is not None else None
    if named is not None and computing and out.suffix == ".gpkg":
        try:
            categories = categories_of(named)
        except ValueError as exc:
            raise typer.BadParameter(str(exc), param_hint="--land-cover-classes") from exc
    no_cells = "--shares" if classes is not None and not computing else None
    reservoir = _reservoir_source(
        reservoirs, categories is not None, known, root,
        _fetch_run(reservoirs or "", domain, box, out_crs, cache), why=no_cells,
    )  # fmt: skip
    preferred = _prefer_sources(
        prefer or [], classes, known, root,
        lambda value: _fetch_run(value, domain, box, out_crs, cache),
    )  # fmt: skip
    units_source = _reservoir_source(
        catchments, categories is not None, known, root,
        _fetch_run(catchments or "", domain, box, out_crs, cache), "--catchments", no_cells,
    )  # fmt: skip
    if catchments is not None and domain is None:
        raise typer.BadParameter(
            "needs --domain: the units whose inner point lies in the domain are the sub-catchments",
            param_hint="--catchments",
        )
    for flag, shift, applies, needs in (
        ("--catchments-max-shift", catchments_max_shift, catchments is not None, "--catchments"),
        ("--domain-max-shift", domain_max_shift, simplify_domain, "--simplify-domain"),
    ):
        if shift is not None and not applies:
            raise typer.BadParameter(f"applies only with {needs}", param_hint=flag)
        if shift is not None and not (math.isfinite(shift) and shift >= 0):
            raise typer.BadParameter(f"must be finite, >= 0, got {shift}", param_hint=flag)
    if catchments is not None and triangles is None and catchments_max_shift is None:
        raise typer.BadParameter(
            "needs --triangles, whose budget sets the farthest a border may shift, or "
            "--catchments-max-shift",
            param_hint="--catchments",
        )
    if simplify_domain and domain is None:
        raise typer.BadParameter("needs --domain", param_hint="--simplify-domain")
    if simplify_domain and triangles is None and domain_max_shift is None:
        raise typer.BadParameter(
            "needs --triangles, whose budget sets the farthest the outline may shift, or "
            "--domain-max-shift",
            param_hint="--simplify-domain",
        )
    regulated = None  # each triangle's share under the reservoir polygons
    catchment_ids = None  # each cell's unit, with --catchments
    overrides = {
        "repair_m": edges_repair,
        "max_shift_m": edges_max_shift,
        "outline_snap_m": edges_outline_snap,
    }
    for (flag, top), metres in zip(OPTION_TOPS, overrides.values(), strict=True):
        if metres is not None and not (math.isfinite(metres) and 0 <= metres < top):
            under = "" if top == math.inf else f" and under {top:g}"
            raise typer.BadParameter(f"must be finite, >= 0{under}, got {metres}", param_hint=flag)
    if domain is None and domain_crs is not None:
        raise typer.BadParameter("applies only with --domain", param_hint="--domain-crs")
    near_lines = _tolerance_lines(tolerance, tolerance_near, tolerance_near_crs, tolerance_ramp)
    steep = _tolerance_slope(tolerance, tolerance_slope)
    if domain is not None and bbox is not None:
        raise typer.BadParameter("--bbox and --domain exclude each other", param_hint="--bbox")
    if out.suffix not in MESH_SUFFIXES:
        raise typer.BadParameter(
            f"{out} has no known suffix; use {' or '.join(MESH_SUFFIXES)}", param_hint="--out"
        )
    if out.suffix == ".vtk" and out_edges is not None:
        raise typer.BadParameter(
            "a .vtk file already carries the constraint edges", param_hint="--out-edges"
        )
    if out.suffix == ".gpkg" and out_edges is not None:
        raise typer.BadParameter(
            "the edges are the lines layer of the GeoPackage", param_hint="--out-edges"
        )
    if out.suffix == ".gpkg" and not binary:
        raise typer.BadParameter(
            "applies to .vtk and .ply; a GeoPackage is an SQLite file", param_hint="--ascii"
        )

    unstated = False  # heights read as metres without a word from anyone
    if stride is not None and stride < 1:
        raise typer.BadParameter(f"must be a positive integer, got {stride}", param_hint="--stride")
    if tolerance is not None and not (math.isfinite(tolerance) and tolerance >= 0):
        raise typer.BadParameter(
            f"must be finite and >= 0, got {tolerance}", param_hint="--tolerance"
        )
    if domain is not None and stride is not None:
        raise typer.BadParameter(
            "a domain is meshed from its boundary alone, not a stride grid",
            param_hint="--stride",
        )
    if tolerance is not None and triangles is not None:
        raise typer.BadParameter(
            "give one of --tolerance and --triangles", param_hint="--triangles"
        )
    if max_triangles is not None and (triangles is None or max_triangles <= triangles):
        raise typer.BadParameter(
            "needs --triangles, and must be over it", param_hint="--max-triangles"
        )
    window = None
    if triangles is not None:
        window = Window(
            triangles, max_triangles or max(math.floor(1.05 * triangles), triangles + 1)
        )
    refining = tolerance is not None or window is not None
    if domain is not None and not refining:
        raise typer.BadParameter(
            "--domain needs --tolerance or --triangles", param_hint="--tolerance"
        )
    if rms and tolerance is None and window is None:
        raise typer.BadParameter("--rms needs --tolerance", param_hint="--rms")
    if start_min_angle is not None and not (
        math.isfinite(start_min_angle) and 0 <= start_min_angle <= 35
    ):
        raise typer.BadParameter(
            f"must be finite, >= 0 and <= 35, got {start_min_angle}",
            param_hint="--start-min-angle",
        )
    if start_min_angle is not None and not refining:
        raise typer.BadParameter(
            "--start-min-angle needs --tolerance or --triangles", param_hint="--start-min-angle"
        )
    if start_quality_gain is not None and not (
        math.isfinite(start_quality_gain) and start_quality_gain <= 10
    ):
        raise typer.BadParameter(
            f"must be finite and <= 10, got {start_quality_gain}",
            param_hint="--start-quality-gain",
        )
    if start_quality_gain is not None and not refining:
        raise typer.BadParameter(
            "--start-quality-gain needs --tolerance or --triangles",
            param_hint="--start-quality-gain",
        )
    if no_constraint_feet and not refining:
        raise typer.BadParameter(
            "--no-constraint-feet needs --tolerance or --triangles",
            param_hint="--no-constraint-feet",
        )
    flag = "--pieces" if pieces is not None else "--piece-budget"
    cuts = pieces is not None or (piece_budget is not None and window is None)
    if cuts and (tolerance is None or domain is None):
        raise typer.BadParameter("cutting needs --tolerance and --domain", param_hint=flag)
    try:
        budget = DEFAULT_PIECE_BUDGET if piece_budget is None else parse_bytes(piece_budget)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--piece-budget") from exc
    if tolerance is not None and budget < (one := bytes_per_node(tolerance)):
        raise typer.BadParameter(
            f"{budget} bytes is under one node's {one} at this tolerance",
            param_hint="--piece-budget",
        )
    # A cut run's index at this --out is stale now, whatever this run is.
    Path(f"{cli_kit.destination(out)}.pieces", "index.json").unlink(missing_ok=True)
    paths, cached, entry = _dem_sources(dem, cache_root(cache, os.environ), known)
    caller = _z_unit(z_unit) or (entry.z_unit if entry else None)
    given = None
    if domain is not None:
        with clock.phase("domain read"):
            try:
                given = read_domain(domain, domain_crs)
            except DomainError as exc:
                raise typer.BadParameter(str(exc), param_hint="--domain") from exc
    units: Units | None = None
    if given is not None and units_source is not None:
        with clock.phase("catchments"), _refused("--catchments"):
            units = _units(*units_source[:3], given, triangles, catchments_max_shift)
        given = units.domain
    outline_shift = None
    if given is not None and units is None and simplify_domain:
        with clock.phase("domain simplify"):
            before = _ring_vertices(given.polygon)
            given, outline_shift = sub_catchments.simplify_domain(
                given, triangles, domain_max_shift
            )
            area = sub_catchments.equal_area(given).polygon.area
        why = "--domain-max-shift" if domain_max_shift is not None else (
            f"a triangle's edge at {triangles} triangles, {8 * outline_shift:.0f} m, over 8"
        )  # fmt: skip
        typer.echo(
            f"domain: outline simplified at {outline_shift:.0f} m ({why}): {before} to "
            f"{_ring_vertices(given.polygon)} vertices, {area / 1e6:.2f} km2",
            err=True,
        )
    try:
        planned = _open_dem(paths, box, clock, given, cached, out_crs, res, caller)
        meta, part, grid = planned.meta, None, planned.grid
        if tolerance is not None and domain is not None:  # one cell is the uncut run
            cells = partition(
                meta.cols, meta.rows, tolerance, pieces=pieces or 1, piece_budget=budget
            )
            if cells.nx * cells.ny > 1:
                part = _cut_partition(planned, tolerance, cells, pieces or 1, budget)
        if part is not None:  # each of these describes one mesh, which a cut run has not
            asked = (
                ("--stats", stats), ("--out-edges", out_edges), ("--record", record_path),
                ("--land-cover", land_cover if computing and not shares else None),
                ("--shares", shares or None), ("--out", out if out.suffix == ".gpkg" else None),
                ("--rms", rms or None),
            )  # fmt: skip
            for hint, given_option in asked:
                if given_option is not None:
                    raise typer.BadParameter("not written for a cut run yet", param_hint=hint)
            from tin_engine.basin_run import piece_windows

            with _dem_refusals(paths):  # every piece within the cap before any is read
                plans = {w: piece_plan(planned, w) for w in piece_windows(part).values()}
                for p in sorted(
                    plans.values(), key=lambda p: -p.mosaic.meta.rows * p.mosaic.meta.cols
                ):
                    check_canvas(p.mosaic, res.memory_budget)
                planned.repository.check(planned.mosaic)
            held: Once[DemTile] | RasterMeta = meta
            checks, opened_seams = None, None
        else:
            with clock.phase("decode"), _dem_refusals(paths):
                opened = assemble_dem(planned, res)
            # What is read later comes out of `opened`, which then goes, so the
            # tile's one reference is the holder `_dem_mesh` empties.
            checks, opened_seams, held = opened.checks, opened.seams, Once(opened.tile)
            del opened
    except NotCached as exc:  # an outline past the domain's box needs the window's blocks
        hint = _bbox_fetch(dem[0], exc, cache) if simplify_domain or catchments else ""
        hint = hint or f"{RUN_FETCH}{_fetch_run(dem[0], domain, box, out_crs, cache)}"
        raise typer.BadParameter(f"{exc}{hint}", param_hint="--dem") from exc
    plan, source_crs, dem_domain = planned.mosaic, planned.source_crs, planned.domain
    dem_crs = meta.crs
    # The clean-up scale is the class file's accuracy when the run follows
    # class borders (CORINE's 100 m when it does not), never finer than a DEM cell.
    cell = max(meta.delta_x, abs(meta.delta_y))
    lcc = classes.classes if classes is not None else None
    accuracy = lcc.accuracy_m if borders and lcc is not None and lcc.accuracy_m else 100.0
    derived = derived_cleanup(accuracy, cell).model_dump()
    cleanup = {k: derived[k] if v is None else v for k, v in overrides.items()}
    cleanup["merge_same_class"] = edges_merge_same_class
    # The read region: around the domain, else around the DEM window meshed.
    area = dem_domain or DomainPolygon(polygon=shapely.box(*meta.node_box()), crs=dem_crs)
    land: _LandCover | None = None
    edge_values: dict[str, Value] = {}
    if classes is not None and classes.files:
        land = _land_cover_map(classes, preferred, area, dem_crs, cell, clock, res.threads)
        if len(classes.files) > 1:
            edge_values["land_cover_overlap"] = _overlap_record(land, classes.files, area)
    line_pairs, line_reads = _read_lines(line_sources, area, dem_crs, cell, clock, res.threads)
    found = None
    if edges and dem_domain is not None:
        found = _edges(
            land if borders else None, line_pairs, choice, EdgeCleanup(**cleanup), dem_domain,
            vocabulary, source_region(area, dem_crs, dem_crs, cell=cell), clock,
            sum(r.empty for r in line_reads),
        )  # fmt: skip
        edge_values |= _edge_record(
            found, land if borders else None, classes.files if classes else (), line_sources,
            line_pairs, dem_crs,
        )  # fmt: skip
        if borders:
            edge_values |= {f"edges_{k}": v for k, v in cleanup.items()}
            edge_values["edges_merge_same_class"] = "on" if edges_merge_same_class else "off"
    if units is not None and dem_domain is not None:
        units = units.moved(dem_crs)
        line = sub_catchments.borders([g for _, g in units.units], vocabulary.mask("catchment"))
        found = FeatureSet(features=(line,)) if found is None else found.model_copy(
            update={"features": (*found.features, line)}
        )  # fmt: skip
    # 15 degrees with land cover simplified, unless given.
    simplified = (
        found is not None and found.cover_vertices is not None and cleanup["max_shift_m"] > 0
    )
    angle = LANDCOVER_START_MIN_ANGLE if simplified else DEFAULT_START_MIN_ANGLE
    angle = angle if start_min_angle is None else start_min_angle
    try:
        dem_run = _dem_mesh(
            held,
            ", ".join(map(str, dem)),
            stride,
            delaunay,
            snap_spacing,
            tolerance,
            clock,
            dem_domain,
            domain.name if domain is not None else "",
            angle,
            not no_constraint_feet,
            found,
            grid,
            checks,
            DEFAULT_START_QUALITY_GAIN if start_quality_gain is None else start_quality_gain,
            near_lines,
            steep,
            res.threads,
            vocabulary,
            part,
            window,
            budget if domain is not None else None,
            rms,
        )
    except RunError as exc:
        raise typer.BadParameter(str(exc), param_hint="--dem") from exc
    except typer.BadParameter as exc:
        plain = domain is not None and outline_shift is None and units is None
        if isinstance(exc.__cause__, UnderFloor) and plain:
            exc.message += "; or give --simplify-domain to simplify the outline to the budget"
        raise
    surface_mesh, meta, values = dem_run.trimmed, dem_run.meta, dict(dem_run.values)
    names = [t.name for t in plan.tiles]
    size = f"{meta.cols} columns x {meta.rows} rows"
    if len(names) > 1:
        typer.echo(f"DEM: {len(names)} files, {size}", err=True)
    unstated = meta.vertical_unit_assumed and caller is None
    source = cached.source if cached is not None else escaped_ascii("; ".join(names))
    values |= {
        "crs": dem_crs,
        "dem_source": source,
        "dem_vertical_unit": _z_words(meta, z_unit, caller),
    }
    if meta.file_crs is not None:  # a DEM in feet, read as its metre twin
        feet = BY_NAME[cast(ZUnit, parse_crs(meta.file_crs).axis_info[0].unit_name)][3]
        values["dem_file_crs"] = (
            f"{meta.file_crs}, in {feet}; read as {source_crs or dem_crs}, the same "
            "projection in metres"
        )
    if grid is None:
        dx, dy = meta.delta_x, meta.delta_y
        apart = f"{dx:g} m" if dx == dy else f"{dx:g} m x {dy:g} m"
        values["dem_grid"] = f"{size}, {apart} apart"
    else:  # resampled onto a square grid in --out-crs
        how = transform_description(source_crs, dem_crs)
        values |= {"dem_crs": source_crs, "dem_transform": how}
        values["resampled_grid"] = f"{grid.spacing:g} m square grid in {dem_crs}, {size}"
    if entry is not None:  # the notes the source asks to travel with it
        values["dem_credit"] = escaped_ascii(entry.credit)
        values["licence_note"] = escaped_ascii(entry.licence_note)
        values["cite"] = escaped_ascii("; ".join(entry.cite)) or None
    if cached is not None or len(paths) > 1 or paths[0].is_dir():  # the files used
        values["dem_tiles"] = escaped_ascii("; ".join(names))
        if opened_seams is not None:  # a cut run's are per piece, and reach no file
            seams = opened_seams
            said = "; ".join(s.entry() for s in seams)
            values["dem_seams"] = escaped_ascii(said) or SEAMS_AGREE
    if given is not None:
        how = transform_label(given.crs, dem_crs)
        values |= {"domain_crs": crs_label(given.crs), "domain_transform": how}
        if outline_shift is not None:
            values["domain_max_shift_m"] = outline_shift
    if found is not None and dem_run.feature_counts is not None:
        chains, feature_vertices = dem_run.feature_counts
        totals = f"{plural(chains, 'line', 'lines')}, "
        totals += plural(feature_vertices, "vertex", "vertices")
        values["edges_joined"] = f"{totals} after joining shared edges and adding crossings"
    values |= edge_values
    coded = lcc if borders else None  # codes by the borders' polygons, no shares
    if coded is not None and found is not None:
        palette = coded.legend.colours
        cover = [(f.code, f.polygon) for f in found.features if f.polygon and f.code]
        if dem_run.cut is None and not computing:
            codes, codes_text = _land_cover(dem_run.trimmed, found, coded, snap_spacing, clock)
            values["land_cover_codes"] = codes_text
            values["land_cover_shares"] = "not computed (add --shares)"
    if land is not None and preferred and dem_run.cut is None:
        values |= _prefer_record(land, preferred, dem_run.trimmed, classes)
    if classes is not None and computing:
        table, codes, more = _class_fractions(
            classes, land, dem_run.trimmed, dem_crs, clock, domain, out_crs, cache, res.threads,
        )  # fmt: skip
        codes_text = str(more["land_cover_codes"])
        if lcc is not None:  # names and colours for the .gpkg
            palette = lcc.legend.colours
        elif out.suffix == ".gpkg":
            more["cells"] = "none, the land cover has no land-cover class file"
        if reservoir is not None:
            with clock.phase("reservoirs"):
                xy = np.asarray(dem_run.trimmed.vertices, dtype=np.float64)
                try:
                    found_polygons, repaired = read_polygons(
                        *reservoir[:2],
                        dem_crs,
                        xy[:, :2],
                        functools.partial(_say_repaired, "reservoirs"),
                        reservoir[3],
                    )
                except FeatureError as exc:
                    raise typer.BadParameter(str(exc), param_hint="--reservoirs") from exc
                regulated = overlay(xy[dem_run.trimmed.triangles], found_polygons)
            more |= {"reservoirs": escaped_ascii(reservoir[2]), "reservoirs_repaired": repaired}
        elif categories is not None:
            more["reservoirs"] = "none"
        if units is not None:
            xy = np.asarray(dem_run.trimmed.vertices, dtype=np.float64)[:, :2]
            centroids = xy[dem_run.trimmed.triangles].mean(axis=1)
            catchment_ids, strays = sub_catchments.cell_ids(units.units, centroids)
            if strays:
                typer.echo(
                    f"catchments: {strays} cells had their centroid in no unit and took the "
                    "nearest",
                    err=True,
                )
            more |= {
                "catchments": escaped_ascii(units.name),
                "catchments_units": len(units.units),
                "catchments_max_shift_m": units.max_shift_m,
            }
        values |= more
    if dem_run.cut is not None:
        from tin_engine.io.mesh_index import ConformityError, MeshIndex

        assert tolerance is not None and part is not None
        folder = Path(f"{cli_kit.destination(out)}.pieces")
        folder.mkdir(exist_ok=True)
        write = functools.partial(
            _write_piece, folder, out.suffix, values, found, coded, snap_spacing, clock,
            binary, res.threads, vocabulary,
        )  # fmt: skip

        def load(window: IndexWindow) -> PieceDem:
            with clock.phase("decode"), _dem_refusals(paths):
                return load_piece(planned, plans[window], res)

        try:
            done = dem_run.cut(load, write)
        except RunError as exc:
            raise typer.BadParameter(str(exc), param_hint="--dem") from exc
        except ConformityError as exc:  # two pieces wrote one seam differently
            typer.echo(f"Error: {exc}; the piece files are written, and no index", err=True)
            raise typer.Exit(1) from exc
        credited = escaped_ascii(", ".join(names))  # the credited tile names
        if cached is not None and entry is not None:
            credited = escaped_ascii(f"{cached.source}, {entry.credit}; {credited}")
        index = MeshIndex(
            crs=dem_crs, tolerance=tolerance, source=credited,
            vocabulary=_piece_vocabulary(vocabulary).fingerprint(),
            partition=done.record, pieces=tuple(done.entries),
        )  # fmt: skip
        part_file = folder / "index.json.part"  # written last, whole or not at all
        part_file.write_text(index.model_dump_json(indent=1), encoding="utf-8")
        part_file.replace(folder / "index.json")
        typer.echo(f"{folder / 'index.json'}")
        typer.echo(f"{len(index.pieces)} pieces and their index.json in {folder}", err=True)
        return
    build = refined_record if refining else stride_record
    record = build(triangles=len(surface_mesh.triangles), threads=res.threads, **values)

    vertices = surface_mesh.vertices
    surface = cli_kit.destination(out)
    targets = [surface]
    if out_edges is not None:
        # Resolved, because two spellings of one path are still one file. Both
        # writes succeed, the second overwrites the first, the command echoes
        # two paths and exits 0 -- the caller has lost the surface they asked
        # for and nothing said so.
        constraints = cli_kit.destination(out_edges, "--out-edges")
        if constraints == surface:
            raise typer.BadParameter(
                f"--out and --out-edges both resolve to {surface}; "
                "the second would overwrite the first",
                param_hint="--out-edges",
            )
        targets.append(constraints)
    report_target = _report_target(stats, targets)
    record_target = _record_target(record_path, [*targets, report_target])
    typer.echo(summary(record), err=True)
    if unstated:
        typer.echo(
            "DEM heights read as metres; the file does not say. If they are feet, give "
            "--z-unit ft or us-ft.",
            err=True,
        )

    fields, comments = file_fields(record), field_comments(file_fields(record))
    encoders: list[Callable[[], bytes]] = []
    if out.suffix == ".gpkg":
        with clock.phase("write: geopackage"):
            _write_geopackage(
                surface, surface_mesh, codes, palette, cover, record, vocabulary, table, categories,
                regulated, catchment_ids,
            )  # fmt: skip
        typer.echo(f"{surface}")
    elif out.suffix == ".vtk":
        encoders = [
            lambda: write_vtk(
                vertices,
                triangles=surface_mesh.triangles,
                edges=surface_mesh.edges,
                edge_masks=surface_mesh.edge_masks,
                vocabulary=vocabulary,
                fields=fields,
                binary=binary,
                triangle_codes=codes,
                land_cover_codes=codes_text,
                class_fractions=table,
            )
        ]
    else:
        encoders = [
            lambda: write_ply(
                vertices,
                faces=surface_mesh.triangles,
                ascii=not binary,
                comments=[*comments, *([f"land_cover_codes {codes_text}"] if codes_text else [])],
                face_codes=codes,
                face_fractions=table,
            ),
            lambda: write_ply(
                vertices,
                edges=surface_mesh.edges,
                edge_properties=surface_mesh.edge_masks,
                ascii=not binary,
                comments=comments,
                vocabulary=vocabulary,
            ),
        ][: len(targets)]
    for target, encode in zip(targets[: len(encoders)], encoders, strict=True):
        with clock.phase("write: encode"):
            data = encode()
        with clock.phase("write: disk"):
            target.write_bytes(data)
        typer.echo(f"{target}")
    if stats is not None:
        _write_report(
            clock, report_target, surface_mesh, dem_run, targets, seams, record, res.threads
        )
    if record_target is not None:  # after the mesh and the report
        record_target.write_text(as_json(record, installed_version(), _command()), "ascii")
        typer.echo(f"{record_target}")


def _write_geopackage(
    target: Path,
    mesh: Trimmed,
    codes: npt.NDArray[np.int32] | None,
    palette: Palette | None,
    cover: list[tuple[int, BaseGeometry]],
    record: RunRecord,
    vocabulary: EdgeVocabulary,
    fractions: ClassFractions | None = None,
    categories: Mapping[int, Category] | None = None,
    regulated: npt.NDArray[np.float64] | None = None,
    catchment_ids: npt.NDArray[np.int64] | None = None,
) -> None:
    """``target`` as a GeoPackage, written as ``<target>.part``
    and renamed onto it: an existing file is replaced, never appended to, and
    a failed run leaves nothing under the asked name. Class ``fractions`` add
    their long-form table, and with ``categories`` the ``cells`` layer, its lake
    share split by ``regulated``, each triangle's share under reservoirs, and
    each cell's ``catchment_ids`` (1 when None)."""
    srs = srs_entry(next((e.value for e in record.entries if e.name == "crs"), None))
    part = target.with_name(target.name + ".part")
    part.unlink(missing_ok=True)  # a stale .part would be appended to, not replaced
    conn = create_sqlite(part)
    try:
        create_geopackage(conn, srs)
        write_triangles(conn, mesh.vertices, mesh.triangles, codes, palette, srs.srs_id)
        if len(mesh.edges):
            write_lines(conn, mesh.vertices, mesh.edges, mesh.edge_masks, vocabulary, srs.srs_id)
        if codes is not None:
            # A coded source may leave no polygon inside the domain: then there
            # is no land_cover layer to draw, and no style row for it either.
            if cover:
                write_land_cover(conn, cover, palette, srs.srs_id)
            write_classes(conn)
            if palette is not None:
                write_styles(conn, palette, ["triangles", *(["land_cover"] if cover else [])])
        if fractions is not None:
            write_fractions(conn, fractions)
            if categories is not None:
                cells = cell_table(mesh.vertices, mesh.triangles, fractions, categories, regulated)
                write_cells(conn, cells, srs.srs_id, catchment_ids)
        rows = [("rasputin_version", "rasputin version", installed_version(), None)]
        rows.append(("command", "command", _command(), None))
        write_run(conn, [*rows, *((e.name, e.wording, e.value, e.number) for e in record.entries)])
        conn.commit()
    except BaseException:
        conn.close()
        part.unlink(missing_ok=True)
        raise
    conn.close()
    os.replace(part, target)


@dataclass(frozen=True, slots=True)
class _Vector:
    """A vector source as checked: its file, layer, record name and
    catalogue entry, if any; ``key`` the --prefer or --lines key, and
    ``code`` the class a --prefer source replaces; ``fetch`` the fetch that
    copies a cached source for this run's area."""

    path: Path
    layer: str | None
    name: str
    source: RemoteSource | None
    key: str = ""
    code: int = 0
    fetch: str = ""


@dataclass(frozen=True, slots=True)
class _ClassSource:
    """``--land-cover`` as checked: a class raster, a catalogue entry with
    its year (None for a source without years) and cache root, or a file
    (``path``); or class-polygon files (``files``), with the CRS given for
    them; and the land-cover class file that names its classes, if any."""

    source: RemoteSource | None
    year: int | None
    root: Path | None
    path: Path | None
    classes: LandCoverClasses | None
    files: tuple[_Vector, ...] = ()
    crs: str | None = None


#: The words before a refusal's fetch command, which the MCP server's fetch tool takes
#: as arguments: the child the server runs cannot tell who reads the refusal.
RUN_FETCH = "; over MCP, call the fetch tool with this source and options as arguments; run: "


def _fetch_run(
    key: str,
    domain: Path | None,
    bbox: CrsBox | None,
    out_crs: str | None,
    cache: Path | None,
) -> str:
    """The ``rasputin fetch`` that copies `key` for this run's area into its cache."""
    area = f"--domain {shlex.quote(str(domain))}" if domain is not None else ""
    if bbox is not None:
        area = _bbox_args(bbox.box, bbox.crs)
    frame = f"--out-crs {shlex.quote(out_crs)}" if out_crs else ""
    store = f"--cache {shlex.quote(str(cache))}" if cache else ""
    flags = " ".join(filter(None, (area, frame, store)))
    return f"rasputin fetch {shlex.quote(key)} {flags}".rstrip()


def _bbox_fetch(key: str, exc: NotCached, cache: Path | None) -> str:
    """The refusal's tail naming the ``rasputin fetch`` of `exc`'s window
    box, in its CRS; empty when it has no box."""
    if exc.box is None or exc.crs is None:
        return ""
    given = f" --cache {shlex.quote(str(cache))}" if cache else ""
    return f"{RUN_FETCH}rasputin fetch {shlex.quote(key)} {_bbox_args(exc.box, exc.crs)}{given}"


def _bbox_args(box: tuple[float, ...], crs: str) -> str:
    """``--bbox X0 Y0 X1 Y1 --bbox-crs C``, each number as it reads back."""
    numbers = " ".join(repr(v).removesuffix(".0") for v in box)
    return f"--bbox {numbers} --bbox-crs {shlex.quote(crs)}"


def _crs_box(bbox: tuple[float, float, float, float] | None, crs: str | None) -> CrsBox | None:
    """``--bbox`` with its ``--bbox-crs``; each refusal names its option."""
    if bbox is None:
        if crs is not None:
            raise typer.BadParameter("applies only with --bbox", param_hint="--bbox-crs")
        return None
    if crs is None:
        raise typer.BadParameter(NO_BOX_CRS, param_hint="--bbox-crs")
    try:
        parse_crs(crs)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--bbox-crs") from exc
    try:
        return CrsBox(crs=crs, box=bbox)
    except ValidationError as exc:
        raise typer.BadParameter(cli_kit.words(exc), param_hint="--bbox") from exc


def _ring_vertices(polygon: Polygon) -> int:
    """The vertices of `polygon`'s rings, each closing repeat left out."""
    return sum(len(r.coords) - 1 for r in (polygon.exterior, *polygon.interiors))


def _cached_features(source: RemoteSource, root: Path | None, run: str, flag: str) -> Path:
    """A catalogue source of features' file: its own path, else its copy in the
    cache under `root`; an uncached copy's refusal names `run`, the fetch."""
    if source.path is not None:
        return source.path
    if root is None:
        raise typer.BadParameter(NO_CACHE, param_hint=flag)
    path = root / source.id / "features.gpkg"
    if not path.is_file():
        raise typer.BadParameter(
            f"{source.id}/features.gpkg is not in the cache{RUN_FETCH}{run}", param_hint=flag
        )
    return path


def _reservoir_source(
    value: str | None,
    wanted: bool,
    known: Mapping[str, RemoteSource],
    root: Path | None,
    run: str,
    flag: str = "--reservoirs",
    why: str | None = None,
) -> tuple[Path, str | None, str, tuple[str, Callable[[list[object]], list[int]]] | None] | None:
    """``flag`` (``--reservoirs`` or ``--catchments``) checked before any work:
    the file, layer and record name of a cached catalogue source of features
    (its key), or of a file rasputin reads (its name), and for an entry with
    `named_by` that rule and the names its source gives, read whole.
    ``wanted`` is whether the run writes cells: --land-cover with a land-cover
    class file, its shares and a .gpkg output; ``why``, when given, is the
    option a run that writes none lacks. An uncached source's refusal names `run`."""
    if value is None:
        return None
    _refuse_named(known.get(value), value, flag)
    if not wanted and why is not None:
        raise typer.BadParameter(
            f"needs {why}: the cells come from the class shares, which this run does not compute",
            param_hint=flag,
        )
    if not wanted:
        raise typer.BadParameter(
            "needs --land-cover with a land-cover class file and a .gpkg output", param_hint=flag
        )
    found = _features_source(value, flag, known, root, run)
    rule = known[value].named_by if value in known else None
    if rule is None:
        return (*found, None)
    names = known.get(rule.source)
    if names is None or names.path is None or names.path.suffix.lower() != ".geojson":
        raise typer.BadParameter(
            f"{value} keeps the features {rule.source} names, and {rule.source} is not a "
            "local .geojson source of features in the catalogue",
            param_hint=flag,
        )
    try:  # a GeoJSON file is read whole: the box is not asked for
        read = read_source(names.path, None, rule.field, lambda _: (0.0, 0.0, 0.0, 0.0))
    except FeatureError as exc:
        raise typer.BadParameter(str(exc), param_hint=flag) from exc
    names_given = [v for *_, v in read.rows]
    return (
        *found,
        (rule.key, functools.partial(kept_named, names=names_given, prefix=rule.prefix)),
    )


@dataclass(frozen=True)
class Units:
    """The sub-catchment units as (id, polygon), the domain they join into,
    the farthest their borders moved when simplified, and the source's name."""

    units: list[tuple[int, BaseGeometry]]
    domain: DomainPolygon
    max_shift_m: float
    name: str

    def moved(self, crs: str) -> Units:
        """The units in `crs`, each vertex moved as the domain's are."""
        if same_crs(self.domain.crs, crs):
            return self
        move = reprojector(self.domain.crs, crs)
        units = [(i, shapely.transform(g, move)) for i, g in self.units]
        return replace(self, units=units, domain=self.domain.to_crs(crs))


def _units(
    path: Path,
    layer: str | None,
    name: str,
    given: DomainPolygon,
    triangles: int | None,
    tolerance: float | None,
) -> Units:
    """The units of `path` whose inner point lies in `given`, repaired, their
    borders simplified within `tolerance`, else the band `triangles` sets, in
    their own CRS; said on stderr. A refusal is a ValueError."""
    read = read_source(path, layer, None, lambda own: given.to_crs(own).polygon.bounds)
    if parse_crs(read.crs).is_geographic:
        raise ValueError(
            f"{name} is in degrees ({crs_label(read.crs)}); the border shift is in metres, so "
            "give the units in a projected CRS"
        )
    rows = [
        (fid, shapely.force_2d(g))
        for fid, g, _ in read.rows
        if g is not None and not g.is_empty and g.geom_type in ("Polygon", "MultiPolygon")
    ]
    inside = given.to_crs(read.crs).polygon
    picked = [rows[k] for k in sub_catchments.select([g for _, g in rows], inside)]
    if not picked:
        raise ValueError(f"no unit of {name} has its inner point inside the domain")
    ids = sub_catchments.check_ids([fid for fid, _ in picked], name)
    repaired = sub_catchments.repair([g for _, g in picked])
    area = join(repaired, read.crs).polygon.area
    if tolerance is not None:
        band, why = tolerance, "--catchments-max-shift"
    else:
        assert triangles is not None
        band = sub_catchments.band_m(area, triangles)
        why = f"a triangle's edge at {triangles} triangles, {8 * band:.0f} m, over 8"
    kept = sub_catchments.simplify(repaired, band) if band > 0 else repaired
    union = orient(shapely.coverage_union_all(kept), sign=1.0)
    vertices = int(shapely.get_num_coordinates(np.asarray(kept, dtype=object)).sum())
    typer.echo(
        f"catchments: {len(kept)} units of {name}, borders simplified at {band:.0f} m ({why}): "
        f"{vertices} vertices, {area / 1e6:.2f} km2",
        err=True,
    )
    domain = DomainPolygon(polygon=union, crs=read.crs)
    return Units(list(zip(ids, kept, strict=True)), domain, band, name)


def _refuse_named(source: RemoteSource | None, value: str, flag: str) -> None:
    """Only --reservoirs reads an entry that keeps what another source names."""
    if source is not None and source.named_by is not None and flag != "--reservoirs":
        raise typer.BadParameter(
            f"{value} keeps only the features another source names; only --reservoirs reads it",
            param_hint=flag,
        )


def _features_source(
    value: str,
    flag: str,
    known: Mapping[str, RemoteSource],
    root: Path | None,
    run: str,
    suffixes: tuple[str, ...] = SUFFIXES,
) -> tuple[Path, str | None, str]:
    """`flag`'s value: the file, layer and record name of a cached catalogue
    source of features (its key), or of a file with one of `suffixes` (its
    name). An uncached source's refusal names `run`, the fetch that copies it."""
    source = known.get(value)
    if source is not None and source.values == "features":
        _refuse_named(source, value, flag)
        if Path(value).is_file():
            raise typer.BadParameter(
                f"{value} is both a catalogue source and a file here ({Path(value).resolve()}); "
                "rename or move the file, or give its path as ./" + value,
                param_hint=flag,
            )
        return _cached_features(source, root, run, flag), source.layer, value
    path = Path(value)
    if path.suffix.lower() not in suffixes or not path.is_file():
        raise typer.BadParameter(
            f"{value} is neither a catalogue source of features nor a {', '.join(suffixes)} file",
            param_hint=flag,
        )
    return path, None, path.name


def _class_source(
    values: Sequence[str],
    year: int | None,
    classes_value: str | None,
    layer: str | None,
    crs: str | None,
    root: Path | None,
    known: Mapping[str, RemoteSource],
    run: Callable[[str], str],
) -> _ClassSource | None:
    """``--land-cover`` and the options that qualify it checked before any
    work: one class raster (a catalogue key of classes, or a file whose header
    is read here), or class-polygon files and catalogue keys under one
    land-cover class file. Every refusal is a usage error naming the option;
    an uncached source's names ``run(value)``, the fetch that copies it."""
    if not values:
        given = ((year, "--land-cover-year"), (classes_value, "--land-cover-classes"))
        for value, flag in (*given, (layer, "--land-cover-layer"), (crs, "--land-cover-crs")):
            if value is not None:
                raise typer.BadParameter("needs --land-cover", param_hint=flag)
        return None
    rasters, files, names = [], [], []
    for value in values:
        source = known.get(value)
        # Features with no class file are no classes: refused as a raster below.
        if source is not None and source.values == "features" and (source.classes or classes_value):
            for given_option, flag in ((layer, "--land-cover-layer"), (crs, "--land-cover-crs")):
                if given_option is not None:
                    raise typer.BadParameter(
                        f"applies to a --land-cover file, not the catalogue source {value}",
                        param_hint=flag,
                    )
            path, own, name = _features_source(value, "--land-cover", known, root, run(value))
            files.append(_Vector(path, own, name, source))
            names += [source.classes] if source.classes else []
        elif source is not None or Path(value).suffix.lower() not in SUFFIXES:
            rasters.append((value, source))
        else:
            path = Path(value)
            if not path.is_file():
                raise typer.BadParameter(
                    f"{value} is neither a catalogue source ({', '.join(known)}) nor a file",
                    param_hint="--land-cover",
                )
            if layer is not None and path.suffix.lower() != ".gpkg":
                raise typer.BadParameter(
                    f"applies to a .gpkg only, not {path.name}", param_hint="--land-cover-layer"
                )
            files.append(_Vector(path, layer, path.name, None))
    if rasters and len(values) > 1:
        raise typer.BadParameter(
            f"{rasters[0][0]} is a class raster, which is the only --land-cover of a run",
            param_hint="--land-cover",
        )
    if len(set(names)) > 1 or (classes_value and names and {classes_value} != set(names)):
        shown = ", ".join(dict.fromkeys([*([classes_value] if classes_value else []), *names]))
        raise typer.BadParameter(
            f"the land-cover sources' class files differ ({shown}); one map has one",
            param_hint="--land-cover",
        )
    if files:
        if year is not None:
            raise typer.BadParameter("class polygons have no years", param_hint="--land-cover-year")
        name = classes_value or (names[0] if names else "")
        if not name:
            raise typer.BadParameter(
                "class polygons need a land-cover class file", param_hint="--land-cover-classes"
            )
        return _ClassSource(None, None, None, None, _load_classes(name), tuple(files), crs)
    value, source = rasters[0]
    for given_option, flag in ((layer, "--land-cover-layer"), (crs, "--land-cover-crs")):
        if given_option is not None:
            raise typer.BadParameter(f"applies to class polygons, not {value}", param_hint=flag)
    file: Path | None = None
    chosen = None
    if source is not None:
        try:
            check_values(source, "classes", "--land-cover")
        except ValueError as exc:
            raise typer.BadParameter(str(exc), param_hint="--land-cover") from exc
        chosen = year if year is not None else (source.years[-1] if source.years else None)
        try:
            source.url_for(chosen)
        except ValueError as exc:
            raise typer.BadParameter(str(exc), param_hint="--land-cover-year") from exc
        if root is None:
            raise typer.BadParameter(NO_CACHE, param_hint="--land-cover")
    else:
        if year is not None:
            raise typer.BadParameter("a file has no years", param_hint="--land-cover-year")
        file = Path(value)
        if not file.is_file():
            raise typer.BadParameter(
                f"{value} is neither a catalogue source ({', '.join(known)}) nor a file",
                param_hint="--land-cover",
            )
        try:  # the header and one cell: a class raster in uint8 with a CRS
            corner = IndexWindow(row0=0, col0=0, rows=1, cols=1)
            TiffDemRepository([file], z_unit="metre").load_class_window(file.name, corner)
        except (OSError, ValueError) as exc:
            raise typer.BadParameter(str(exc), param_hint="--land-cover") from exc
    name = classes_value or (source.classes if source is not None else "")
    found = _load_classes(name) if name else None
    return _ClassSource(source, chosen, None if file else root, file, found)


def _load_classes(name: str) -> LandCoverClasses:
    try:
        return load_classes(name)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--land-cover-classes") from exc


def _line_sources(
    values: Sequence[str],
    known: Mapping[str, RemoteSource],
    root: Path | None,
    run: Callable[[str], str],
) -> list[_Vector]:
    """Each ``--lines KEY=SOURCE`` checked before any work: a line kind,
    one source per kind, a file or cached source of features."""
    found: list[_Vector] = []
    for value in values:
        key, _, given = value.partition("=")
        if not given:
            raise typer.BadParameter(
                f"{value} names no source; give KEY=SOURCE", param_hint="--lines"
            )
        if key not in LINE_KEYS:
            area = key in get_args(Category)
            why = (
                "an area kind; give it to --prefer"
                if area
                else f"no line kind; use {', '.join(LINE_KEYS)}"
            )
            raise typer.BadParameter(f"{key} is {why}", param_hint="--lines")
        if any(s.key == key for s in found):
            raise typer.BadParameter(
                f"{key} is given twice; one source per kind", param_hint="--lines"
            )
        path, layer, name = _features_source(given, "--lines", known, root, run(given))
        found.append(_Vector(path, layer, name, known.get(given), key))
    return found


def _edge_choice(
    selectors: Sequence[str], classes: _ClassSource | None, keys: Sequence[str]
) -> EdgeChoice:
    """``--edges`` resolved: ``all``, ``water``, a category of the class
    file, a class code, a ``--lines`` key; refused naming ``--edges``."""
    lcc = classes.classes if classes is not None and classes.files else None
    categories = sorted(set(lcc.categories.values())) if lcc is not None else []
    every, codes, lines = False, set[int](), set[str]()
    for s in selectors:
        if s == "all":
            every, lines = lcc is not None, lines | set(keys)
            continue
        if s in keys:
            lines.add(s)
            continue
        if s in LINE_KEYS:
            raise typer.BadParameter(
                f"{s} needs its lines: add --lines {s}=SOURCE", param_hint="--edges"
            )
        known = known_codes(lcc) - {0} if lcc is not None else frozenset[int]()
        area = s == "water" or s in get_args(Category) or s.isdigit()
        if not area or (lcc is not None and s.isdigit() and int(s) not in known):
            accepted = ", ".join(["all", "water", *categories, *keys])
            raise typer.BadParameter(
                f"{s} is no selector; use {accepted} or a class code", param_hint="--edges"
            )
        if lcc is None:
            raise typer.BadParameter(
                f"{s} selects classes, and no class polygons come from --land-cover",
                param_hint="--edges",
            )
        if s == "water":
            codes |= lcc.water_codes
        elif s.isdigit():
            codes.add(int(s))
        else:
            codes |= {c for c, named in lcc.categories.items() if named == s}
    return EdgeChoice(every_class=every, codes=frozenset(codes), lines=frozenset(lines))


def _class_fractions(
    classes: _ClassSource,
    land: _LandCover | None,
    trimmed: Trimmed,
    dem_crs: str,
    clock: PhaseClock,
    domain: Path | None,
    out_crs: str | None,
    cache: Path | None,
    threads: int = 1,
) -> tuple[ClassFractions, npt.NDArray[np.int32], dict[str, Value]]:
    """The exact class fractions from the cached year, the class raster file
    or the map of class polygons (``land``), timed as the phase ``land cover
    fractions``; the table, the codes and the record entries. An object the
    cache lacks is refused naming the ``fetch`` that gets it."""
    source, year = classes.source, classes.year
    repository: ClassRaster | None = None
    values: dict[str, Value] = {}
    if land is not None:
        label = "; ".join(f.source.id if f.source else escaped_ascii(f.name) for f in classes.files)
    elif source is None:
        assert classes.path is not None
        name = classes.path.name
        label = escaped_ascii(name)
        repository = TiffDemRepository([classes.path], z_unit="metre", threads=threads)
    else:
        assert classes.root is not None
        cached_classes = CacheRepository(
            classes.root, source.id, nodata=source.nodata, threads=threads, z_unit="metre"
        )
        repository = cached_classes
        url = source.url_for(year)
        name = url.rsplit("/", 1)[-1].split("?")[0].rsplit(".", 1)[0]
        label = source.id
        try:
            cached = {f.name for f in cached_classes.footprints()}
        except NotCached:
            cached = set()
        if name not in cached:
            words = [
                f"rasputin fetch {source.id}",
                f"--year {year}" if year is not None else "",
                f"--domain {shlex.quote(str(domain))}" if domain is not None else "",
                f"--out-crs {shlex.quote(out_crs)}" if out_crs else "",
                f"--cache {shlex.quote(str(cache))}" if cache else "",
            ]
            run = " ".join(filter(None, words))
            raise typer.BadParameter(
                f"{source.id}/{name} is not in the cache{RUN_FETCH}{run}",
                param_hint="--land-cover",
            )
        fetched = cached_classes.fetched(name) or datetime.date.today()
    with clock.phase("land cover fractions"):
        try:
            if repository is None:
                assert land is not None and classes.classes is not None
                xy = np.asarray(trimmed.vertices, dtype=np.float64)[:, :2]
                listed = known_codes(classes.classes)
                found, repaired = list(land.map.areas), land.base.repaired
                assigned = assign_polygon_fractions(xy, trimmed.triangles, found, known=listed)
            else:
                assigned = assign_fractions(
                    trimmed.vertices,
                    trimmed.triangles,
                    mesh_crs=dem_crs,
                    repository=repository,
                    name=name,
                    threads=threads,
                )
        except (ValueError, NotCached) as exc:
            raise typer.BadParameter(cli_kit.words(exc), param_hint="--land-cover") from exc
    typer.echo(
        f"land cover: {len(assigned.table.classes)} class shares over "
        f"{assigned.table.triangles} triangles, {len(np.unique(assigned.table.classes))} classes",
        err=True,
    )
    values |= {
        "land_cover_codes": f"{label} class number; 0 = no data, and every constraint line",
        "land_cover_source": label if year is None else f"{label}, {year}",
    }
    if repository is None:
        values["land_cover_repaired"] = repaired
        credits = [_credits(f.source, f.path) for f in classes.files if f.source is not None]
        if credits:  # a file carries no credit, licence or citation
            values |= {
                "land_cover_credit": "; ".join(c for c, _, _ in credits),
                "land_cover_licence_note": "; ".join(n for _, n, _ in credits),
                "land_cover_cite": "; ".join(c for *_, c in credits if c) or None,
            }
    elif source is not None:
        values |= {
            "land_cover_credit": escaped_ascii(source.credit_for(url=url, date=fetched)),
            "land_cover_licence_note": escaped_ascii(source.licence_note),
            "land_cover_cite": escaped_ascii("; ".join(source.cite)) or None,
        }
    return assigned.table, assigned.codes, values


def _credits(source: RemoteSource, path: Path) -> tuple[str, str, str]:
    """The credit, licence note and citation of polygons read from ``path``."""
    url, fetched = _fetched(source, path)
    return (
        escaped_ascii(source.credit_for(url=url, date=fetched)),
        escaped_ascii(source.licence_note),
        escaped_ascii("; ".join(source.cite)),
    )


def _prefer_sources(
    values: Sequence[str],
    classes: _ClassSource | None,
    known: Mapping[str, RemoteSource],
    root: Path | None,
    run: Callable[[str], str],
) -> list[_Vector]:
    """Each ``--prefer`` checked before any work: class polygons from
    --land-cover, a key its class file's ``[cells.prefer]`` names, and a
    source of polygons. An uncached source's refusal names ``run(source)``."""
    if values and classes is None:
        raise typer.BadParameter("needs --land-cover", param_hint="--prefer")
    found = []
    for value in values:
        assert classes is not None
        key, _, given = value.partition("=")
        if key in LINE_KEYS:
            raise typer.BadParameter(
                f"{key} is a line kind; give it to --lines", param_hint="--prefer"
            )
        if classes.classes is None or not classes.files:
            raise typer.BadParameter(
                "needs class polygons from --land-cover; a raster --land-cover has none",
                param_hint="--prefer",
            )
        prefer, name = classes.classes.prefer, classes.classes.name
        if not given:
            raise typer.BadParameter(
                f"{value} names no source; give KEY=SOURCE", param_hint="--prefer"
            )
        if not prefer:
            raise typer.BadParameter(
                f"land-cover class file {name} names no class a preferred source replaces",
                param_hint="--prefer",
            )
        if key not in prefer:
            raise typer.BadParameter(
                f"{key} is no key of land-cover class file {name}; use {', '.join(prefer)}",
                param_hint="--prefer",
            )
        path, layer, shown = _features_source(given, "--prefer", known, root, run(given))
        code = next(c for k, c in prefer.items() if k == key)
        source = known.get(given)
        fetch = run(given) if source is not None and source.path is None else ""
        found.append(_Vector(path, layer, shown, source, key, code, fetch))
    return found


def _fetched(source: RemoteSource | None, path: Path) -> tuple[str, datetime.date]:
    """The URL and date to credit for polygons read from ``path``: what was
    copied, for a fetched catalogue source; else the path and today."""
    url, fetched = str(path), datetime.date.today()
    if source is not None and source.path is None:
        copied = CacheWriter(path.parent.parent, source.id).feature_manifest()
        url, fetched = (copied.url, copied.first_fetched) if copied else (url, fetched)
    return url, fetched


@dataclass(frozen=True, slots=True)
class _LandCover:
    """The map of class polygons, the map before ``--prefer`` (the same
    without it), and what each ``--land-cover`` file read."""

    map: LandCoverMap
    base: LandCoverMap
    reads: tuple[ClassPolygons, ...]


def _land_cover_map(
    classes: _ClassSource,
    preferred: Sequence[_Vector],
    area: DomainPolygon,
    dem_crs: str,
    cell: float,
    clock: PhaseClock,
    threads: int,
) -> _LandCover:
    """The --land-cover files read near ``area`` and joined, cut by the
    --prefer sources: the phase ``land cover read``."""
    assert classes.classes is not None
    with clock.phase("land cover read"):
        reads = []
        for f in classes.files:
            try:
                reads.append(
                    read_class_polygons(
                        f.path,
                        f.layer,
                        classes.classes,
                        dem_crs,
                        lambda own: source_region(area, dem_crs, own, cell=cell).bounds,
                        crs=classes.crs,
                        threads=threads,
                    )
                )
            except FeatureError as exc:
                hint = "--land-cover-crs" if "given CRS" in str(exc) else "--land-cover"
                raise typer.BadParameter(str(exc), param_hint=hint) from exc
        xy = shapely.get_coordinates(area.polygon.exterior)
        region = source_region(area, dem_crs, dem_crs, cell=cell)
        chosen = []
        for p in preferred:
            try:
                polygons, _ = read_polygons(
                    p.path, p.layer, dem_crs, xy,
                    functools.partial(_say_repaired, f"--prefer {p.key}"),
                    lines="lines go to --lines",
                )  # fmt: skip
            except FeatureError as exc:
                raise typer.BadParameter(str(exc), param_hint="--prefer") from exc
            if not np.any(shapely.intersects(polygons, region)):
                box = " ".join(f"{v:.0f}" for v in region.bounds)
                why = f"; the cache holds only what was fetched{RUN_FETCH}{p.fetch}"
                raise typer.BadParameter(
                    f"{p.name} has no polygon meeting the read region, box {box} in {dem_crs}"
                    + (why if p.fetch else f"; drop --prefer {p.key} if it has none"),
                    param_hint="--prefer",
                )
            chosen.append((p.code, polygons))
        try:
            base = build([r.areas for r in reads], [], classes.classes, region, threads=threads)
        except FeatureError as exc:
            raise typer.BadParameter(str(exc), param_hint="--land-cover") from exc
        again = functools.partial(build, fix=False, outside=base.outside)
        built = again([base.areas], chosen, classes.classes, region) if chosen else base
        _say_repaired("land cover", base.repaired, sum(len(r.areas) for r in reads))
    for r in reads:
        if r.scanned:
            typer.echo(
                f"{r.layer}: this layer has no spatial index, so every row was read", err=True
            )
    return _LandCover(built, base, tuple(reads))


def _read_lines(
    sources: Sequence[_Vector],
    area: DomainPolygon,
    dem_crs: str,
    cell: float,
    clock: PhaseClock,
    threads: int,
) -> tuple[list[tuple[str, LineString]], list[SourceLines]]:
    """The --lines files read near ``area``: the phase ``lines read``."""
    pairs: list[tuple[str, LineString]] = []
    reads: list[SourceLines] = []
    if not sources:
        return pairs, reads
    with clock.phase("lines read"):
        for s in sources:
            try:
                read = read_linestrings(
                    s.path, s.layer, dem_crs,
                    lambda own: source_region(area, dem_crs, own, cell=cell).bounds,
                    threads=threads,
                )  # fmt: skip
            except FeatureError as exc:
                raise typer.BadParameter(str(exc), param_hint="--lines") from exc
            reads.append(read)
            pairs += [(s.key, line) for line in read.lines]
    for r in reads:
        if r.scanned:
            typer.echo(
                f"{r.layer}: this layer has no spatial index, so every row was read", err=True
            )
    return pairs, reads


def _edges(
    land: _LandCover | None,
    pairs: Sequence[tuple[str, LineString]],
    choice: EdgeChoice,
    cleanup: EdgeCleanup,
    domain: DomainPolygon,
    vocabulary: EdgeVocabulary,
    region: Polygon,
    clock: PhaseClock,
    empty: int,
) -> FeatureSet:
    """`edges_of`, timed as ``edges`` and ``edges: clean-up``, with its stderr
    line; ``empty`` counts the rows the readers skipped as empty."""
    t0 = time.perf_counter()
    try:
        found = edges_of(
            None if land is None else land.map, pairs, choice, cleanup, domain, vocabulary,
            region=region,
        )  # fmt: skip
    except FeatureError as exc:
        raise typer.BadParameter(str(exc), param_hint="--edges") from exc
    clock.add("edges", time.perf_counter() - t0)
    if found.cover_vertices is not None:
        clock.add("edges: clean-up", found.cleanup_seconds)
    empty += found.empty + (sum(r.empty for r in land.reads) if land is not None else 0)
    typer.echo(
        f"edges: {len(found.features)} kept ({found.clipped} cut at the domain outline), "
        f"{found.outside} outside the domain, {empty} empty",
        err=True,
    )
    return found


def _edge_record(
    found: FeatureSet,
    land: _LandCover | None,
    files: Sequence[_Vector],
    sources: Sequence[_Vector],
    pairs: Sequence[tuple[str, LineString]],
    dem_crs: str,
) -> dict[str, Value]:
    """The record's entries for the edges: one ``edges`` entry per source
    of edges (`found.counts`), the land-cover ``files``' CRS and credit, the
    line kinds' lines, and the clean-up's vertices."""
    texts: list[str] = []
    entries: list[tuple[str, str]] = []
    if land is not None:
        layers = ", ".join(dict.fromkeys(r.layer for r in land.reads if r.layer))
        names = " + ".join(f.name for f in files) + (f" layer {layers}" if layers else "")
        entries.append((names, f"classes {land.map.classes.name}"))
    by_key = {s.key: s for s in sources}
    entries += [(by_key[k].name, f"lines {k}") for k in dict.fromkeys(k for k, _ in pairs)]
    start = 0
    for (what, of), n in zip(entries, found.counts, strict=True):
        kept = [g for f in found.features[start : start + n] for g in f.lines]
        start += n
        coords = sum(int(shapely.get_num_coordinates(g)) for g in kept)
        mine = f"{plural(len(kept), 'line', 'lines')}, {plural(coords, 'vertex', 'vertices')}"
        texts.append(f"{what}, {of}: {plural(n, 'feature', 'features')}, {mine}")
    values: dict[str, Value] = {"edges": escaped_ascii("; ".join(texts))}
    rows = []
    for s in sources:
        lines = [g for f in found.features if str(f.fid).startswith(f"{s.key} ") for g in f.lines]
        counted = plural(len(lines), "line", "lines")
        km = sum(g.length for g in lines) / 1e3
        rows.append(f"{s.key} from {escaped_ascii(s.name)}: {counted}, {km:.3f} km")
    values["lines"] = "; ".join(rows) or None
    if land is not None:
        lcc = land.map.classes
        values["land_cover_crs"] = "; ".join(crs_label(r.crs) for r in land.reads)
        values["land_cover_transform"] = "; ".join(
            transform_label(r.crs, dem_crs) for r in land.reads
        )
        values["land_cover_notice"] = lcc.notice or None
        values["land_cover_licence"] = escaped_ascii(lcc.licence) or None
    if found.cover_vertices is not None:
        before, after = found.cover_vertices
        values["land_cover_vertices"] = f"{before} after the clip, {after} after clean-up"
        values["land_cover_area_moved_m2"] = f"{found.area_changed:.1f}"
    return values


def _overlap_record(land: _LandCover, files: Sequence[_Vector], area: DomainPolygon) -> str:
    """Per file that kept an overlap, the area inside ``area`` it kept, m2."""
    shapely.prepare(area.polygon)
    kept: dict[int, float] = {}
    for i, g in land.base.overlap:
        kept[i] = kept.get(i, 0.0) + float(shapely.intersection(g, area.polygon).area)
    words = [f"{a:.1f} m2 kept by {escaped_ascii(files[i].name)}" for i, a in kept.items() if a]
    return "; ".join(words) or "none"


def _prefer_record(
    land: _LandCover, preferred: Sequence[_Vector], trimmed: Trimmed, classes: _ClassSource | None
) -> dict[str, Value]:
    """The record entries of the --prefer sources, their areas inside the
    mesh, and a warning per key whose source gives less than the land cover did."""
    assert land.map.replaced is not None and classes is not None
    label = " + ".join(escaped_ascii(f.name) for f in classes.files)
    xy = np.asarray(trimmed.vertices, dtype=np.float64)[:, :2]
    codes = list(dict.fromkeys(p.code for p in preferred))
    base = [(c, g) for c, g in land.base.areas if c in codes]
    areas = areas_inside(base, land.map.replaced, codes, xy[trimmed.triangles])
    rows, credits, licences, cites = [], [], [], []
    for key in dict.fromkeys(p.key for p in preferred):
        mine = [p for p in preferred if p.key == key]
        names = " and ".join(escaped_ascii(p.name) for p in mine)
        a = areas[mine[0].code]
        given, replaced, refilled = (x / 1e6 for x in (a.given, a.replaced, a.refilled))
        rows.append(
            f"{key} from {names}: {replaced:.4f} km2 replaced, {given:.4f} km2 given, "
            f"{refilled:.4f} km2 refilled"
        )
        if a.given < a.replaced * (1 - 1e-9):
            what = "water" if key == "lake" else key
            typer.echo(
                f"{key}: {names} gives {given:.4f} km2 where {label} gave {replaced:.4f} km2; "
                f"{refilled:.4f} km2 of {label} {what} had no {names} {key} and took a "
                "neighbouring class",
                err=True,
            )
        for p in mine:
            if p.source is None:
                credits.append(f"{key}: {escaped_ascii(p.name)}")
                continue
            credit, licence, _ = _credits(p.source, p.path)
            credits.append(f"{key}: {credit}")
            licences.append(f"{key}: {licence}")
            cites += [f"{key}: {escaped_ascii(c)}" for c in p.source.cite]
    return {
        "land_cover_prefer": "; ".join(rows),
        "land_cover_prefer_credit": "; ".join(credits),
        "land_cover_prefer_licence_note": "; ".join(licences) or None,
        "land_cover_prefer_cite": "; ".join(cites) or None,
    }


def _say_repaired(what: str, repaired: int, read: int) -> None:
    """One stderr line when the polygons read had self-crossing rings."""
    if repaired:
        typer.echo(
            f"{what}: {repaired} of {read} polygons crossed itself and was repaired", err=True
        )


def _piece_vocabulary(vocabulary: EdgeVocabulary) -> EdgeVocabulary:
    """The run's vocabulary plus `seam`, for a cut run's piece files."""
    seam = [p for p in PIECE_VOCABULARY.properties if p.name == "seam"]
    return EdgeVocabulary(properties=(*vocabulary.properties, *seam))


def _write_piece(
    where: Path,
    suffix: str,
    values: Mapping[str, Value],
    found: FeatureSet | None,
    coded: LandCoverClasses | None,
    spacing: float,
    clock: PhaseClock,
    binary: bool,
    threads: int,
    vocabulary: EdgeVocabulary,
    piece: Piece,
) -> PieceEntry:
    """``piece``, if it has triangles, as its own ``suffix`` file in ``where``,
    its record the run's ``values`` and its own; its entry."""
    from tin_engine.io.mesh_index import PieceCounts, PieceEntry

    file, digest, made = None, None, piece.trimmed
    if len(made.triangles):
        codes, codes_text, figures = None, "", dict(values) | dict(piece.values)
        if coded is not None and found is not None:
            codes, codes_text = _land_cover(made, found, coded, spacing, clock)
            figures["land_cover_codes"] = codes_text
        told = refined_record(triangles=len(made.triangles), threads=threads, **figures)
        typer.echo(f"piece {'-'.join(map(str, piece.id))}: {summary(told)}", err=True)
        if suffix == ".vtk":
            data = write_vtk(
                made.vertices, triangles=made.triangles, edges=made.edges,
                edge_masks=made.edge_masks, vocabulary=_piece_vocabulary(vocabulary),
                fields=file_fields(told), binary=binary, triangle_codes=codes,
                land_cover_codes=codes_text,
            )  # fmt: skip
        else:
            coded_text = [f"land_cover_codes {codes_text}"] if codes_text else []
            data = write_ply(
                made.vertices, faces=made.triangles, ascii=not binary,
                comments=[*field_comments(file_fields(told)), *coded_text], face_codes=codes,
            )  # fmt: skip
        file = "-".join(map(str, piece.id)) + suffix
        digest = hashlib.sha256(data).hexdigest()
        (where / file).write_bytes(data)
    counts = PieceCounts(
        triangles=len(made.triangles), vertices=len(made.vertices), on_frozen=piece.on_frozen
    )
    return PieceEntry(id=piece.id, file=file, sha256=digest, window=piece.window, counts=counts)


def _land_cover(
    trimmed: Trimmed, found: FeatureSet, rules: LandCoverClasses, spacing: float, clock: PhaseClock
) -> tuple[npt.NDArray[np.int32], str]:
    """A code per triangle from the coded polygons, timed as the phase
    ``land cover``, with its stderr line; the codes and their text."""
    polygons = [(f.polygon, f.code) for f in found.features if f.polygon is not None and f.code]
    with clock.phase("land cover"):
        mesh = (trimmed.vertices, trimmed.triangles, trimmed.edges)
        labels = label_triangles(*mesh, polygons=polygons, margin=2 * spacing)
    typer.echo(
        f"land cover: {labels.regions} areas between lines; {labels.outside} in no polygon, "
        f"{labels.overlapped} in more than one (the smallest wins), "
        f"{labels.thin} too narrow to label with certainty",
        err=True,
    )
    splits = "".join(
        f", class {s.code} split by {s.attribute} into "
        + ", ".join(map(str, sorted(s.codes.values())))
        for s in rules.splits
    )
    text = f"{rules.codes}, attribute {rules.attribute}{splits}, legend {rules.name}; "
    return labels.codes, text + "0 = in no polygon, and every constraint line"


def palette(
    name: Annotated[
        str,
        typer.Argument(help=f"A legend: {', '.join(legend_names())}, or a .toml file."),
    ],
    out: Annotated[
        Path | None,
        typer.Option("--out", writable=True, help="Write the preset here, not to stdout."),
    ] = None,
) -> None:
    """Write a ParaView colour preset for a mesh's land_cover_code.

    Import it in ParaView's Colour Map Editor (Choose Preset, Import), then
    colour by land_cover_code; tick Interpret Values As Categories if the
    preset does not.
    """
    try:
        legend = load_legend(name)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    text = json.dumps(paraview_preset(legend.colours, legend.preset), indent=1) + "\n"
    if out is None:
        typer.echo(text, nl=False)
        return
    cli_kit.destination(out).write_text(text)
    typer.echo(f"{out}")


def _report_target(stats: Path | None, meshes: list[Path]) -> Path | None:
    """Where ``--stats`` writes, None for stdout (``-``) or no report,
    checked before any mesh file is written."""
    if stats is None or stats == Path("-"):
        return None
    target = cli_kit.destination(stats, "--stats")
    if target in meshes:
        raise typer.BadParameter(
            f"resolves to {target}; the report would overwrite the mesh", param_hint="--stats"
        )
    return target


def _record_target(record: Path | None, taken: list[Path | None]) -> Path | None:
    """Where ``--record`` writes, checked before any file is written; never
    standard output (``--stats -``'s) nor a file this run also writes."""
    if record is None:
        return None
    if record == Path("-"):
        raise typer.BadParameter("standard output is --stats -'s", param_hint="--record")
    target = cli_kit.destination(record, "--record")
    if target in taken:
        raise typer.BadParameter(
            f"resolves to {target}, which this run also writes", param_hint="--record"
        )
    return target


def _command() -> str:
    """The command line as given, for ``--stats`` and ``--record``."""
    return shlex.join([Path(sys.argv[0]).name, *sys.argv[1:]])


def _write_report(
    clock: PhaseClock,
    target: Path | None,
    trimmed: Trimmed,
    dem_run: _DemMesh,
    files: list[Path],
    seams: tuple[Seam, ...],
    record: RunRecord,
    threads: int,
) -> None:
    """Build the report from what ran and write it, or print it for ``-``.
    The total stops here; the quality pass is timed on its own line.
    ``threads`` is the count the run used."""
    total = clock.elapsed()
    t0 = time.perf_counter_ns()
    meta = dem_run.meta
    sizes = Sizes(
        output_vertices=len(trimmed.vertices),
        output_triangles=len(trimmed.triangles),
        constraint_edges=len(trimmed.edges),
        files=[(f.name, f.stat().st_size) for f in files],
        dem_nodes=(meta.rows, meta.cols),
        dem_spacing=(meta.delta_x, meta.delta_y),
        domain_vertices=dem_run.domain_vertices,
        domain_holes=dem_run.domain_holes,
        start_vertices=dem_run.start_vertices,
        start_triangles=dem_run.start_triangles,
        resampled=any(e.name == "resampled_grid" for e in record.entries),
    )
    rows = [((e.wording, e.value, e.name), e.in_inputs) for e in record.entries]
    measured = quality(trimmed.vertices, trimmed.triangles)
    text = render(
        Report(
            command=_command(),
            sizes=sizes,
            quality=measured,
            phases=clock.phases(),
            total=total,
            stats_seconds=(time.perf_counter_ns() - t0) / 1e9,
            bounds_checks=bounds_checks(),
            threads=threads,
            seams=[s.cells() for s in seams],
            inputs=[row for row, inputs in rows if inputs],
            result=[row for row, inputs in rows if not inputs],
        )
    )
    if target is None:
        typer.echo(text, nl=False)
        return
    target.write_text(text, encoding="utf-8")
    typer.echo(f"{target}")


def fetch(
    source: Annotated[str, typer.Argument(help=f"A catalogue source: {', '.join(SOURCES)}.")],
    domain: Annotated[
        Path | None,
        typer.Option("--domain", help="Fetch what a mesh of this polygon reads, as mesh reads it."),
    ] = None,
    domain_crs: Annotated[
        str | None, typer.Option("--domain-crs", help="The --domain file's CRS, as for mesh.")
    ] = None,
    bbox: Annotated[
        tuple[float, float, float, float] | None,
        typer.Option("--bbox", metavar="BOX", help=f"Or {BBOX_HELP}."),
    ] = None,
    bbox_crs: BBOX_CRS = None,
    out_crs: Annotated[
        str | None,
        typer.Option(
            "--out-crs",
            help="The mesh's CRS, which --domain or --bbox is moved into, as mesh moves it.",
        ),
    ] = None,
    cache: Annotated[
        Path | None, typer.Option("--cache", writable=True, help=f"The tile cache. {CACHE_DEFAULT}")
    ] = None,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Read the headers, print the plan, write nothing.")
    ] = False,
    refresh: Annotated[
        bool, typer.Option("--refresh", help="Discard the source's cache and fetch anew.")
    ] = False,
    connections: Annotated[
        int, typer.Option("--connections", min=1, help="Parallel requests, for byte-range sources.")
    ] = 8,
    year: Annotated[
        int | None,
        typer.Option(
            "--year", help="For a source with one file per year: which. Default: the last."
        ),
    ] = None,
) -> None:
    """Copy what a mesh of a domain reads from a remote source into the tile cache.

    Then rasputin mesh --dem SOURCE runs offline. A source is the package's or
    one in your sources file (see rasputin sources). A source of map-service
    features is cached as one GeoPackage, which rasputin sources lists by path
    for the vector options.
    """
    import asyncio

    from tin_engine.fetch.features import fetch_features
    from tin_engine.fetch.http import FetchError, RangeClient
    from tin_engine.fetch.plan import FetchRequest
    from tin_engine.fetch.run import fetch as run
    from tin_engine.io.repository import CacheError, CacheWriter

    known = _known_sources(os.environ)
    if source not in known:
        raise typer.BadParameter(f"{source} is not a catalogue source ({', '.join(known)})")
    if known[source].kind == "local":
        raise typer.BadParameter(f"{source} is a local source, read in place; it needs no fetch")
    if year is not None:  # refused here, before any request
        try:
            known[source].url_for(year)
        except ValueError as exc:
            raise typer.BadParameter(str(exc), param_hint="--year") from exc
    box = _crs_box(bbox, bbox_crs)
    if (domain is None) == (bbox is None):
        raise typer.BadParameter("give exactly one of --domain and --bbox", param_hint="--bbox")
    root = cache_root(cache, os.environ)
    if root is None:
        raise typer.BadParameter(
            NO_CACHE,
            param_hint="--cache",
        )
    try:
        given = read_domain(domain, domain_crs) if domain is not None else None
        request = FetchRequest(
            source=source, domain=given, box=box, out_crs=out_crs, connections=connections,
            dry_run=dry_run, refresh=refresh, year=year,
        )  # fmt: skip
    except DomainError as exc:
        raise typer.BadParameter(str(exc), param_hint="--domain") from exc
    except ValueError as exc:
        raise typer.BadParameter(cli_kit.words(exc)) from exc
    tenths = [0]

    def progress(done: int, total: int) -> None:
        if total and done * 10 // total > tenths[0]:
            tenths[0] = done * 10 // total
            typer.echo(f"{done:,} of {total:,} bytes", err=True)

    client, writer = RangeClient(), CacheWriter(root, source)
    head = "would fetch" if dry_run else "fetched"
    try:
        if known[source].kind == "arcgis-features":
            got = fetch_features(request, known[source], client.get_text, writer)
        else:
            report = asyncio.run(run(request, known[source], client, writer, progress=progress))
    except BoxError as exc:
        raise typer.BadParameter(str(exc), param_hint="--bbox-crs") from exc
    except (FetchError, CacheError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1) from exc
    if known[source].kind == "arcgis-features":
        typer.echo(
            f"{source}: {got.in_box:,} features in the box, {got.present:,} present, {head} "
            f"{got.fetched:,} in {got.requests:,} requests, {got.bytes:,} bytes, "
            f"{got.seconds:.1f} s into {writer.directory}"
        )
        return
    for plan in report.plans:
        typer.echo(f"{plan.object_id}: {len(plan.blocks):,} blocks, {plan.bytes:,} bytes to fetch")
    typer.echo(
        f"{source}: {report.objects} objects, {report.needed:,} blocks needed, "
        f"{report.present:,} present, {head} {report.fetched:,} ({report.empty} empty) in "
        f"{report.requests:,} requests, {report.bytes:,} bytes, {report.seconds:.1f} s "
        f"into {writer.directory}"
    )
    if report.no_tile:
        typer.echo(f"no tile (sea): {', '.join(report.no_tile)}")


def fetch_stations(
    source: Annotated[str, typer.Argument(help=f"A station list: {', '.join(STATION_SOURCES)}.")],
    out_dir: Annotated[
        Path, typer.Option("--out-dir", writable=True, help="Where to write the files.")
    ],
    refresh: Annotated[
        bool, typer.Option("--refresh", help="Fetch again even if the files are there.")
    ] = False,
) -> None:
    """Fetch a list's stations, catchment polygons and nearby river lines into --out-dir.

    NOTICE.txt and a manifest are written beside them. Files already there are
    kept, and nothing is requested, unless --refresh.
    """
    from tin_engine.fetch.http import FetchError, RangeClient
    from tin_engine.fetch.nve import FILES, fetch_station_set

    if source not in STATION_SOURCES:
        raise typer.BadParameter(
            f"{source} is not a station list ({', '.join(STATION_SOURCES)})", param_hint="SOURCE"
        )
    names = (*FILES, "manifest.json")
    if not refresh and all((out_dir / name).is_file() for name in names):
        typer.echo(f"{out_dir}: already fetched; --refresh fetches again", err=True)
        return
    try:
        files = fetch_station_set(STATION_SOURCES[source], RangeClient().get_text)
    except KeyError as exc:  # a field the fetch reads, left out of the answer
        typer.echo(f"Error: the NVE service's answer is missing the field {exc.args[0]}", err=True)
        raise typer.Exit(1) from exc
    except (FetchError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1) from exc
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, data in files.items():  # the manifest last
        (out_dir / name).write_bytes(data)
    typer.echo(f"{source}: {len(files)} files written to {out_dir}", err=True)


def sources(
    bbox: Annotated[
        tuple[float, float, float, float] | None,
        typer.Option(
            "--bbox",
            metavar="BOX",
            help=f"{BBOX_HELP}: report the share of it each source covers. A point is a "
            "box with XMIN = XMAX and YMIN = YMAX.",
        ),
    ] = None,
    bbox_crs: BBOX_CRS = None,
    cache: Annotated[
        Path | None, typer.Option("--cache", help=f"The tile cache read. {CACHE_DEFAULT}")
    ] = None,
) -> None:
    """List the elevation and land-cover sources this machine knows, as JSON.

    Each, local or remote, with its coverage (a longitude/latitude box), cell
    size in metres, CRS, height unit, tile count, description and credit, read
    from the headers; and the sources file and cache read.
    """
    from tin_engine.source_report import listing

    box = _crs_box(bbox, bbox_crs)
    try:
        found = listing(os.environ, cache_root(cache, os.environ), box)
    except BoxError as exc:
        raise typer.BadParameter(str(exc), param_hint="--bbox-crs") from exc
    typer.echo(found.model_dump_json())


def gml_to_gpkg(
    source: Annotated[Path, typer.Argument(help="The GML 3.2 file to convert.")],
    out: Annotated[
        Path, typer.Option("--out", writable=True, help="The GeoPackage to write (.gpkg).")
    ],
    crs: Annotated[
        str | None, typer.Option("--crs", help="The CRS, when the geometries name none.")
    ] = None,
) -> None:
    """Convert a GML 3.2 file of surfaces to a GeoPackage, once, for --land-cover and the catalogue.

    One table of the file's feature type: a row per feature with a surface,
    its gml:id and its plain values as columns, in the file's CRS. The file
    is written as OUT.part and renamed onto OUT when finished.
    """
    if out.suffix.lower() != ".gpkg":
        raise typer.BadParameter(f"{out} must end in .gpkg", param_hint="--out")
    target = cli_kit.destination(out)
    part = target.with_name(target.name + ".part")
    part.unlink(missing_ok=True)  # a stale .part would be appended to, not replaced
    try:
        with source.open("rb") as stream, closing(create_sqlite(part)) as conn:
            done = gml_to_geopackage(stream, conn, crs=crs)
    except (OSError, ValueError) as exc:
        part.unlink(missing_ok=True)
        raise typer.BadParameter(f"{source}: {exc}") from exc
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    os.replace(part, target)
    typer.echo(
        f"{out}: {done.features} features of {done.layer} in {crs_label(done.crs)},"
        f" {done.without_surface} without a surface; columns {', '.join(done.columns)}",
        err=True,
    )
    typer.echo(f"{out}")


def cache_root(option: Path | None, environ: Mapping[str, str]) -> Path | None:
    """The tile cache: ``--cache`` if given, else ``$RASPUTIN_DATA/cache``, else
    ``$XDG_CACHE_HOME/rasputin``, else ``$HOME/.cache/rasputin``, else None.
    The XDG Base Directory rule: an empty or relative variable is ignored.
    Pure in its two arguments (no ``Path.home()``); ``mesh`` passes ``os.environ``."""
    if option is not None:
        return option
    if (root := data_root(environ)) is not None:
        return root / "cache"
    base = xdg_base(environ, "XDG_CACHE_HOME", ".cache")
    return None if base is None else base / "rasputin"


def _known_sources(environ: Mapping[str, str]) -> Mapping[str, RemoteSource]:
    """The catalogue's sources with the user's from their sources file
    (`user_source_file`) and the plugins'; a user file it refuses is a usage
    error naming it, a plugin file skipped a line on stderr."""
    try:
        found = catalogue(user_source_file(environ), plugin_source_files())
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    for line in found.refused:
        typer.echo(f"rasputin: {line}", err=True)
    return found.sources


def _dem_sources(
    dem: list[str], root: Path | None, known: Mapping[str, RemoteSource]
) -> tuple[tuple[Path, ...], CachedSource | None, RemoteSource | None]:
    """``--dem`` as given, with the catalogue entry when it names one: one key
    alone is that source, a local one as its path, a remote one as its cache;
    anything else is a path, so ``./glo30`` is a file."""
    keys = [d for d in dem if d in known]
    if not keys:
        return tuple(Path(d) for d in dem), None, None
    if len(dem) > 1:
        raise typer.BadParameter(
            f"{keys[0]} is a catalogue source; give it alone, not with other --dem",
            param_hint="--dem",
        )
    try:
        check_values(known[keys[0]], "heights", "--dem")  # a class source is not a DEM
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--dem") from exc
    entry = known[keys[0]]
    if entry.path is not None:
        return (entry.path,), None, entry
    if root is None:
        raise typer.BadParameter(
            NO_CACHE,
            param_hint="--dem",
        )
    return (), CachedSource(source=keys[0], cache=root), entry


def _dem_repository(
    dem: list[str], cache: Path | None, res: Resources, z_unit: ZSpelling | None
) -> DemRepository:
    """``--dem`` of the catchment commands as a repository, by path or by
    source name; the height unit is ``--z-unit``, else the source's."""
    known = _known_sources(os.environ)
    paths, cached, entry = _dem_sources(dem, cache_root(cache, os.environ), known)
    caller = _z_unit(z_unit) or (entry.z_unit if entry else None)
    try:
        return repository_for(paths, cached=cached, threads=res.threads, z_unit=caller)
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(
            cli_kit.words(exc) if isinstance(exc, ValueError) else str(exc), param_hint="--dem"
        ) from exc


def _open_dem(
    dem: tuple[Path, ...],
    bbox: CrsBox | None,
    clock: PhaseClock,
    domain: DomainPolygon | None = None,
    cached: CachedSource | None = None,
    target_crs: str | None = None,
    res: Resources | None = None,
    z_unit: ZUnit | None = None,
) -> DemPlan:
    """``--dem`` and ``--bbox`` or the read ``--domain`` planned, from headers;
    every refusal, the reader's, the mosaic's or the domain's extent, is a
    usage error in its own words."""
    with clock.phase("decode"), _dem_refusals(dem):
        request = DemRequest(
            sources=dem,
            cached=cached,
            box=bbox,
            domain=domain,
            target_crs=target_crs,
            z_unit=z_unit,
        )
        return plan_dem(request, res)


@contextmanager
def _dem_refusals(dem: tuple[Path, ...]) -> Iterator[None]:
    """A refusal of the DEM's plan or read made a usage error in its own words."""
    try:
        yield
    except NotCached:
        raise
    except OutCrsRequiredError as exc:  # one unwrapped line, outside the panel, to paste
        typer.echo(f"--out-crs {shlex.quote(exc.suggestion.proj)}", err=True)
        raise typer.BadParameter(
            f"{exc.head}; the suggested --out-crs is the line above", param_hint="--dem"
        ) from exc
    except OSError as exc:
        where = exc.filename or ", ".join(map(str, dem))
        raise typer.BadParameter(
            f"cannot read {where}: {exc.strerror or exc}", param_hint="--dem"
        ) from exc
    except DomainError as exc:
        raise typer.BadParameter(str(exc), param_hint="--domain") from exc
    except BoxError as exc:
        raise typer.BadParameter(str(exc), param_hint="--bbox-crs") from exc
    except NoAreaError as exc:
        raise typer.BadParameter(str(exc), param_hint="--bbox") from exc
    except TargetRegionError as exc:
        raise typer.BadParameter(str(exc), param_hint="--out-crs") from exc
    except ValueError as exc:
        raise typer.BadParameter(cli_kit.words(exc), param_hint="--dem") from exc


def _cut_partition(
    planned: DemPlan, tolerance: float, cells: Partition, pieces: int, budget: int
) -> Partitioned:
    """The partition of a window of more than one cell, on its lattice."""
    from tin_engine import basin_run

    window, grid = planned.mosaic.window, planned.grid
    ask = basin_run.CutAsk(pieces, budget, planned.mosaic.reference, window.row0, window.col0)
    if grid is not None:  # node (R, K) of a resampled grid is at (K h, -R h)
        ask = replace(ask, origin=(0.0, 0.0), row0=grid.row0, col0=grid.col0)
    try:
        return basin_run.cut_partition(planned.meta, tolerance, cells, ask)
    except RunError as exc:
        raise typer.BadParameter(str(exc), param_hint="--dem") from exc


def _z_unit(option: ZSpelling | None) -> ZUnit | None:
    """``--z-unit``'s spelling to the unit's name."""
    return next((row[0] for row in UNITS.values() if row[2] == option), None)


def _z_words(meta: RasterMeta, option: str | None, caller: ZUnit | None) -> str:
    """The record's ``dem_vertical_unit``: metres, converted from what, and who
    said so: the file or the catalogue entry (no words), --z-unit, or nobody."""
    unit = BY_NAME[meta.z_unit][3]
    words = "metres" if meta.z_unit == "metre" else f"metres, converted from {unit}"
    if not meta.vertical_unit_assumed or (caller is not None and option is None):
        return words
    return words + (" (given with --z-unit)" if option else " (assumed: the DEM file does not say)")


def _tolerance_lines(
    far: float | None,
    near: tuple[Path, float] | None,
    crs: str | None,
    ramp: tuple[float, float] | None,
) -> ToleranceLines | None:
    """The three line-tolerance flags checked, as one spec or None."""
    if near is None:
        if crs is not None:
            raise typer.BadParameter(
                "applies only with --tolerance-near", param_hint="--tolerance-near-crs"
            )
        if ramp is not None:
            raise typer.BadParameter(
                "--tolerance-ramp needs --tolerance-near", param_hint="--tolerance-ramp"
            )
        return None
    hint = "--tolerance-near"
    if far is None:
        raise typer.BadParameter("--tolerance-near needs --tolerance", param_hint=hint)
    if ramp is None:
        raise typer.BadParameter("--tolerance-near needs --tolerance-ramp", param_hint=hint)
    (path, n), (start, end) = near, ramp
    for flag, values in ((hint, (n,)), ("--tolerance-ramp", ramp)):
        if not all(math.isfinite(v) and v >= 0 for v in values):
            got = " ".join(f"{v:g}" for v in values)
            raise typer.BadParameter(f"must be finite and >= 0, got {got}", param_hint=flag)
    if n > far:
        raise typer.BadParameter(f"N {n:g} is above --tolerance {far:g}", param_hint=hint)
    if start > end:
        raise typer.BadParameter(
            f"START {start:g} is above END {end:g}", param_hint="--tolerance-ramp"
        )
    return ToleranceLines(path=path, crs=crs, near_m=n, start_m=start, end_m=end)


def _tolerance_slope(
    far: float | None, given: tuple[float, float, float] | None
) -> ToleranceSlope | None:
    """The slope-tolerance flag checked in its order, as one spec;
    None without it, and when N equals --tolerance (then stderr says so)."""
    if given is None:
        return None
    hint = "--tolerance-slope"
    if far is None:
        raise typer.BadParameter("--tolerance-slope needs --tolerance", param_hint=hint)
    n, start, end = given
    if not all(math.isfinite(v) for v in given):
        got = " ".join(f"{v:g}" for v in given)
        raise typer.BadParameter(f"must be finite, got {got}", param_hint=hint)
    for broken, words in (
        (n <= 0, f"N must be above 0, got {n:g}"),
        (n > far, f"N {n:g} is above --tolerance {far:g}"),
        (start < 0, f"START {start:g} is below 0"),
        (start > end, f"START {start:g} is above END {end:g}"),
        (end >= 90, f"END {end:g} must be below 90"),
    ):
        if broken:
            raise typer.BadParameter(words, param_hint=hint)
    if n == far:
        typer.echo(f"{hint}: N equals --tolerance, so slope changes nothing", err=True)
        return None
    return ToleranceSlope(near_m=n, start_deg=start, end_deg=end)


def _line_field(
    lines: ToleranceLines, tile: DemTile, far: float
) -> tuple[LineTolerance | None, dict[str, Value]]:
    """The line tolerance field on ``tile``'s grid, None when no line is in
    reach (then stderr says so), and the record's three input entries."""
    m = tile.meta
    window = (
        m.x_min,
        m.y_max - (m.rows - 1) * m.delta_y,
        m.x_min + (m.cols - 1) * m.delta_x,
        m.y_max,
    )
    try:
        segs, margin = line_segments(lines, window, m.crs)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--tolerance-near") from exc
    values: dict[str, Value] = {
        "tolerance_near_m": lines.near_m,
        "tolerance_ramp_m": f"{lines.start_m:g} to {lines.end_m:g}",
        "tolerance_lines": f"{lines.path.name}: {plural(len(segs), 'segment', 'segments')}",
    }
    if len(segs) == 0:
        typer.echo(
            f"no tolerance lines within {lines.end_m:g} m of the domain; "
            f"every triangle is held to {far:g} m",
            err=True,
        )
        return None, values
    near, start, end = lines.near_m, lines.start_m, lines.end_m
    try:
        return LineTolerance(to_core(tile), segs, near, far, start, end, margin), values
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--tolerance-near") from exc


@dataclass(frozen=True, slots=True)
class _DemMesh:
    """What ``_dem_mesh`` made: the mesh, then what ``--stats`` reports about
    the run (None where a row does not apply), and the record entries it knows
    (``run_record``'s builder arguments). A cut run's meshes are its pieces',
    in ``cut``, and ``trimmed`` is then empty."""

    trimmed: Trimmed
    meta: RasterMeta
    domain_vertices: int | None
    domain_holes: int | None
    start_vertices: int | None
    start_triangles: int | None
    values: Mapping[str, Value]
    feature_counts: tuple[int, int] | None = None
    cut: RunCut | None = None


#: A cut run's `_DemMesh.trimmed`: its meshes are its pieces'.
_NO_MESH = Trimmed(np.zeros((0, 3)), np.zeros((0, 3), np.uint32), np.zeros((0, 2), np.uint32),
                   np.zeros(0, np.uint32), 0)  # fmt: skip


@dataclass(frozen=True, slots=True)
class _Tuning:
    """What every refinement call of a run takes besides ``How``: the start's
    quality gain, the line and slope tolerances, and the thread count."""

    gain: float = -1.0
    lines: ToleranceLines | None = None
    steep: ToleranceSlope | None = None
    threads: int = 0


def _dem_mesh(
    held: Once[DemTile] | RasterMeta,
    dem: str,
    stride: int | None,
    delaunay: bool,
    spacing: float,
    tolerance: float | None,
    clock: PhaseClock,
    domain: DomainPolygon | None = None,
    domain_name: str = "",
    min_angle: float = 0.0,
    feet: bool = False,
    features: FeatureSet | None = None,
    grid: TargetGrid | None = None,
    checks: Iterator[Block] | None = None,
    gain: float = -1.0,
    lines: ToleranceLines | None = None,
    steep: ToleranceSlope | None = None,
    threads: int = 0,
    vocabulary: EdgeVocabulary = DEFAULT_VOCABULARY,
    cut: Partitioned | None = None,
    budget: Window | None = None,
    piece_budget: int | None = None,
    rms: bool = False,
) -> _DemMesh:
    """Subsample, triangulate, sample or refine, and trim ``held``'s tile
    (with ``cut``, ``held`` is the window's grid alone, and the cut run is
    returned to be run).

    Without ``tolerance``, z is sampled bilinearly at the stride grid. With it,
    the stride grid is the start mesh, refined against the DEM's nodes; with
    ``domain`` (already in the DEM's CRS; ``domain_name`` is its file's) the
    start is the polygon's rings and ``features``' lines instead.
    ``min_angle`` > 0 improves the start's angles first, a point accepted only
    for a gain of ``gain`` degrees when ``gain`` >= 0; ``feet`` inserts a
    node's foot on a nearby constraint segment. With ``grid`` and its
    ``checks``, the refined mesh is checked against the source's nodes after
    the tile is dropped. ``lines`` and ``steep`` build the line and slope
    tolerances on the tile's grid, passed to every refinement call. With
    ``cut``, a partition of more than one cell, the start has the partition
    lines too, and the run is cut into pieces. With ``budget`` the tolerance
    is searched for (``_budget_refined``, which checks ``piece_budget``).
    With ``rms``, or a budget, the error summary of the mesh is measured.
    Returns the mesh, the record entries it knows, and the ``--stats`` sizes;
    ``clock`` gets the phases.
    ``dem`` names the source in messages. Every refusal is a usage error in the
    engine's own words, and no file is written.
    """
    tile = held.take() if isinstance(held, Once) else None
    tolerance = FLOOR_TOLERANCE if budget is not None else tolerance  # each trial sets its own
    meta = held if isinstance(held, RasterMeta) else cast(DemTile, tile).meta
    values: dict[str, Value] = {}
    domain_vertices = domain_holes = None
    feature_counts = None
    if domain is not None:
        started = start_chains(domain, features.features if features else (), vocabulary)
        chains = [(indices, ROLES[role], mask) for indices, role, mask in started.chains]
        rings = (domain.polygon.exterior, *domain.polygon.interiors)
        domain_vertices, domain_holes = sum(len(r.coords) - 1 for r in rings), len(rings) - 1
        shape = f"{plural(domain_holes, 'hole', 'holes')}, {domain_vertices} vertices"
        values["domain"] = f"{domain_name}: 1 outline, {shape}"
        values["start_mesh"] = "the domain outline"
        vertices = started.vertices
        if cut is not None:  # the partition lines, as breaklines with the seam bit
            seam = PIECE_VOCABULARY.mask("seam")
            for xy in cut_lines(domain.polygon, cut.xs, cut.ys):
                at = list(range(len(vertices), len(vertices) + len(xy)))
                chains.append((at, ChainRole.Breakline, seam))
                vertices = np.concatenate([vertices, xy])
        run = _engine(vertices, chains, delaunay, spacing, clock)
        if features is not None:
            values["start_mesh"] = "the domain outline and the feature lines"
            feature_counts = (
                len(started.chains) - len(rings),
                len(started.vertices) - domain_vertices,
            )
            noded = len(run.noded.vertices) if run.noded is not None else 0
            typer.echo(
                f"lines: {len(started.vertices)} vertices read, {noded} after joining "
                "shared edges and adding crossings",
                err=True,
            )
    else:
        if stride is not None:
            step = stride
        else:
            step = default_stride(meta) if tolerance is None else refine_start_stride(meta)
        xy, ring = subsample(meta, step)
        values["start_mesh"] = ordinal(step)
        run = _engine(xy, [(ring, ChainRole.Outer, 0)], delaunay, spacing, clock)
    if run.mesh is None or run.noded is None:
        raise typer.BadParameter(f"{dem} has no mesh to write: {run.status}. {run.message}")

    with clock.phase("start mesh: constraint edges"):
        edges, masks = _constraint_arrays(run.mesh, run.noded)
    refined = tolerance is not None
    if tolerance is None:
        mesh_xy = np.asarray(run.mesh.vertices)
        assert tile is not None  # only a cut run comes without its tile
        with clock.phase("sample"):
            z, valid = sample(to_core(tile), mesh_xy)
        with clock.phase("trim"):
            trimmed = trim(
                vertices=mesh_xy,
                triangles=np.asarray(run.mesh.triangles),
                edges=edges,
                edge_masks=masks,
                z=z,
                valid=valid,
            )
    else:
        how = How(dem, tolerance, min_angle, feet, grid, checks, budget, rms and budget is None)
        refine_with = functools.partial(_refined, tuning=_Tuning(gain, lines, steep, threads))
        if cut is not None:
            from tin_engine.basin_run import cut_run

            start, part = run.mesh, cut

            def run_cut(
                load: Callable[[IndexWindow], PieceDem], write: Callable[[Piece], PieceEntry]
            ) -> CutRun:
                return cut_run(
                    load, start, edges, masks, part, spacing, how, clock, refine_with, write
                )

            return _DemMesh(
                _NO_MESH, meta, domain_vertices, domain_holes, None, None, values,
                feature_counts, run_cut,
            )  # fmt: skip
        assert tile is not None
        if budget is not None:
            made = _budget_refined(
                tile, run.mesh, edges, masks, how, clock, refine_with, piece_budget, threads
            )
        else:
            once = Once(tile)
            del tile  # `_refined` drops the last reference before the final check
            made = refine_with(once, run.mesh, edges, masks, how, clock)
        trimmed = made.trimmed
        values |= made.values
    if len(trimmed.triangles) == 0:
        raise typer.BadParameter(
            f"{dem} has no data under any triangle; nothing to write", param_hint="--dem"
        )
    values["nodata_vertices_removed"] = trimmed.dropped
    return _DemMesh(
        trimmed=trimmed,
        meta=meta,
        domain_vertices=domain_vertices,
        domain_holes=domain_holes,
        start_vertices=len(run.mesh.vertices) if refined else None,
        start_triangles=len(run.mesh.triangles) if refined else None,
        values=values,
        feature_counts=feature_counts,
    )


def _budget_refined(
    tile: DemTile,
    start: IndexedMesh2,
    edges: npt.NDArray[np.uint32],
    masks: npt.NDArray[np.uint32],
    how: How,
    clock: PhaseClock,
    refine_with: Callable[..., Refined],
    piece_budget: int | None,
    threads: int = 0,
) -> Refined:
    """The mesh of the tolerance whose written count is in ``how.budget``:
    each trial is one ``refine_with`` run at its tolerance, the check points
    filed once for all. The trial with the most triangles not over the cap
    is kept, so the chosen one is never run twice. With ``piece_budget``, a
    trial whose tolerance would cut the window into pieces refuses the run
    before it is refined. The kept mesh's error summary is measured once."""
    window = how.budget
    assert window is not None
    if how.grid is not None and how.checks is not None and not isinstance(how.checks, CheckPoints):
        how = replace(how, checks=final_check.stored(how.grid, how.checks, clock))
    kept: Refined | None = None

    def trial(tolerance: float) -> Trial:
        nonlocal kept
        if piece_budget is not None:
            try:
                cells = partition(
                    tile.meta.cols, tile.meta.rows, tolerance, piece_budget=piece_budget
                )
            except ValueError as exc:
                raise RunError(str(exc)) from exc
            if cells.nx * cells.ny > 1:
                t = f"{tolerance:.3g}"
                raise typer.BadParameter(
                    f"--triangles {window.asked} needs a tolerance of {t} m here, which would cut "
                    f"the window into {cells.nx} \N{MULTIPLICATION SIGN} {cells.ny} pieces under "
                    f"the piece budget of {piece_budget} bytes; the triangle budget does not run "
                    f"on pieces yet: use --tolerance {t} with --domain, or raise --piece-budget",
                    param_hint="--triangles",
                )
        made = refine_with(
            Once(tile), start, edges, masks, replace(how, tolerance=tolerance), clock
        )
        count = len(made.trimmed.triangles)
        if count <= window.cap and (kept is None or count > len(kept.trimmed.triangles)):
            kept = made
        at = made.values.get("dem_nodes_at_vertices_max_error_m") or 0.0
        return Trial(
            tolerance, count, max(float(cast(float, made.values["max_error_m"])), float(at))
        )

    try:
        found = search(trial, window)
    except UnderFloor as exc:
        f, e = exc.floor.triangles, exc.floor.max_error
        raise typer.BadParameter(
            f"the fewest triangles this DEM gives is {f} triangles, at a largest error of "
            f"{e:.3g} m; ask for at least {f}",
            param_hint="--triangles",
        ) from exc
    assert kept is not None
    tried = "; ".join(f"{t.tolerance:.3g} m: {t.triangles}" for t in found.trials)
    store = how.checks if how.grid is not None and isinstance(how.checks, CheckPoints) else None
    with clock.phase("error summary"):
        measured = error_values(kept.final, to_core(tile) if store is None else store, threads)
    values = (
        dict(kept.values)
        | measured
        | {
            "tolerance_m": found.chosen.tolerance,
            "triangles_asked": window.asked,
            "triangles_cap": window.cap,
            "budget_trials": tried,
            "budget_window_missed": "yes" if found.missed else None,
        }
    )
    return replace(kept, values=values)


def _refined(
    held: Once[DemTile],
    start: IndexedMesh2,
    edges: npt.NDArray[np.uint32],
    masks: npt.NDArray[np.uint32],
    how: How,
    clock: PhaseClock,
    frozen: int = 0,
    heights: Mapping[int, float] | None = None,
    *,
    tuning: _Tuning,
) -> Refined:
    """Refine from ``start``, then the edge strip, or the final check with a
    grid, which runs after the tile is dropped; then trim. No vertex goes on
    an edge whose mask meets ``frozen``, and no strip point lies on one;
    ``heights`` are a piece's seam records' z, written over the refined z of
    those vertices."""
    tile = held.take()
    meta = tile.meta
    tolerance, feet, threads = how.tolerance, how.feet, tuning.threads
    values: dict[str, Value] = {}
    field = None
    if tuning.lines is not None:
        field, line_values = _line_field(tuning.lines, tile, tolerance)
        values |= line_values
    slope, steep = None, tuning.steep
    if steep is not None:
        with clock.phase("slope"):
            ramp = (steep.near_m, tolerance, steep.start_deg, steep.end_deg)
            slope = SlopeTolerance(to_core(tile), *ramp, threads=threads)
    t0 = time.perf_counter_ns()
    out = refine(
        to_core(tile),
        start,
        edges,
        masks,
        tolerance=tolerance,
        min_angle_deg=how.min_angle,
        constraint_feet=feet,
        frozen_mask=frozen,
        min_gain_deg=tuning.gain,
        field=field,
        slope=slope,
        threads=threads,
    )
    _refine_phases(clock, (time.perf_counter_ns() - t0) / 1e9, out)
    if not out.ok():
        raise typer.BadParameter(f"{how.dem}: {out.message}", param_hint="--dem")
    strip = edge_strip.generate(to_core(tile), out, clock, frozen)  # while the tile is held
    if how.grid is None or how.checks is None:
        final = edge_strip.run(
            to_core(tile), strip, out, tolerance, clock, feet, field, slope, threads, frozen
        )
        if how.rms and final.ok():
            with clock.phase("error summary"):
                values |= error_values(final, to_core(tile), threads)
        del tile
        # Refine's maximum and the strip run's make an upper bound.
        max_error = max(out.max_error, final.max_error)
        near_error = max(out.max_error_near, final.max_error_near)
        share = max(out.max_slope_share, final.max_slope_share)
        inside, tight = out.slope_nodes
        values["line_check_dem_nodes_inserted"] = final.nodes_inserted
    else:
        del tile  # the final check runs without the target tile
        store = how.checks
        if not isinstance(store, CheckPoints):
            store = final_check.stored(how.grid, store, clock)
        final, n = final_check.run(
            out, store, tolerance, clock, strip, feet, field, slope, threads, frozen
        )
        if how.rms and final.ok():
            with clock.phase("error summary"):
                values |= error_values(final, store, threads)
        max_error, near_error = final.max_error, final.max_error_near
        share = final.max_slope_share
        inside, tight = final.slope_nodes
        values |= {
            "resampled_grid_max_error_m": out.max_error,
            "dem_nodes_checked": n,
            "dem_check_points_inserted": final.inserted,
            "dem_check_rounds": final.rounds,
        }
    if not final.ok():
        raise typer.BadParameter(f"{how.dem}: {final.message}", param_hint="--dem")
    values |= {
        "dem_nodes_at_vertices": final.coincident,
        "dem_nodes_at_vertices_max_error_m": final.coincident_max_error,
        "line_points_checked": strip.size,
        "line_max_error_m": final.strip_max_error,
        "line_points_on_nodata": strip.no_data,
        "line_points_refused": final.strip_refused,
        "line_points_refused_max_error_m": final.strip_refused_max_error,
        "line_points_inserted": final.strip_inserted,
        "line_points_duplicate": strip.duplicates,
        "triangles_flat_as_written": final.flat_as_written,
        "edges_not_delaunay_as_written": final.illegal_as_written,
    }
    z = np.asarray(final.z)
    if heights:
        z = z.copy()
        z[list(heights)] = list(heights.values())
    with clock.phase("trim"):
        trimmed = trim(
            vertices=final.vertices,
            triangles=final.triangles,
            edges=final.edges,
            edge_masks=final.masks,
            z=z,
            valid=final.valid,
        )
    values |= {
        "start_min_angle_deg": how.min_angle,
        "start_quality_gain_deg": tuning.gain,
        "snap_to_lines": "on" if feet else "off",
        "tolerance_m": tolerance,
        "max_error_m": max_error,
        "max_error_near_lines_m": near_error if field is not None else None,
        "dem_nodes_outside_mesh": final.uncovered,
        "refinement_rounds": out.rounds,
        "points_inserted": out.inserted,
        "points_inserted_on_nodata": out.carved,
        "edge_flips": out.flips,
        "start_quality_points_inserted": out.quality_inserted,
        "start_quality_points_skipped": out.quality_skipped,
        "start_quality_points_without_gain": out.quality_no_gain,
        "start_quality_points_snapped_to_lines": out.quality_feet,
        "start_quality_lines_split": out.quality_line_splits,
        "points_snapped_to_lines": out.feet,
        "snaps_refused": out.feet_refused,
        "final_check_points_snapped_to_lines": final.feet,
        "final_check_snapped_points_added_anyway": final.feet_fallback,
        "start_vertices_between_dem_nodes": _off_node(np.asarray(start.vertices), meta),
        "nodata_vertices_removed": trimmed.dropped,
    }
    if steep is not None:
        values |= {
            "tolerance_slope_m": steep.near_m,
            "tolerance_slope_deg": f"{steep.start_deg:g} to {steep.end_deg:g}",
            "slope_nodes_tightened": (
                f"{tight} of {inside} DEM nodes ({100 * tight / max(inside, 1):.1f} %)"
            ),
            "max_error_slope_share_of_tolerance": share,
        }
    return Refined(out, final, z, trimmed, values, final.on_frozen)


def _refine_phases(clock: PhaseClock, seconds: float, out: RefineOutcome) -> None:
    """``refine`` and its sub-rows; setup + output is the remainder."""
    inner = (
        ("refine: legalise start", out.legalise_seconds),
        ("refine: start quality", out.quality_seconds),
        ("refine: scan (parallel)", out.scan_seconds),
        ("refine: split + flip (serial)", out.split_seconds),
    )
    clock.add("refine", seconds)
    for name, part in inner:
        clock.add(name, part)
    clock.add("refine: setup + output", max(0.0, seconds - sum(p for _, p in inner)))


def _off_node(xy: npt.NDArray[np.float64], meta: RasterMeta) -> int:
    """How many of ``xy`` are not a DEM node bit for bit, as ``refine`` classifies."""
    x, y = meta.node_xy(*np.round(meta.index_of(xy[:, 0], xy[:, 1])))
    node = (x == xy[:, 0]) & (y == xy[:, 1])
    return int(np.count_nonzero(~node))


def _placed(
    rivers: Path,
    repository: DemRepository,
    seed: tuple[float, float],
    seed_crs: str,
    request: BatchRequest,
) -> tuple[tuple[Placement, str | None], CatchmentRequest]:
    """The gauge placed on the river file's nearest line (no watercourse
    number), with its line's river name, and the catchment request seeded by
    its reach in the file's CRS, which `delineate` moves into the DEM's; `U`'s
    floor is three of the DEM's cells at the station."""
    segments, crs, _ = _segments(rivers)
    ((x, y),) = reprojector(seed_crs, crs)([seed])
    cell = cell_at(repository, seed, seed_crs)
    try:
        placement, _, made = seed_for(Gauge(x=x, y=y), segments, crs, request, cell_m=cell)
    except NoRiverLine as exc:
        raise typer.BadParameter(str(exc), param_hint="--rivers") from exc
    assert placement is not None  # without lakes, `seed_for` places or refuses
    name = next(s.name for s in segments if s.objectid == placement.objectid)
    return (placement, name), made


def _reach_crs(crs: str, repository: DemRepository) -> None:
    """`check_reach_crs`, its refusal naming --rivers."""
    try:
        check_reach_crs(crs, repository)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--rivers") from exc


@contextmanager
def _writing(path: Path) -> Iterator[None]:
    """A write to `path` whose `OSError` becomes a refusal naming --out-dir."""
    try:
        yield
    except OSError as exc:
        raise typer.BadParameter(
            f"cannot write {path}: {exc.strerror or exc}", param_hint="--out-dir"
        ) from exc


def _segments(rivers: Path) -> tuple[tuple[RiverSegment, ...], str, int]:
    """The river file's segments, its CRS and the exact copies dropped, said
    on stderr; a file that cannot be read is refused naming --rivers."""
    try:
        segments, crs, dropped = read_segments(rivers)
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(str(exc), param_hint="--rivers") from exc
    copies = plural(dropped, "exact copy", "exact copies")
    typer.echo(
        f"rivers: {plural(len(segments) + dropped, 'segment', 'segments')} read, "
        f"{copies} dropped (same river, same vertices to 1 cm)",
        err=True,
    )
    return segments, crs, dropped


def _placement_report(
    placed: tuple[Placement, str | None], gauge: GaugeResult
) -> dict[str, object]:
    """Say where the gauge went, and return the placement and sensitivity
    properties for the catchment file."""
    (p, name), s, causes = placed, gauge.sensitivity, gauge.causes
    verdict = f"uncertain ({', '.join(causes)})" if causes else "well defined"
    typer.echo(
        f"placed on the river line {p.distance_m:.0f} m from the station "
        f"({'' if name is None else f'river {name}, '}line {p.objectid}), moved "
        f"{gauge.node_offset_m:.0f} m onto the DEM's valley floor; {s.a0:.4g} km2 drain "
        f"through it, and the area changes by {100 * s.swing:.1f} % within "
        f"{p.reach.uncertainty:.0f} m up and down the river: {verdict}",
        err=True,
    )
    return {**p.report_fields(), **gauge.report_fields()}


#: What the suffix of ``catchment --out`` may be: GeoJSON, which ``--domain`` reads.
CATCHMENT_SUFFIXES = (".geojson", ".json")


def catchment(
    out: Annotated[
        Path,
        typer.Option("--out", writable=True, help="Where to write the polygon: .geojson or .json."),
    ],
    dem: Annotated[list[str], typer.Option("--dem", help=f"{DEM_HELP}.")],
    seed: Annotated[
        tuple[float, float],
        typer.Option("--seed", metavar="X Y", help="A point in the lake, in --seed-crs."),
    ],
    seed_crs: SEED_CRS = None,
    lakes: Annotated[
        Path | None,
        typer.Option(
            "--lakes",
            help="Any polygon file (.gpkg, .geojson or .json, in any CRS; lines and points "
            "in it are skipped): the lake polygon under the seed is the seed. Without it, "
            "the seed is the DEM node nearest the point.",
        ),
    ] = None,
    lakes_layer: Annotated[
        str | None, typer.Option("--lakes-layer", help="The GeoPackage's features table.")
    ] = None,
    outline_tolerance: Annotated[
        float | None,
        typer.Option(
            "--outline-tolerance",
            metavar="METRES",
            help="Reduce the outline, keeping its area, to within this of the fine one; "
            "0 writes the fine outline without its collinear vertices. Default: twice the "
            "DEM's cell.",
        ),
    ] = None,
    rivers: Annotated[
        Path | None,
        typer.Option(
            "--rivers",
            help="River lines (GeoJSON, in any CRS): the seed is a gauge, placed on the "
            "nearest line and on the DEM's flow path along it.",
        ),
    ] = None,
    map_radius: Annotated[
        float | None,
        typer.Option(
            "--map-radius",
            metavar="METRES",
            help="How far a river line may be from the gauge. Default 500.",
        ),
    ] = None,
    reach_up: Annotated[
        float | None,
        typer.Option(
            "--reach-up",
            metavar="METRES",
            help="How much river above the gauge is burnt into the DEM. Default 1000.",
        ),
    ] = None,
    cache: CACHE = None,
    threads: THREADS = None,
    memory_budget: MEMORY_BUDGET = None,
    z_unit: Z_UNIT = None,
) -> None:
    """Write the catchment of a lake or a gauge, from the DEM, as a GeoJSON polygon.

    The polygon is in the DEM's CRS, and mesh --domain reads it. The fine
    outline is drawn between DEM nodes, holes filled, then reduced to
    --outline-tolerance keeping its area, with the seed inside. A catchment cut
    by the data's edge or by NoData is refused, and nothing is written. With
    --rivers the seed is a gauge on a mapped river.
    """
    if seed_crs is None:
        raise typer.BadParameter(NO_SEED_CRS, param_hint="--seed-crs")
    if lakes is None and lakes_layer is not None:
        raise typer.BadParameter("applies only with --lakes", param_hint="--lakes-layer")
    for name, value in (("--map-radius", map_radius), ("--reach-up", reach_up)):
        if value is not None and rivers is None:
            raise typer.BadParameter("applies only with --rivers", param_hint=name)
        if value is not None and not (math.isfinite(value) and value > 0.0):
            raise typer.BadParameter("must be a finite number above 0", param_hint=name)
    if rivers is not None and lakes is not None:
        raise typer.BadParameter(
            "cannot be combined with --lakes: a lake is already the seed", param_hint="--rivers"
        )
    if outline_tolerance is not None and not (
        math.isfinite(outline_tolerance) and outline_tolerance >= 0.0
    ):
        raise typer.BadParameter("must be finite and at least 0", param_hint="--outline-tolerance")
    if out.suffix.lower() not in CATCHMENT_SUFFIXES:
        raise typer.BadParameter(f"use {' or '.join(CATCHMENT_SUFFIXES)}", param_hint="--out")
    target = cli_kit.destination(out)
    try:
        found = None if lakes is None else read_lake_polygons(lakes, lakes_layer, seed, seed_crs)
    except FeatureError as exc:
        raise typer.BadParameter(str(exc), param_hint="--lakes") from exc
    res = cli_kit.cli_resources(threads, memory_budget)
    repository = _dem_repository(dem, cache, res, z_unit)
    try:
        placement: tuple[Placement, str | None] | None = None
        if rivers is not None:
            batch = BatchRequest(
                map_radius=map_radius or 500.0,
                reach_up=reach_up or 1000.0,
                outline_tolerance=outline_tolerance,
            )
            placement, request = _placed(rivers, repository, seed, seed_crs, batch)
        else:
            request = CatchmentRequest(
                seed=seed,
                seed_crs=seed_crs,
                lakes=None if found is None else found[0],
                lakes_crs=None if found is None else found[1],
                outline_tolerance=outline_tolerance,
            )
        result = delineate(request, repository, res)
    except NotCached as exc:
        raise typer.BadParameter(
            f"{exc}{_bbox_fetch(dem[0], exc, cache)}", param_hint="--dem"
        ) from exc
    except OSError as exc:
        raise typer.BadParameter(f"cannot read {exc.filename}: {exc}", param_hint="--dem") from exc
    except ValueError as exc:
        hint = "--lakes" if isinstance(exc, LakeError) else "--dem"
        raise typer.BadParameter(cli_kit.words(exc), param_hint=hint) from exc
    degrees = result.meta.geographic
    for k, w in enumerate(result.windows, 1):
        b = w.bounds
        grown = f"window widened to the {', '.join(w.grown)}"
        where = (
            f"longitude {b.x_min:.5f} to {b.x_max:.5f}, latitude {b.y_min:.5f} to {b.y_max:.5f}"
            if degrees
            else f"x {b.x_min:.0f}-{b.x_max:.0f}, y {b.y_min:.0f}-{b.y_max:.0f}"
        )
        typer.echo(
            f"window {k}: {where}, {w.rows} x {w.cols} nodes, searched in {w.seconds:.2f} s, "
            + (grown if w.grown else "catchment inside the window"),
            err=True,
        )
    extra: dict[str, object] = {}
    if placement is not None and result.gauge is not None:
        extra = _placement_report(placement, result.gauge)
    elif result.lake_area is None:
        typer.echo(
            f"start: the outlet node at {tuple(map(float, result.seed))} (an outlet must lie "
            "on the flow line; it is not moved there)",
            err=True,
        )
    else:
        typer.echo(
            f"start: lake of {result.lake_area / 1e6:.6f} km2, {result.seed_nodes} DEM nodes",
            err=True,
        )
    typer.echo(
        f"catchment: {result.nodes} DEM nodes, {result.nodes_area / 1e6:.6f} km2",
        err=True,
    )
    vertices = len(result.fine.exterior.coords) - 1
    typer.echo(
        f"outline along DEM cells: {vertices} vertices, {result.fine_area / 1e6:.6f} km2, "
        f"{plural(result.rings_dropped, 'separate patch', 'separate patches')} left out "
        f"({result.dropped_nodes} nodes), "
        f"{plural(result.holes_filled, 'enclosed gap', 'enclosed gaps')} filled "
        f"({result.holes_area / 1e6:.6f} km2), "
        f"traced in {result.trace_seconds:.2f} s",
        err=True,
    )
    reduced = result.reduced
    kept = len(reduced.exterior.coords) - 1
    change = result.reduced_area - result.fine_area
    typer.echo(
        f"reduced outline: {kept} vertices, {result.reduced_area / 1e6:.6f} km2, difference "
        f"{change:.3g} m2 ({change / result.fine_area:.2g} relative), tolerance "
        f"{result.tolerance:g} m, {result.reduce_seconds:.2f} s",
        err=True,
    )
    properties = {"seed": list(seed), "seed_crs": seed_crs, **_outline_properties(result), **extra}
    target.write_bytes(catchment_geojson(reduced, result.crs, properties))
    typer.echo(f"{target}")


def _outline_properties(result: Catchment) -> dict[str, object]:
    """The properties of a catchment file: the counts, the areas,
    the tolerance and the windows."""
    return {
        "nodes": result.nodes,
        "fine_vertices": len(result.fine.exterior.coords) - 1,
        "fine_area_m2": result.fine_area,
        "reduced_vertices": len(result.reduced.exterior.coords) - 1,
        "reduced_area_m2": result.reduced_area,
        "outline_tolerance_m": result.tolerance,
        "windows": [[w.rows, w.cols] for w in result.windows],
    }


@contextmanager
def _refused(hint: str, prefix: str = "") -> Iterator[None]:
    """A `ValueError` made a usage error naming `hint`, its words after `prefix`."""
    try:
        yield
    except ValueError as exc:
        raise typer.BadParameter(f"{prefix}{exc}", param_hint=hint) from exc


def network_catchment(
    out: Annotated[
        Path,
        typer.Option("--out", writable=True, help="Where to write the polygon: .geojson or .json."),
    ],
    network: Annotated[
        str,
        typer.Option(
            "--network",
            help="The units' polygons: a catalogue source of features or a .gpkg, .geojson "
            "or .json file, read in place.",
        ),
    ],
    id_attribute: Annotated[str, typer.Option("--id-attribute", help="The unit's own ID.")],
    down_attribute: Annotated[
        str, typer.Option("--down-attribute", help="The ID of the unit it drains to.")
    ],
    network_layer: Annotated[
        str | None, typer.Option("--network-layer", help="The GeoPackage's features table.")
    ] = None,
    upstream_area_attribute: Annotated[
        str | None,
        typer.Option(
            "--upstream-area-attribute",
            help="The area upstream of the unit's outlet, its own included; the joined area "
            "must be within 1 % of the seed unit's.",
        ),
    ] = None,
    area_unit: Annotated[
        Literal["m2", "km2"] | None,
        typer.Option("--area-unit", help="The unit of --upstream-area-attribute."),
    ] = None,
    name_attribute: Annotated[
        str | None, typer.Option("--name-attribute", help="The units' names, for the report.")
    ] = None,
    unit: Annotated[
        str | None,
        typer.Option("--unit", help="Seed: this unit and every unit upstream of it."),
    ] = None,
    seed: Annotated[
        tuple[float, float] | None,
        typer.Option("--seed", metavar="X Y", help="Seed: a point at a unit's outlet."),
    ] = None,
    seed_crs: SEED_CRS = None,
    lakes: Annotated[
        str | None,
        typer.Option(
            "--lakes",
            help="Any polygon source or file: the lake under the seed point is the seed.",
        ),
    ] = None,
    lakes_layer: Annotated[
        str | None, typer.Option("--lakes-layer", help="The GeoPackage's features table.")
    ] = None,
    outlets: Annotated[
        str | None,
        typer.Option(
            "--outlets", help="The units' outlet points, one per unit, keyed by --id-attribute."
        ),
    ] = None,
    outlets_layer: Annotated[
        str | None, typer.Option("--outlets-layer", help="The GeoPackage's features table.")
    ] = None,
    outlet_radius: Annotated[
        float | None,
        typer.Option(
            "--outlet-radius",
            metavar="METRES",
            help="How near the seed must be to its unit's outlet point. Default 200.",
        ),
    ] = None,
    check_dem: Annotated[
        str | None,
        typer.Option(
            "--check-dem",
            help=f"Also run the terrain route (catchment) from the seed on this DEM ({DEM_HELP}) "
            "and report where the two differ.",
        ),
    ] = None,
    cache: CACHE = None,
    threads: THREADS = None,
    memory_budget: MEMORY_BUDGET = None,
    z_unit: Z_UNIT = None,
) -> None:
    """Write the catchment joined from a published network of sub-catchment units.

    Each unit names its downstream unit; the catchment is the seed's unit and
    every unit upstream of it, joined as published, in the network's CRS, and
    mesh --domain reads it. A point seed must lie, or its lake must, within
    --outlet-radius of its unit's outlet point: a point inside a unit is
    refused. --check-dem compares the result with the terrain route.
    """
    if (unit is None) == (seed is None):
        raise typer.BadParameter("give exactly one of --unit and --seed", param_hint="--unit")
    if seed is not None and seed_crs is None:
        raise typer.BadParameter(NO_SEED_CRS, param_hint="--seed-crs")
    for name, value, other, given in (
        ("--seed-crs", seed_crs, "--seed", seed),
        ("--lakes", lakes, "--seed", seed),
        ("--outlets", outlets, "--seed", seed),
        ("--outlet-radius", outlet_radius, "--seed", seed),
        ("--check-dem", check_dem, "--seed", seed),
        ("--area-unit", area_unit, "--upstream-area-attribute", upstream_area_attribute),
        ("--upstream-area-attribute", upstream_area_attribute, "--area-unit", area_unit),
        ("--lakes-layer", lakes_layer, "--lakes", lakes),
        ("--outlets-layer", outlets_layer, "--outlets", outlets),
    ):
        if value is not None and given is None:
            raise typer.BadParameter(f"needs {other}", param_hint=name)
    radius = 200.0 if outlet_radius is None else outlet_radius
    if not (math.isfinite(radius) and radius >= 0.0):
        raise typer.BadParameter("must be finite and at least 0", param_hint="--outlet-radius")
    if out.suffix.lower() not in CATCHMENT_SUFFIXES:
        raise typer.BadParameter(f"use {' or '.join(CATCHMENT_SUFFIXES)}", param_hint="--out")
    target = cli_kit.destination(out)
    known, root = _known_sources(os.environ), cache_root(cache, os.environ)

    def read_in_place(value: str, layer: str | None, flag: str) -> tuple[Path, str | None]:
        """`flag`'s file and its layer: the one given, else the entry's."""
        run = _fetch_run(value, None, None, None, cache)
        path, entry_layer, _ = _features_source(value, flag, known, root, run, VECTOR_SUFFIXES)
        return path, layer or entry_layer

    network_file, network_layer = read_in_place(network, network_layer, "--network")
    lakes_file = outlets_file = None
    if lakes is not None:
        lakes_file, lakes_layer = read_in_place(lakes, lakes_layer, "--lakes")
    if outlets is not None:
        outlets_file, outlets_layer = read_in_place(outlets, outlets_layer, "--outlets")
    attributes = (upstream_area_attribute, area_unit, name_attribute)
    with _refused("--network"):
        crs, links = read_links(
            network_file, network_layer, id_attribute, down_attribute, *attributes
        )
    with _refused("--network", f"{network_file.name}, {id_attribute}: "):
        net = Network.build(links)
    ends = ", ".join(f"{n} to {v or 'null'}" for v, n in net.ends.items())
    typer.echo(f"network: {len(links)} units in {crs}; draining out of the file: {ends}", err=True)
    seed_unit, place = unit or "", dict[str, object]({"unit": unit})
    if seed_unit and seed_unit not in net.links:
        raise typer.BadParameter(f"{network_file.name} has no unit {unit}", param_hint="--unit")
    request: CatchmentRequest | None = None
    if seed is not None and seed_crs is not None:  # both, as refused above
        ((x, y),) = reprojector(seed_crs, crs)([seed])
        with _refused("--seed"):
            covering = units_at(network_file, network_layer, id_attribute, Point(x, y))
            seed_unit = home_unit(Point(x, y), covering)
        with _refused("--lakes"):
            found = (
                None
                if lakes_file is None
                else read_lake_polygons(lakes_file, lakes_layer, seed, seed_crs)
            )
            request = CatchmentRequest(
                seed=seed,
                seed_crs=seed_crs,
                lakes=None if found is None else found[0],
                lakes_crs=None if found is None else found[1],
            )
            geometry = Point(x, y) if found is None else _lake(request, crs)
        link, own_km2 = net.links[seed_unit], area_m2(dict(covering)[seed_unit], crs) / 1e6
        named = f"{seed_unit}{'' if link.name is None else f' ({link.name})'}"
        choices = (
            f"give --unit {seed_unit} for the whole unit, up to {own_km2:.2f} km2 of it below the "
            "seed, or use the terrain route (rasputin catchment), which cuts at the point"
        )
        if outlets_file is None:
            raise typer.BadParameter(
                f"the seed is in the unit {named}; give --outlets to check that it is at the "
                f"unit's outlet, or {choices}",
                param_hint="--seed",
            )
        with _refused("--outlets"):
            outlet = outlet_point(outlets_file, outlets_layer, id_attribute, seed_unit, crs)
        if outlet is None:
            raise typer.BadParameter(
                f"{outlets_file.name} has no outlet point for the seed's unit {named}; {choices}",
                param_hint="--seed",
            )
        distance = outlet_distance_m(geometry, outlet, crs)
        if not distance <= radius:
            raise typer.BadParameter(
                f"the seed is {distance:.0f} m from the outlet of its unit {named}, "
                f"{own_km2:.2f} km2, beyond --outlet-radius {radius:g} m; a point inside a unit "
                f"is not joined: {choices}",
                param_hint="--seed",
            )
        typer.echo(
            f"seed unit: {named}, {own_km2:.6f} km2, the seed {distance:.0f} m from its outlet",
            err=True,
        )
        place = {"seed": seed, "seed_crs": seed_crs, "home_unit": seed_unit, "distance_m": distance}
    with _refused("--network"):
        ids = net.upstream(seed_unit)
        geometries = unit_geometries(network_file, network_layer, id_attribute, ids)
        joined = join(list(geometries.values()), crs)
    typer.echo(
        f"joined: {plural(len(ids), 'unit', 'units')}, {joined.area_m2 / 1e6:.6f} km2, "
        f"{plural(joined.holes, 'hole', 'holes')} ({joined.holes_m2 / 1e6:.6f} km2, filled)",
        err=True,
    )
    stated = net.links[seed_unit].upstream_area_m2
    if stated is not None:
        with _refused("--network"):
            check_area(joined.area_m2, stated)
        typer.echo(f"area check: stated upstream area {stated / 1e6:.6f} km2, agrees", err=True)
    properties: dict[str, object] = {
        **place,
        "units": sorted(ids),
        "unit_count": len(ids),
        "area_km2": joined.area_m2 / 1e6,
        "stated_upstream_area_km2": None if stated is None else stated / 1e6,
        "holes": joined.holes,
        "holes_km2": joined.holes_m2 / 1e6,
        "network_ends": dict(net.ends),
    }
    if check_dem is not None and request is not None:
        properties["check_dem"] = _terrain_check(
            joined.polygon,
            geometries,
            net,
            crs,
            request,
            check_dem,
            cache,
            (threads, memory_budget, z_unit),
        )
    target.write_bytes(catchment_geojson(joined.polygon, crs, properties))
    typer.echo(f"{target}")


def _terrain_check(
    polygon: Polygon,
    units: Mapping[str, BaseGeometry],
    net: Network,
    crs: str,
    request: CatchmentRequest,
    dem: str,
    cache: Path | None,
    how: tuple[int | None, str | None, ZSpelling | None],
) -> dict[str, object]:
    """The terrain route from `request`'s seed on `dem`, compared with the
    network's `polygon`, said on stderr; its failure is reported, not refused."""
    res = cli_kit.cli_resources(how[0], how[1])
    try:
        result = delineate(request, _dem_repository([dem], cache, res, how[2]), res)
    except (OSError, ValueError, typer.BadParameter) as exc:
        text = exc.format_message() if isinstance(exc, typer.BadParameter) else str(exc)
        typer.echo(f"check: the terrain route failed, nothing compared: {text}", err=True)
        return {"dem": dem, "failed": text}
    terrain = shapely.transform(result.reduced, reprojector(result.crs, crs))
    check = compare(polygon, terrain, units, crs)
    outside = [
        {"id": u, "name": net.links[u].name, "area_km2": a} for u, a in check.units_mostly_outside
    ]
    total = sum(a for _, a in check.units_mostly_outside)
    typer.echo(
        f"check: terrain {check.terrain_area_km2:.2f} km2, shared {check.shared_km2:.2f}, "
        f"network only {check.network_only_km2:.2f}, terrain only {check.terrain_only_km2:.2f}; "
        f"{plural(len(outside), 'unit', 'units')} mostly outside the terrain catchment, "
        f"{total:.2f} km2",
        err=True,
    )
    return {
        "dem": dem,
        **asdict(check),
        "units_mostly_outside": outside,
        "units_mostly_outside_km2": total,
    }


@dataclass
class _DirectorySink:
    """`station-catchments`' `BatchSink`: a station's catchment file is written
    when its row arrives (the row holds the placement), and each row is said
    on stderr and written to `table` (results.csv) at once, flushed, so a run
    stopped by a bug keeps every finished row."""

    out: Path
    table: TextIO
    pending: dict[str, tuple[Station, Catchment]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        names = [f.name for f in fields(StationResult)]
        csv.writer(self.table).writerow(["class" if n == "station_class" else n for n in names])

    def catchment(self, station: Station, result: Catchment) -> None:
        self.pending[station.station] = (station, result)

    def row(self, row: StationResult) -> None:
        with _writing(self.out / "results.csv"):
            csv.writer(self.table).writerow([_cell(getattr(row, f.name)) for f in fields(row)])
            self.table.flush()
        if row.station in self.pending:
            station, result = self.pending.pop(row.station)
            names = [f.name for f in fields(StationResult)]
            gauge = names[names.index("placed_on") : names.index("causes") + 1]
            properties = {
                "station": station.station,
                "name": station.name,
                "series": list(station.series),
                **_outline_properties(result),
                **{k: getattr(row, k) for k in gauge},
            }
            # A station number is digits and dots (`Station`), so a safe name.
            path = self.out / f"{row.station}.geojson"
            with _writing(path):
                path.write_bytes(catchment_geojson(result.reduced, result.crs, properties))
        word: str = row.station_class or "well defined, no reference"
        if row.seeded_by == "lake":
            lake = row.lake_name or (None if row.lake_number is None else str(row.lake_number))
            word += f" (seeded by {'its lake' if lake is None else f'the lake {lake}'})"
        if row.station_class == "refused":
            word += f": {row.refusal_message}"
        elif row.causes:
            word += f" ({', '.join(row.causes)})"
        if row.reference_in_ours is not None and row.ours_in_reference is not None:
            word += (
                f"; reference in ours {100 * row.reference_in_ours:.1f} %, ours in reference "
                f"{100 * row.ours_in_reference:.1f} %, area ratio {row.area_ratio:.3f}"
            )
        typer.echo(
            f"{row.station} {row.name}: {word}" if row.name else f"{row.station}: {word}", err=True
        )


def _cell(value: object) -> object:
    """A table cell: a tuple joined by `;`, None empty, anything else as is."""
    if isinstance(value, tuple):
        return ";".join(map(str, value))
    return "" if value is None else value


def _read_beside[T](
    read: Callable[[Path], tuple[T, str]], path: Path, hint: str, what: str, crs: str
) -> T:
    """`read(path)`'s content, refused naming `hint` when it cannot be read or
    its CRS is not the river file's, `crs`."""
    try:
        content, file_crs = read(path)
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(str(exc), param_hint=hint) from exc
    if not same_crs(file_crs, crs):
        raise typer.BadParameter(
            f"the {what} file's CRS, {file_crs}, is not the river file's, {crs}; "
            "polygons are not reprojected",
            param_hint=hint,
        )
    return content


def station_catchments(
    dem: Annotated[list[str], typer.Option("--dem", help=f"{DEM_HELP}.")],
    stations: Annotated[
        Path, typer.Option("--stations", help="The stations (GeoJSON points, with a crs member).")
    ],
    rivers: Annotated[
        Path, typer.Option("--rivers", help="River lines (GeoJSON, in the DEM's CRS).")
    ],
    out_dir: Annotated[
        Path, typer.Option("--out-dir", writable=True, help="Where to write the files.")
    ],
    reference: Annotated[
        Path | None,
        typer.Option(
            "--reference",
            help="Reference catchment polygons (GeoJSON, in the river file's CRS), one per "
            "station number: each catchment is compared with its own and classed.",
        ),
    ] = None,
    lakes: Annotated[
        Path | None,
        typer.Option(
            "--lakes",
            help="The lakes.geojson that rasputin fetch-stations writes, or GeoJSON like it: "
            "polygons only, with a crs member naming the river file's CRS. A gauge in a "
            "lake, or on a lake line within 30 m of it, is seeded with the whole lake.",
        ),
    ] = None,
    map_radius: Annotated[
        float,
        typer.Option("--map-radius", metavar="METRES", help="How far a river line may be."),
    ] = 500.0,
    reach_up: Annotated[
        float,
        typer.Option(
            "--reach-up", metavar="METRES", help="River above each gauge burnt into the DEM."
        ),
    ] = 1000.0,
    only: Annotated[
        list[str] | None,
        typer.Option("--only", metavar="ID", help="Run only this station (repeat --only)."),
    ] = None,
    outline_tolerance: Annotated[
        float | None,
        typer.Option(
            "--outline-tolerance",
            metavar="METRES",
            help="Reduce each outline to within this of the fine one; default twice the cell.",
        ),
    ] = None,
    cache: CACHE = None,
    threads: THREADS = None,
    memory_budget: MEMORY_BUDGET = None,
    z_unit: Z_UNIT = None,
) -> None:
    """Write a catchment per station into --out-dir, as catchment --rivers does for one.

    Beside them: results.csv, a row per station, and summary.json. With
    --reference each catchment is compared with its reference polygon and
    classed: match, close, miss, uncertain or refused.
    """
    for name, value in (("--map-radius", map_radius), ("--reach-up", reach_up)):
        if not (math.isfinite(value) and value > 0.0):
            raise typer.BadParameter("must be a finite number above 0", param_hint=name)
    if outline_tolerance is not None and not (
        math.isfinite(outline_tolerance) and outline_tolerance >= 0.0
    ):
        raise typer.BadParameter("must be finite and at least 0", param_hint="--outline-tolerance")
    target = cli_kit.destination(out_dir, "--out-dir")
    try:
        station_list, stations_crs = read_stations(stations)
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(str(exc), param_hint="--stations") from exc
    unknown = set(only or ()) - {s.station for s in station_list}
    if unknown:
        raise typer.BadParameter(
            f"not in the stations file: {', '.join(sorted(unknown))}", param_hint="--only"
        )
    segments, crs, dropped = _segments(rivers)
    references = None
    if reference is not None:
        references = _read_beside(read_references, reference, "--reference", "reference", crs)
    lake_list = None
    if lakes is not None:
        lake_list = _read_beside(read_nve_lakes, lakes, "--lakes", "lakes", crs)
    res = cli_kit.cli_resources(threads, memory_budget)
    repository = _dem_repository(dem, cache, res, z_unit)
    try:
        geographic = repository.footprints()[0].meta.geographic
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(
            cli_kit.words(exc) if isinstance(exc, ValueError) else str(exc), param_hint="--dem"
        ) from exc
    if geographic:
        raise typer.BadParameter(
            "station-catchments needs a DEM in metres (a projected CRS); for a DEM in "
            "degrees, run catchment for each station",
            param_hint="--dem",
        )
    _reach_crs(crs, repository)
    request = BatchRequest(
        map_radius=map_radius,
        reach_up=reach_up,
        outline_tolerance=outline_tolerance,
        only=tuple(only or ()),
    )
    with _writing(target):
        target.mkdir(exist_ok=True)
    with _writing(target / "results.csv"):
        table = (target / "results.csv").open("w", encoding="utf-8", newline="")
    # Outside `table`, so the header write and the final flush and close are
    # refused too; `run_batch`'s own `OSError` is caught inside, as --dem's.
    with _writing(target / "results.csv"), table:
        sink = _DirectorySink(target, table)
        try:
            summary = asyncio.run(
                run_batch(
                    request,
                    repository,
                    station_list,
                    stations_crs,
                    segments,
                    crs,
                    references,
                    sink,
                    lakes=lake_list,
                    resources=res,
                )
            )
        except OSError as exc:
            raise typer.BadParameter(
                f"cannot read {exc.filename}: {exc}", param_hint="--dem"
            ) from exc
        except ValueError as exc:
            typer.echo(f"Error: {cli_kit.words(exc)}", err=True)
            raise typer.Exit(1) from exc
    dump = {**summary.model_dump(mode="json"), "river_copies_dropped": dropped}
    with _writing(target / "summary.json"):
        (target / "summary.json").write_text(json.dumps(dump, indent=2) + "\n", encoding="utf-8")
    counts = ", ".join(f"{k} {v}" for k, v in summary.classes.items())
    typer.echo(f"summary: {plural(summary.stations, 'station', 'stations')}: {counts}", err=True)
    if summary.known_refusals.line is not None:
        typer.echo(summary.known_refusals.line, err=True)
    typer.echo(f"{target}")


app.callback()(main)
for panel, commands in (
    (MESHING, (mesh, palette)),
    (CATCHMENTS, (catchment, network_catchment, station_catchments)),
    (DATA, (fetch, fetch_stations, sources, gml_to_gpkg)),
    (OTHER, (version,)),
):
    for command in commands:
        app.command(rich_help_panel=panel)(command)
app.command(rich_help_panel=OTHER, help=MCP_HELP)(mcp)
