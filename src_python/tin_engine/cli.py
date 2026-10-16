"""Command-line interface for the rasputin terrain engine.

This module is the single composition root: joining the file system to the
engine is its whole job. Paths from the command line also reach the domain,
DEM-input and repository modules, which read the files they name; ``viz/`` is
written against protocols and never names a core type, while the core never
sees a file, a path or a CRS. Everything that has to know both sides lives here.
"""

from __future__ import annotations

import asyncio
import csv
import datetime
import json
import math
import os
import shlex
import sys
import time
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Annotated, Final, TextIO, cast

import numpy as np
import numpy.typing as npt
import typer
from pydantic import ValidationError
from shapely.geometry.base import BaseGeometry

from tin_engine import edge_strip, final_check, installed_version
from tin_engine._core import (
    ChainRole,
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
from tin_engine.catchment import (
    Catchment,
    CatchmentRequest,
    GaugeResult,
    LakeError,
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
from tin_engine.chains import start_chains
from tin_engine.class_maps import derived_cleanup, run_vocabulary, shipped_names
from tin_engine.class_maps import load as load_class_map
from tin_engine.cli_kit import THREADS
from tin_engine.cli_kit import cli_resources as _resources
from tin_engine.cli_kit import destination as _destination
from tin_engine.cli_kit import words as _words
from tin_engine.crs import (
    crs_label,
    parse_crs,
    reprojector,
    same_crs,
    srs_entry,
    transform_description,
    transform_label,
)
from tin_engine.dem_input import (
    CachedSource,
    DemInput,
    DemRequest,
    OutCrsRequiredError,
    open_dem,
    repository_for,
)
from tin_engine.domain import DomainError, DomainPolygon
from tin_engine.elevation import Trimmed, trim
from tin_engine.feature_input import (
    ClassMap,
    FeatureError,
    FeatureRequest,
    FeatureSet,
    FeatureSource,
    open_features,
    read_lake_polygons,
)
from tin_engine.features import DEFAULT_VOCABULARY, EdgeVocabulary
from tin_engine.gauge import Gauge, Placement
from tin_engine.grid_domain import default_stride, refine_start_stride, subsample
from tin_engine.hydrography import RiverSegment, Station
from tin_engine.io.cog import NotCached
from tin_engine.io.domain_file import read_domain
from tin_engine.io.geojson import catchment_geojson
from tin_engine.io.geopackage_write import (
    Palette,
    create_geopackage,
    write_classes,
    write_land_cover,
    write_lines,
    write_run,
    write_styles,
    write_triangles,
)
from tin_engine.io.models import (
    BY_NAME,
    UNITS,
    Bounds,
    ClassFractions,
    DemTile,
    IndexWindow,
    RasterMeta,
    ZSpelling,
    ZUnit,
)
from tin_engine.io.ply import field_comments, write_ply
from tin_engine.io.repository import (
    CacheRepository,
    DemRepository,
    TiffDemRepository,
    create_sqlite,
)
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_nve_lakes, read_references, read_stations
from tin_engine.io.vtk_legacy import write_vtk
from tin_engine.land_cover_fractions import ClassRaster
from tin_engine.land_cover_fractions import assign as assign_fractions
from tin_engine.landcover import label_triangles
from tin_engine.mosaic import Seam
from tin_engine.palettes import paraview_preset
from tin_engine.raster import to_core
from tin_engine.resources import Resources
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
    RemoteSource,
    catalogue,
    check_values,
    data_root,
    user_source_file,
    xdg_base,
)
from tin_engine.stats import PhaseClock, Report, Sizes, quality, render
from tin_engine.target_grid import Block, TargetGrid
from tin_engine.tolerance_field import ToleranceLines, ToleranceSlope, line_segments
from tin_engine.viz.fixtures import GALLERY, Fixture
from tin_engine.viz.protocols import PslgLike
from tin_engine.viz.scene import build_scene
from tin_engine.viz.style import PropertyStroke, SvgStyle
from tin_engine.viz.svg import render_svg

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

#: The fixtures' vocabulary, mapped onto the enum. A role outside it is a
#: ``KeyError`` here rather than a silently unconstrained chain downstream.
ROLES = {
    "outer": ChainRole.Outer,
    "hole": ChainRole.Hole,
    "breakline": ChainRole.Breakline,
}

#: Which roles name chains that are rings, and therefore carry a closing edge no
#: graph stores. ``scene.py`` cannot work this out -- it may not name the enum --
#: and without it every ring closure is drawn as a false alarm.
#:
#: TWO vocabularies, because there are two kinds of source. A fixture's roles are
#: strings and a core graph's are ``ChainRole`` members, and the two compare
#: unequal; which one applies is a property of the source, so it is carried by
#: :class:`Attempt` beside it rather than chosen by a caller.
FIXTURE_CLOSED_ROLES: tuple[object, ...] = ("outer", "hole")
CORE_CLOSED_ROLES: tuple[object, ...] = (ChainRole.Outer, ChainRole.Hole)

#: The **gallery's** draw precedence, in DESCENDING priority: first match wins,
#: and an edge carrying several properties is drawn with exactly one token.
#:
#: Named by FEATURE rather than derived from the vocabulary's numbering, because
#: those are two different questions -- water is drawn over infrastructure
#: however a vocabulary chose to number them -- and named *here* because
#: `style.py` and `svg.py` hold no vocabulary and `viz/` may not import one at
#: all. It is this module's to supply, exactly as :data:`ROLES`,
#: :data:`FIXTURE_CLOSED_ROLES` and :data:`CORE_CLOSED_ROLES` are.
#:
#: It is the gallery's list and not all of :data:`DEFAULT_VOCABULARY`. The
#: legend is derived from it, and a token `svg.py`'s stylesheet has no rule for
#: draws identically to the row above it, so declaring all nine features would
#: put six legend rows on every picture that a reader cannot tell apart --
#: measured by drawing the gallery. It grows when the stylesheet does.
#: `coastline` is listed because `bridge-over-lake`'s shoreline is an area
#: feature: without an entry its token is never emitted, and the shoreline would
#: draw as a plain breakline. It sits between the other two, water still over
#: infrastructure, and `svg.py` has its CSS rule.
_PRECEDENCE = ("river", "coastline", "road")

_BIT_OF = {prop.name: prop.bit for prop in DEFAULT_VOCABULARY.properties}

PROPERTY_STROKES = tuple(PropertyStroke(bit=_BIT_OF[name], token=name) for name in _PRECEDENCE)

DEFAULT_LABEL_LIMIT = 500

#: Metres. The gallery's coordinates are UTM-shaped (``viz.fixtures.ORIGIN``), so
#: this is a millimetre on the ground.
#:
#: A policy choice made here on the record, because there is no default spacing
#: anywhere in C++ and three of the nine noder statuses point at this number --
#: two of them in opposite directions. Coarser merges features that are genuinely
#: distinct (a road four millimetres from a river becomes one edge carrying both
#: bits); finer runs out of lattice, ``kMaxGridIndex`` being ``2**51``, so at
#: spacing ``s`` the representable coordinate is ``s * 2.25e15`` and a UTM
#: easting overflows well before ``1e-12``. A millimetre sits four orders of
#: magnitude clear of that and shows the gallery what its input said.
DEFAULT_SNAP_SPACING = 1e-3
# The start mesh's minimum angle, degrees.
DEFAULT_START_MIN_ANGLE = 25.0
# The start mesh's minimum angle with land cover simplified, degrees.
LANDCOVER_START_MIN_ANGLE = 15.0
# The gain a start-quality point must bring, degrees; a negative gain adds every
# point.
DEFAULT_START_QUALITY_GAIN = 0.0
#: The four land-cover flags have their own help panel, so the long
#: --[no-]features-merge-same-class does not narrow the main one's names.
CLEAN_UP = "Land-cover clean-up"
#: The three clean-up distances default to values derived from the map.
DERIVED = "Default: from the class map's positional accuracy and the DEM cell (CORINE: {})."
#: Their flags and the bound each given value stays under.
OPTION_TOPS = (
    ("--features-repair", math.inf),
    ("--features-tolerance", math.inf),
    ("--features-outline-snap", 100.0),
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
    f"Needs the mcp extra: {_INSTALL_MCP_HELP}"
)


def mcp(
    jobs: Annotated[
        int, typer.Option("--jobs", min=1, help="Tool calls run at once; one beyond this waits.")
    ] = 1,
) -> None:
    """Serve the commands to an agent over the Model Context Protocol (``--help`` says how)."""
    try:
        from tin_engine.mcp import server
    except ImportError:
        typer.echo(f"Error: the mcp extra is not installed; run {INSTALL_MCP}", err=True)
        raise typer.Exit(2) from None
    from tin_engine import plugins  # the hook imports this module

    asyncio.run(server.serve(plugins.served or plugins.application(), jobs))


@dataclass(frozen=True, slots=True)
class Attempt:
    """What the engine made of a fixture, and what to draw it from.

    ``source`` is **the most-processed graph that exists**, with no branch on the
    CDT's status: no ``Pslg`` means the fixture, a ``Pslg`` but no ``NodedPslg``
    still means the fixture (a ``Pslg``'s indices *are* the fixture's, so nothing
    is gained by switching), and a ``NodedPslg`` means the noded graph whether or
    not the triangulation then succeeded -- when the noder succeeded and the
    backend did not, the noded graph is the truer picture, because it shows the
    splits the backend was actually handed.

    ``closed_roles`` is the vocabulary ``source``'s roles speak, and it is a
    field rather than a module constant because it is a property of the source:
    a fixture's roles are strings and a ``NodedPslg``'s are ``ChainRole``
    members, which compare unequal, so a source drawn with the other's
    vocabulary closes no ring and draws every ring edge as an alarm.

    Both are set together at every exit point below, which is the only form of
    this that a later edit cannot undo by accident.
    """

    source: PslgLike
    closed_roles: tuple[object, ...]
    mesh: IndexedMesh2 | None
    ok: bool
    status: str
    message: str


@dataclass(frozen=True, slots=True)
class _Run:
    """What ``build_pslg`` -> ``node`` -> ``triangulate`` made of some chains.

    ``noded`` is None when the validator or the noder refused; ``mesh`` is None
    when anything refused. Shared by the fixture path and the DEM path, so the
    DEM path reaches the engine without going through a ``Fixture``.
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

    ``clock`` records the three calls as ``start mesh: ...`` rows; ``draw``
    passes none and gets a fresh one that is discarded."""
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


def _triangulated(fixture: Fixture, delaunay: bool, spacing: float) -> Attempt:
    """Run the engine on a fixture and reduce what it said to one record.

    ``status`` and ``message`` are plain text for the header band. THREE distinct
    refusals reach this point and none of them is an exception: the validator can
    reject the constraint set, in which case there is no ``Pslg`` and the words
    come from every ``PslgDiagnostic`` rather than only the first; the noder can
    refuse to resolve it at this spacing; or the backend can refuse to
    triangulate what the noder produced.

    The two-part message is one string on purpose -- ``describe(status)`` and
    then the engine's own words, and nothing this module wrote. A sentence
    authored here would be a second authority for a C++ fact, drifting the first
    time a row is reworded.
    """
    chains = [
        ([int(i) for i in fixture.indices_of(c)], ROLES[chain.role], int(chain.properties))
        for c, chain in enumerate(fixture.chains)
    ]
    run = _engine(fixture.vertices, chains, delaunay, spacing)
    if run.noded is None:
        return Attempt(
            source=fixture,
            closed_roles=FIXTURE_CLOSED_ROLES,
            mesh=None,
            ok=False,
            status=run.status,
            message=run.message,
        )
    return Attempt(
        source=run.noded,
        closed_roles=CORE_CLOSED_ROLES,
        mesh=run.mesh,
        ok=run.ok,
        status=run.status,
        message=run.message,
    )


def _band(sentence: str, detail: str) -> str:
    """The engine's prose and the engine's detail, as ONE row of the band.

    One string rather than two, because ``svg.py`` joins rows with no separator
    and a sentence split across two would be printed with its words run
    together.
    """
    return " ".join(part for part in (sentence, detail) if part)


def _snap_spacing(value: float) -> float:
    """A usage error, exited before anything is drawn.

    NOT duplication of the engine's ``InvalidSnapSpacing``: the two answer
    different questions to different callers. This is a mis-typed option, and
    ``typer`` refuses it the same way a bad ``--out`` is refused, without drawing
    a picture of a failure that is not about the terrain. The C++ check is a
    library precondition with callers that are not this CLI, and neither can be
    deleted in favour of the other.
    """
    if not math.isfinite(value) or value <= 0.0:
        raise typer.BadParameter(f"snap-spacing must be finite and positive, got {value}")
    return value


def draw(
    name: Annotated[str | None, typer.Argument(help="Gallery fixture to draw.")] = None,
    list_: Annotated[bool, typer.Option("--list", help="List the gallery and exit.")] = False,
    out: Annotated[Path | None, typer.Option("--out", help="Where to write the SVG.")] = None,
    out_parent: Annotated[
        Path | None,
        typer.Option("--out-parent", help="Refuse any --out that resolves outside this."),
    ] = None,
    title: Annotated[str, typer.Option("--title", help="Free text for the header band.")] = "",
    labels: Annotated[bool, typer.Option("--labels", help="Draw vertex indices.")] = False,
    label_limit: Annotated[
        int, typer.Option("--label-limit", help="Refuse --labels above this many vertices.")
    ] = DEFAULT_LABEL_LIMIT,
    show_vertices: Annotated[
        bool, typer.Option("--vertices", help="Draw a dot per vertex.")
    ] = False,
    delaunay: Annotated[
        bool,
        typer.Option(
            "--delaunay/--no-delaunay",
            help="Triangulate with or without the Delaunay property, to compare the two.",
        ),
    ] = True,
    snap_spacing: Annotated[
        float,
        typer.Option(
            "--snap-spacing",
            callback=_snap_spacing,
            help="Snap grid spacing for the noder, in the input's units.",
        ),
    ] = DEFAULT_SNAP_SPACING,
) -> None:
    """Draw a gallery fixture as an SVG.

    A fixture the engine refuses is drawn too, with the engine's status and
    words in the header band.
    """
    if list_:
        for fixture in GALLERY.values():
            typer.echo(f"{fixture.name}: {fixture.description}")
        return
    if name is None:
        raise typer.BadParameter("no fixture named; run 'rasputin draw --list' to see the gallery")
    if name not in GALLERY:
        raise typer.BadParameter(f"unknown fixture {name}; the gallery is: {', '.join(GALLERY)}")

    fixture = GALLERY[name]
    attempt = _triangulated(fixture, delaunay=delaunay, spacing=snap_spacing)
    scene = build_scene(
        attempt.source, attempt.mesh, ok=attempt.ok, closed_roles=attempt.closed_roles
    )
    if labels and len(scene.vertices) > label_limit:
        raise typer.BadParameter(
            f"{len(scene.vertices)} vertices exceeds the label limit of {label_limit}; "
            f"raise --label-limit or drop --labels"
        )

    target = _destination(out, out_parent, name)
    document = render_svg(
        scene,
        SvgStyle(show_vertices=show_vertices, property_strokes=PROPERTY_STROKES),
        title=title,
        status=attempt.status,
        message=attempt.message,
        labels=labels,
    )
    target.write_text(document, encoding="utf-8")
    typer.echo(f"{target}")


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

    A pair reached by two chains takes the UNION, as ``viz.scene`` does: the
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

    Computed from ``constrained_edges`` and the noded graph alone, never from
    ``viz.scene``: that module is the renderer, and a mesh writer reaching into
    it for a topology join would make ``viz/`` load-bearing for a path that
    draws nothing.

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
        typer.Option("--out", help="Where to write: .vtk for ParaView, .ply or .gpkg for QGIS."),
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
            help="Mesh only XMIN YMIN XMAX YMAX in the DEM's CRS, or in --out-crs's when "
            "given, snapped outward to the DEM's nodes.",
        ),
    ] = None,
    stride: Annotated[
        int | None,
        typer.Option("--stride", help="Every Nth DEM node. Default: <= 256 a side."),
    ] = None,
    out_edges: Annotated[
        Path | None,
        typer.Option("--out-edges", help="With .ply, also write the constraint edges as a PLY."),
    ] = None,
    out_parent: Annotated[
        Path | None,
        typer.Option("--out-parent", help="Refuse any output path resolving outside this."),
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
            help="Mesh in this CRS: a DEM in another CRS (required for a "
            "geographic one) is resampled onto a square grid here, and with --tolerance "
            "checked against its own nodes. --bbox is then in this CRS.",
        ),
    ] = None,
    start_min_angle: Annotated[
        float | None,
        typer.Option(
            "--start-min-angle",
            help="With --tolerance, add DEM nodes to the start mesh until its triangles "
            "have no angle under this many degrees, or a stated reason prevents it; "
            "0 is off, at most 35. Default: 25, or 15 with land cover simplified "
            "(--features-tolerance above 0).",
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
    features: Annotated[
        list[Path] | None,
        typer.Option(
            "--features",
            help="With --domain, polygons and lines as constraints (.geojson, .json, .gpkg "
            "or .gml), clipped to the domain, each edge carrying its class map's bits. "
            "Repeatable: each --features is one source, paired by position with the Nth "
            "of --features-crs, --features-layer and --features-map.",
        ),
    ] = None,
    features_crs: Annotated[
        list[str] | None,
        typer.Option("--features-crs", help="The Nth --features file's CRS, any pyproj reads."),
    ] = None,
    features_layer: Annotated[
        list[str] | None,
        typer.Option("--features-layer", help="The Nth --features GeoPackage's features table."),
    ] = None,
    features_map: Annotated[
        list[str] | None,
        typer.Option(
            "--features-map",
            help="The Nth --features' map from attributes to edge bits: "
            f"{', '.join([*shipped_names(), 'property'])}, or a class-map file (.toml). "
            "Default: property.",
        ),
    ] = None,
    features_repair: Annotated[
        float | None,
        typer.Option(
            "--features-repair",
            rich_help_panel=CLEAN_UP,
            show_default=False,
            help="Metres: land-cover borders of one file closer than this become one "
            "border, and gaps narrower than it are filled; 0 is off. " + DERIVED.format(1),
        ),
    ] = None,
    features_merge_same_class: Annotated[
        bool,
        typer.Option(
            help="Drop borders between neighbouring land-cover polygons of one class. Default: on.",
            rich_help_panel=CLEAN_UP,
            show_default=False,
        ),
    ] = True,
    features_tolerance: Annotated[
        float | None,
        typer.Option(
            "--features-tolerance",
            rich_help_panel=CLEAN_UP,
            show_default=False,
            help="Metres: simplify land-cover borders, each moved at most this far from its "
            "border after the repair and the outline rule, each class keeping its area; the "
            "simplification can bring a border back within the outline-snap distance of the "
            "outline, but not closer than the repair distance. 0 is off. " + DERIVED.format(50),
        ),
    ] = None,
    features_outline_snap: Annotated[
        float | None,
        typer.Option(
            "--features-outline-snap",
            rich_help_panel=CLEAN_UP,
            show_default=False,
            help="Metres, under 100: move land-cover borders this close to the domain "
            "outline onto it; 0 is off. " + DERIVED.format(5),
        ),
    ] = None,
    land_cover: Annotated[
        str | None,
        typer.Option(
            "--land-cover",
            help="The exact share of each land-cover class in every triangle, "
            "from a cached catalogue source of classes "
            f"({', '.join(k for k, s in SOURCES.items() if s.values == 'classes')}) or a uint8 "
            "GeoTIFF, in the mesh file.",
        ),
    ] = None,
    land_cover_year: Annotated[
        int | None,
        typer.Option(
            "--land-cover-year",
            help="With --land-cover, the year of the classes. Default: the source's last year.",
        ),
    ] = None,
    land_cover_map: Annotated[
        str | None,
        typer.Option(
            "--land-cover-map",
            help="With --land-cover, a class map with codes naming its classes in a .gpkg, as "
            "--features-map. Default: the source's.",
        ),
    ] = None,
    threads: THREADS = None,
    memory_budget: MEMORY_BUDGET = None,
    z_unit: Z_UNIT = None,
    stats: Annotated[
        str | None,
        typer.Option(
            "--stats",
            help="Also write sizes, quality and timings as Markdown to this path (.md "
            "recommended); - prints it on stdout after the path line(s).",
        ),
    ] = None,
    record_path: Annotated[
        str | None,
        typer.Option(
            "--record",
            help="Also write the run's record as JSON to this path, for reproducibility.",
        ),
    ] = None,
) -> None:
    """Mesh a DEM into a TIN and write it as .vtk, .ply or .gpkg, chosen by the suffix of --out.

    The DEM is one GeoTIFF, several, a directory of tiles, or a source: the
    package's or one in your sources file (see rasputin sources). --bbox or
    --domain cuts it. Without --tolerance every --stride-th DEM node is
    triangulated; with it the mesh is refined until every triangle is within
    that many metres of the DEM. --features adds polygons and lines as
    constraint edges; a class map with codes, or --land-cover, gives every
    triangle a land_cover_code, coloured by rasputin palette.

    .vtk is one file for ParaView. .ply is the surface for QGIS, with the
    constraint lines in a second file given by --out-edges. .gpkg is a GeoPackage
    with triangle and line layers. Each carries what the run did; --stats and
    --record write it in full. A refused run exits non-zero and writes no file.
    """
    clock = PhaseClock()
    res = _resources(threads, memory_budget)
    codes: npt.NDArray[np.int32] | None = None
    codes_text = ""
    table: ClassFractions | None = None
    palette: Palette | None = None
    cover: list[tuple[int, BaseGeometry]] = []
    seams: tuple[Seam, ...] = ()
    feature_paths = features or []
    feature_crss = features_crs or []
    feature_layers = features_layer or []
    feature_maps = features_map or []
    # One length rule per paired option. A list longer than --features is
    # refused naming the offending flag and the counts (none given with --features
    # absent included); a shorter list takes per-source defaults.
    for flag, value in (
        ("--features-crs", feature_crss),
        ("--features-layer", feature_layers),
        ("--features-map", feature_maps),
    ):
        if len(value) > len(feature_paths):
            raise typer.BadParameter(
                f"{len(value)} {flag} given for {len(feature_paths)} --features; "
                "each pairs with one source by position",
                param_hint=flag,
            )
    if feature_paths and domain is None:
        raise typer.BadParameter("needs --dem, --domain and --tolerance", param_hint="--features")
    # Every map is read here, before the DEM, so every file error comes first.
    maps: list[ClassMap] = []
    where = "--features"
    try:
        for i, path in enumerate(feature_paths):
            where = f"--features {i + 1} ({path.name})"
            maps.append(load_class_map(feature_maps[i] if i < len(feature_maps) else "property"))
        vocabulary = run_vocabulary(maps)
    except ValueError as exc:
        raise typer.BadParameter(f"{where}: {exc}", param_hint="--features-map") from exc
    known = _known_sources(os.environ)
    classes = _class_source(
        land_cover, land_cover_year, land_cover_map, cache_root(cache, os.environ), known
    )
    if classes is not None:
        for coded_map in (m for m in maps if m.codes):
            raise typer.BadParameter(
                f"one mesh has one land_cover_code system, and --features-map {coded_map.name} "
                "already gives every triangle a code",
                param_hint="--land-cover",
            )
    overrides = {
        "repair_m": features_repair,
        "tolerance_m": features_tolerance,
        "outline_snap_m": features_outline_snap,
    }
    for (flag, top), metres in zip(OPTION_TOPS, overrides.values(), strict=True):
        if metres is not None and not (math.isfinite(metres) and 0 <= metres < top):
            under = "" if top == math.inf else f" and under {top:g}"
            raise typer.BadParameter(f"must be finite, >= 0{under}, got {metres}", param_hint=flag)
    if domain is None and domain_crs is not None:
        raise typer.BadParameter("applies only with --domain", param_hint="--domain-crs")
    lines = _tolerance_lines(tolerance, tolerance_near, tolerance_near_crs, tolerance_ramp)
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
    if domain is not None and tolerance is None:
        raise typer.BadParameter("--domain needs --tolerance", param_hint="--tolerance")
    if start_min_angle is not None and not (
        math.isfinite(start_min_angle) and 0 <= start_min_angle <= 35
    ):
        raise typer.BadParameter(
            f"must be finite, >= 0 and <= 35, got {start_min_angle}",
            param_hint="--start-min-angle",
        )
    if start_min_angle is not None and tolerance is None:
        raise typer.BadParameter(
            "--start-min-angle needs --tolerance", param_hint="--start-min-angle"
        )
    if start_quality_gain is not None and not (
        math.isfinite(start_quality_gain) and start_quality_gain <= 10
    ):
        raise typer.BadParameter(
            f"must be finite and <= 10, got {start_quality_gain}",
            param_hint="--start-quality-gain",
        )
    if start_quality_gain is not None and tolerance is None:
        raise typer.BadParameter(
            "--start-quality-gain needs --tolerance", param_hint="--start-quality-gain"
        )
    if no_constraint_feet and tolerance is None:
        raise typer.BadParameter(
            "--no-constraint-feet needs --tolerance", param_hint="--no-constraint-feet"
        )
    paths, cached, entry = _dem_sources(dem, cache_root(cache, os.environ), known)
    caller = _z_unit(z_unit) or (entry.z_unit if entry else None)
    given = None
    if domain is not None:
        with clock.phase("domain read"):
            try:
                given = read_domain(domain, domain_crs)
            except DomainError as exc:
                raise typer.BadParameter(str(exc), param_hint="--domain") from exc
    try:
        opened = _open_dem(paths, bbox, clock, given, cached, out_crs, res, caller)
    except NotCached as exc:
        area = f"--domain {domain}" if domain is not None else ""
        if bbox is not None:
            area = "--bbox " + " ".join(repr(v).removesuffix(".0") for v in bbox)
        frame = f"--out-crs {shlex.quote(out_crs)}" if out_crs else ""
        flags = " ".join(filter(None, (area, frame, f"--cache {cache}" if cache else "")))
        run = f"rasputin fetch {dem[0]} {flags}".rstrip()
        raise typer.BadParameter(f"{exc}; run: {run}", param_hint="--dem") from exc
    # What is read later comes out of `opened`, which then goes, so the
    # tile's one reference is the holder `_dem_mesh` empties.
    label, plan, grid, source_crs = opened.label, opened.plan, opened.grid, opened.source_crs
    dem_domain, checks, opened_seams = opened.domain, opened.checks, opened.seams
    held, dem_crs = _Once(opened.tile), opened.tile.meta.crs
    # The clean-up scale is the coarsest coded map's accuracy (CORINE's
    # 100 m without one), never finer than a DEM cell.
    cell = max(opened.tile.meta.delta_x, abs(opened.tile.meta.delta_y))
    del opened
    accuracy = max((m.accuracy_m for m in maps if m.codes), default=100.0)
    derived = derived_cleanup(accuracy, cell).model_dump()
    cleanup = {k: derived[k] if v is None else v for k, v in overrides.items()}
    cleanup["merge_same_class"] = features_merge_same_class
    found = None
    sources: tuple[FeatureSource, ...] = ()
    if feature_paths and dem_domain is not None:
        sources = _feature_sources(feature_paths, feature_crss, feature_layers, maps)
        found = _open_features(
            sources, dem_domain, dem_crs, clock, cleanup, vocabulary, res.threads, cell
        )
    # 15 degrees with land cover simplified, unless given.
    simplified = (
        found is not None and found.cover_vertices is not None and cleanup["tolerance_m"] > 0
    )
    angle = LANDCOVER_START_MIN_ANGLE if simplified else DEFAULT_START_MIN_ANGLE
    angle = angle if start_min_angle is None else start_min_angle
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
        lines,
        steep,
        res.threads,
        vocabulary,
    )
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
        seams = opened_seams
        values["dem_tiles"] = escaped_ascii("; ".join(names))
        values["dem_seams"] = escaped_ascii("; ".join(s.entry() for s in seams)) or SEAMS_AGREE
    if given is not None:
        how = transform_label(given.crs, dem_crs)
        values |= {"domain_crs": crs_label(given.crs), "domain_transform": how}
    if found is not None and sources and dem_run.feature_counts is not None:
        # One entry per source, joined. Lines and their vertices are a
        # property of the merged PSLG, not of a source, so every entry
        # carries the run totals; the per-source part is the name, layer,
        # map and feature count (found.counts[i]).
        chains, feature_vertices = dem_run.feature_counts
        totals = f"{plural(chains, 'line', 'lines')}, "
        totals += plural(feature_vertices, "vertex", "vertices")
        texts, crs_texts, transforms, notices = [], [], [], []
        for i, src in enumerate(sources):
            layer = f" layer {found.layers[i]}" if found.layers[i] else ""
            counted = plural(found.counts[i], "feature", "features")
            texts.append(
                f"{src.path.name}{layer}, class map {src.class_map.name}: {counted}, {totals}"
            )
            own = found.crs[i]
            crs_texts.append(crs_label(own))
            transforms.append(transform_label(own, dem_crs))
            if src.class_map.notice and src.class_map.notice not in notices:
                notices.append(src.class_map.notice)
        values["features"] = escaped_ascii("; ".join(texts))
        values["features_crs"] = "; ".join(crs_texts)
        values["features_transform"] = "; ".join(transforms)
        values["features_notice"] = "; ".join(notices) or None
        licences = dict.fromkeys(s.class_map.licence for s in sources if s.class_map.licence)
        values["features_licence"] = escaped_ascii("; ".join(licences)) or None
        values |= {f"features_{k}": v for k, v in cleanup.items()}
        values["features_merge_same_class"] = "on" if features_merge_same_class else "off"
        if found.cover_vertices is not None:
            before, after = found.cover_vertices
            values["land_cover_vertices"] = f"{before} after the clip, {after} after clean-up"
            values["land_cover_area_moved_m2"] = f"{found.area_changed:.1f}"
        # Labelling runs iff a source carries codes (one code system per
        # mesh is enforced in _feature_sources). Every coded source's
        # polygons are labelled together over the merged FeatureSet.
        coded = next((s.class_map for s in sources if s.class_map.codes), None)
        if coded is not None:
            codes, codes_text = _land_cover(dem_run.trimmed, found, coded, snap_spacing, clock)
            values["land_cover_codes"] = codes_text
            palette = coded.palette
            cover = [(f.code, f.polygon) for f in found.features if f.polygon and f.code]
    if classes is not None:
        table, codes, more = _class_fractions(
            classes, dem_run.trimmed, dem_crs, clock, domain, out_crs, cache, res.threads
        )
        codes_text = str(more["land_cover_codes"])
        if classes.class_map is not None:  # names and colours for the .gpkg
            palette = classes.class_map.palette
        values |= more
    build = stride_record if tolerance is None else refined_record
    record = build(triangles=len(surface_mesh.triangles), threads=res.threads, **values)

    vertices = surface_mesh.vertices
    surface = _destination(out, out_parent, label)
    targets = [surface]
    if out_edges is not None:
        # Resolved, because two spellings of one path are still one file. Both
        # writes succeed, the second overwrites the first, the command echoes
        # two paths and exits 0 -- the caller has lost the surface they asked
        # for and nothing said so.
        constraints = _destination(out_edges, out_parent, label, "--out-edges")
        if constraints == surface:
            raise typer.BadParameter(
                f"--out and --out-edges both resolve to {surface}; "
                "the second would overwrite the first",
                param_hint="--out-edges",
            )
        targets.append(constraints)
    report_target = _report_target(stats, out_parent, label, targets)
    record_target = _record_target(record_path, out_parent, label, [*targets, report_target])
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
            _write_geopackage(surface, surface_mesh, codes, palette, cover, record, vocabulary)
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
) -> None:
    """``target`` as a GeoPackage, written as ``<target>.part``
    and renamed onto it: an existing file is replaced, never appended to, and
    a failed run leaves nothing under the asked name."""
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
class _ClassSource:
    """``--land-cover`` as checked: a catalogue entry with its year (None for a
    source without years) and cache root, or a file; and the class map that
    names its classes, if any."""

    source: RemoteSource | None
    year: int | None
    root: Path | None
    path: Path | None
    class_map: ClassMap | None


def _class_source(
    value: str | None,
    year: int | None,
    map_value: str | None,
    root: Path | None,
    known: Mapping[str, RemoteSource],
) -> _ClassSource | None:
    """``--land-cover``, ``--land-cover-year`` and ``--land-cover-map`` checked
    before any work: a catalogue key is that cached source, anything else a
    file, whose header is read here. Every refusal is a usage error naming
    the option."""
    if value is None:
        for given, flag in ((year, "--land-cover-year"), (map_value, "--land-cover-map")):
            if given is not None:
                raise typer.BadParameter("needs --land-cover", param_hint=flag)
        return None
    source, path, chosen = known.get(value), None, None
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
            raise typer.BadParameter(
                NO_CACHE,
                param_hint="--land-cover",
            )
    else:
        if year is not None:
            raise typer.BadParameter("a file has no years", param_hint="--land-cover-year")
        path = Path(value)
        if not path.is_file():
            raise typer.BadParameter(
                f"{value} is neither a catalogue source ({', '.join(known)}) nor a file",
                param_hint="--land-cover",
            )
        try:  # the header and one cell: a class raster in uint8 with a CRS
            corner = IndexWindow(row0=0, col0=0, rows=1, cols=1)
            TiffDemRepository([path], z_unit="metre").load_class_window(path.name, corner)
        except (OSError, ValueError) as exc:
            raise typer.BadParameter(str(exc), param_hint="--land-cover") from exc
    name = map_value or (source.class_map if source is not None else "")
    try:
        class_map = load_class_map(name) if name else None
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--land-cover-map") from exc
    if class_map is not None and not class_map.codes:
        raise typer.BadParameter(
            f"{name} names no class codes; a land-cover raster holds codes",
            param_hint="--land-cover-map",
        )
    return _ClassSource(source, chosen, None if path else root, path, class_map)


def _class_fractions(
    classes: _ClassSource,
    trimmed: Trimmed,
    dem_crs: str,
    clock: PhaseClock,
    domain: Path | None,
    out_crs: str | None,
    cache: Path | None,
    threads: int = 1,
) -> tuple[ClassFractions, npt.NDArray[np.int32], dict[str, Value]]:
    """The exact class fractions from the cached year or the file, timed as
    the phase ``land cover fractions``; the table, the codes and the record
    entries. An object the cache lacks is refused naming the ``fetch`` that
    gets it."""
    source, year = classes.source, classes.year
    repository: ClassRaster
    if source is None:
        assert classes.path is not None
        name = classes.path.name
        repository = TiffDemRepository([classes.path], z_unit="metre", threads=threads)
    else:
        assert classes.root is not None
        cached_classes = CacheRepository(
            classes.root, source.id, nodata=source.nodata, threads=threads, z_unit="metre"
        )
        repository = cached_classes
        url = source.url_for(year)
        name = url.rsplit("/", 1)[-1].split("?")[0].rsplit(".", 1)[0]
        try:
            cached = {f.name for f in cached_classes.footprints()}
        except NotCached:
            cached = set()
        if name not in cached:
            words = [
                f"rasputin fetch {source.id}",
                f"--year {year}" if year is not None else "",
                f"--domain {domain}" if domain is not None else "",
                f"--out-crs {shlex.quote(out_crs)}" if out_crs else "",
                f"--cache {cache}" if cache else "",
            ]
            run = " ".join(filter(None, words))
            raise typer.BadParameter(
                f"{source.id}/{name} is not in the cache; run: {run}", param_hint="--land-cover"
            )
    with clock.phase("land cover fractions"):
        try:
            assigned = assign_fractions(
                trimmed.vertices,
                trimmed.triangles,
                mesh_crs=dem_crs,
                repository=repository,
                name=name,
                threads=threads,
            )
        except (ValueError, NotCached) as exc:
            raise typer.BadParameter(_words(exc), param_hint="--land-cover") from exc
    typer.echo(
        f"land cover: {len(assigned.table.classes)} class shares over "
        f"{assigned.table.triangles} triangles, {len(np.unique(assigned.table.classes))} classes",
        err=True,
    )
    label = escaped_ascii(name) if source is None else source.id
    values: dict[str, Value] = {
        "land_cover_codes": f"{label} class number; 0 = no data, and every constraint line",
        "land_cover_source": label if year is None else f"{label}, {year}",
    }
    if source is not None:  # a file carries no credit, licence or citation
        fetched = cached_classes.fetched(name) or datetime.date.today()
        values |= {
            "land_cover_credit": escaped_ascii(source.credit_for(url=url, date=fetched)),
            "land_cover_licence_note": escaped_ascii(source.licence_note),
            "land_cover_cite": escaped_ascii("; ".join(source.cite)) or None,
        }
    return assigned.table, assigned.codes, values


def _land_cover(
    trimmed: Trimmed, found: FeatureSet, cmap: ClassMap, spacing: float, clock: PhaseClock
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
    text = f"{cmap.codes}, attribute {cmap.attribute}, map {cmap.name}; "
    return labels.codes, text + "0 = in no polygon, and every constraint line"


def palette(
    name: Annotated[
        str,
        typer.Argument(help=f"A class map: {', '.join(shipped_names())}, or a .toml file."),
    ],
    out: Annotated[
        Path | None, typer.Option("--out", help="Write the preset here, not to stdout.")
    ] = None,
) -> None:
    """Write a ParaView colour preset for a mesh's land_cover_code.

    Import it in ParaView's Colour Map Editor (Choose Preset, Import), then
    colour by land_cover_code; tick Interpret Values As Categories if the
    preset does not.
    """
    try:
        cmap = load_class_map(name)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    if cmap.palette is None:
        raise typer.BadParameter(f"{name} has no classes to colour")
    text = json.dumps(paraview_preset(cmap.palette, cmap.preset), indent=1) + "\n"
    if out is None:
        typer.echo(text, nl=False)
        return
    out.write_text(text)
    typer.echo(f"{out}")


def _report_target(
    stats: str | None, out_parent: Path | None, label: str, meshes: list[Path]
) -> Path | None:
    """Where ``--stats`` writes, None for stdout (``-``) or no report,
    checked before any mesh file is written."""
    if stats is None or stats == "-":
        return None
    target = _destination(Path(stats), out_parent, label, "--stats")
    if target in meshes:
        raise typer.BadParameter(
            f"resolves to {target}; the report would overwrite the mesh", param_hint="--stats"
        )
    return target


def _record_target(
    record: str | None, out_parent: Path | None, label: str, taken: list[Path | None]
) -> Path | None:
    """Where ``--record`` writes, checked before any file is written; never
    standard output (``--stats -``'s) nor a file this run also writes."""
    if record is None:
        return None
    if record == "-":
        raise typer.BadParameter("standard output is --stats -'s", param_hint="--record")
    target = _destination(Path(record), out_parent, label, "--record")
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
        typer.Option(
            "--bbox", metavar="BOX", help="Or XMIN YMIN XMAX YMAX in --out-crs, else the source's."
        ),
    ] = None,
    out_crs: Annotated[
        str | None, typer.Option("--out-crs", help="The frame the mesh's box is in.")
    ] = None,
    cache: Annotated[
        Path | None, typer.Option("--cache", help=f"The tile cache. {CACHE_DEFAULT}")
    ] = None,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Read the headers, print the plan, write nothing.")
    ] = False,
    refresh: Annotated[
        bool, typer.Option("--refresh", help="Discard the source's cache and fetch anew.")
    ] = False,
    connections: Annotated[int, typer.Option("--connections", min=1)] = 8,
    year: Annotated[
        int | None,
        typer.Option(
            "--year", help="For a source with one file per year: which. Default: the last."
        ),
    ] = None,
) -> None:
    """Copy what a mesh of a domain reads from a remote source into the tile cache.

    Then rasputin mesh --dem SOURCE runs offline. A source is the package's or
    one in your sources file (see rasputin sources).
    """
    import asyncio

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
        corners = ("x_min", "y_min", "x_max", "y_max")
        box = None if bbox is None else Bounds(**dict(zip(corners, bbox, strict=True)))
        request = FetchRequest(
            source=source, domain=given, box=box, out_crs=out_crs, connections=connections,
            dry_run=dry_run, refresh=refresh, year=year,
        )  # fmt: skip
    except DomainError as exc:
        raise typer.BadParameter(str(exc), param_hint="--domain") from exc
    except ValueError as exc:
        raise typer.BadParameter(_words(exc)) from exc
    tenths = [0]

    def progress(done: int, total: int) -> None:
        if total and done * 10 // total > tenths[0]:
            tenths[0] = done * 10 // total
            typer.echo(f"{done:,} of {total:,} bytes", err=True)

    client, writer = RangeClient(), CacheWriter(root, source)
    try:
        report = asyncio.run(run(request, known[source], client, writer, progress=progress))
    except (FetchError, CacheError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1) from exc
    for plan in report.plans:
        typer.echo(f"{plan.object_id}: {len(plan.blocks):,} blocks, {plan.bytes:,} bytes to fetch")
    head = "would fetch" if dry_run else "fetched"
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
    out_dir: Annotated[Path, typer.Option("--out-dir", help="Where to write the files.")],
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
            help="WEST SOUTH EAST NORTH in degrees (WGS 84): report the share of it each "
            "source covers. A point is a box with WEST = EAST and SOUTH = NORTH.",
        ),
    ] = None,
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

    if bbox is not None and (bbox[0] > bbox[2] or bbox[1] > bbox[3]):
        raise typer.BadParameter("WEST must not exceed EAST, nor SOUTH NORTH", param_hint="--bbox")
    typer.echo(listing(os.environ, cache_root(cache, os.environ), bbox).model_dump_json())


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
    (`user_source_file`); a file it refuses is a usage error naming it."""
    try:
        return catalogue(user_source_file(environ)).sources
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


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
        return repository_for(paths, cached=cached, threads=res.threads, z_unit=caller)[0]
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(
            _words(exc) if isinstance(exc, ValueError) else str(exc), param_hint="--dem"
        ) from exc


def _open_dem(
    dem: tuple[Path, ...],
    bbox: tuple[float, float, float, float] | None,
    clock: PhaseClock,
    domain: DomainPolygon | None = None,
    cached: CachedSource | None = None,
    target_crs: str | None = None,
    res: Resources | None = None,
    z_unit: ZUnit | None = None,
) -> DemInput:
    """``--dem`` and ``--bbox`` or the read ``--domain`` to one tile; every
    refusal, the reader's, the mosaic's or the domain's extent, is a usage
    error in its own words."""
    try:
        bounds = (
            None
            if bbox is None
            else Bounds(**dict(zip(("x_min", "y_min", "x_max", "y_max"), bbox, strict=True)))
        )
    except ValidationError as exc:
        raise typer.BadParameter(_words(exc), param_hint="--bbox") from exc
    try:
        with clock.phase("decode"):
            request = DemRequest(
                sources=dem,
                cached=cached,
                bounds=bounds,
                domain=domain,
                target_crs=target_crs,
                z_unit=z_unit,
            )
            return open_dem(request, res)
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
    except ValueError as exc:
        raise typer.BadParameter(_words(exc), param_hint="--dem") from exc


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


def _feature_sources(
    paths: list[Path], crss: list[str], layers: list[str], maps: list[ClassMap]
) -> tuple[FeatureSource, ...]:
    """The repeatable ``--features`` flags as one source per path, paired
    by position; a shorter paired list defaults its missing entries. A flag that
    cannot apply is a usage error naming the 1-based index and path, because with
    several sources the flag alone no longer identifies which one."""
    built: list[FeatureSource] = []
    for i, path in enumerate(paths):
        layer = layers[i] if i < len(layers) else None
        if layer is not None and path.suffix.lower() != ".gpkg":
            raise typer.BadParameter(
                f"applies to a .gpkg only, not --features {i + 1} ({path.name})",
                param_hint="--features-layer",
            )
        crs = crss[i] if i < len(crss) else None
        built.append(FeatureSource(path=path, class_map=maps[i], layer=layer, crs=crs))
    # One mesh has one land_cover_code system. Two sources whose maps name
    # different code systems (any file can name its own) are refused before
    # any meshing.
    systems = {s.class_map.codes for s in built if s.class_map.codes}
    if len(systems) > 1:
        raise typer.BadParameter(
            f"sources carry different code systems ({', '.join(sorted(systems))}); "
            "one mesh has one land_cover_code system",
            param_hint="--features-map",
        )
    return tuple(built)


def _open_features(
    sources: tuple[FeatureSource, ...],
    domain: DomainPolygon,
    dem_crs: str,
    clock: PhaseClock,
    cleanup: dict[str, float],
    vocabulary: EdgeVocabulary,
    threads: int = 1,
    cell: float = 0.0,
) -> FeatureSet:
    """Read and clip every source into one set, timed as
    ``features read`` and ``features clip``; a refusal is a usage error naming
    the feature."""
    t0 = time.perf_counter()
    try:
        request = FeatureRequest(sources=sources, vocabulary=vocabulary, threads=threads, **cleanup)
        found = open_features(request, domain, dem_crs, cell=cell)
    except FeatureError as exc:
        raise typer.BadParameter(str(exc), param_hint="--features") from exc
    clock.add("features read", time.perf_counter() - t0 - found.clip_seconds)
    clock.add("features clip", found.clip_seconds)
    if found.cover_vertices is not None:
        clock.add("features clip: clean-up", found.cleanup_seconds)
    typer.echo(
        f"features: {len(found.features)} kept ({found.clipped} cut at the domain outline), "
        f"{found.outside} outside the domain, {found.empty} empty",
        err=True,
    )
    for table in found.scanned:
        typer.echo(f"{table}: this layer has no spatial index, so every row was read", err=True)
    return found


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
    (``run_record``'s builder arguments)."""

    trimmed: Trimmed
    meta: RasterMeta
    domain_vertices: int | None
    domain_holes: int | None
    start_vertices: int | None
    start_triangles: int | None
    values: Mapping[str, Value]
    feature_counts: tuple[int, int] | None = None


class _Once[T]:
    """A value handed on once: `take` returns it and forgets it,
    so the holder keeps nothing alive after; a second `take` is a bug."""

    def __init__(self, value: T) -> None:
        self._value: T | None = value

    def take(self) -> T:
        value, self._value = self._value, None
        if value is None:  # explicit: `python -O` strips an assert
            raise RuntimeError("_Once.take called twice")
        return value


def _dem_mesh(
    held: _Once[DemTile],
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
) -> _DemMesh:
    """Subsample, triangulate, sample or refine, and trim ``held``'s tile.

    Without ``tolerance``, z is sampled bilinearly at the stride grid. With it,
    the stride grid is the start mesh, refined against the DEM's nodes; with
    ``domain`` (already in the DEM's CRS; ``domain_name`` is its file's) the
    start is the polygon's rings and ``features``' lines instead.
    ``min_angle`` > 0 improves the start's angles first, a point accepted only
    for a gain of ``gain`` degrees when ``gain`` >= 0; ``feet`` inserts a
    node's foot on a nearby constraint segment. With ``grid`` and its
    ``checks``, the refined mesh is checked against the source's nodes after
    the tile is dropped. ``lines`` and ``steep`` build the line and slope
    tolerances on the tile's grid, passed to every refinement call.
    Returns the mesh, the record entries it knows, and the ``--stats`` sizes;
    ``clock`` gets the phases.
    ``dem`` names the source in messages. Every refusal is a usage error in the
    engine's own words, and no file is written.
    """
    tile = held.take()
    meta = tile.meta
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
        run = _engine(started.vertices, chains, delaunay, spacing, clock)
        if features is not None:
            values["start_mesh"] = "the domain outline and the feature lines"
            feature_counts = (len(chains) - len(rings), len(started.vertices) - domain_vertices)
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
        field = None
        if lines is not None:
            field, line_values = _line_field(lines, tile, tolerance)
            values |= line_values
        slope = None
        if steep is not None:
            with clock.phase("slope"):
                ramp = (steep.near_m, tolerance, steep.start_deg, steep.end_deg)
                slope = SlopeTolerance(to_core(tile), *ramp, threads=threads)
        t0 = time.perf_counter_ns()
        out = refine(
            to_core(tile),
            run.mesh,
            edges,
            masks,
            tolerance=tolerance,
            min_angle_deg=min_angle,
            constraint_feet=feet,
            min_gain_deg=gain,
            field=field,
            slope=slope,
            threads=threads,
        )
        _refine_phases(clock, (time.perf_counter_ns() - t0) / 1e9, out)
        if not out.ok():
            raise typer.BadParameter(f"{dem}: {out.message}", param_hint="--dem")
        strip = edge_strip.generate(to_core(tile), out, clock)  # while the tile is held
        if grid is None or checks is None:
            final = edge_strip.run(
                to_core(tile), strip, out, tolerance, clock, feet, field, slope, threads
            )
            del tile
            # Refine's maximum and the strip run's make an upper bound.
            max_error = max(out.max_error, final.max_error)
            near_error = max(out.max_error_near, final.max_error_near)
            share = max(out.max_slope_share, final.max_slope_share)
            inside, tight = out.slope_nodes
            values["line_check_dem_nodes_inserted"] = final.nodes_inserted
        else:
            del tile  # the final check runs without the target tile
            final, n = final_check.run(
                out, grid, checks, tolerance, clock, strip, feet, field, slope, threads
            )
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
            raise typer.BadParameter(f"{dem}: {final.message}", param_hint="--dem")
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
        }
        with clock.phase("trim"):
            trimmed = trim(
                vertices=final.vertices,
                triangles=final.triangles,
                edges=final.edges,
                edge_masks=final.masks,
                z=final.z,
                valid=final.valid,
            )
        values |= {
            "start_min_angle_deg": min_angle,
            "start_quality_gain_deg": gain,
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
            "start_vertices_between_dem_nodes": _off_node(np.asarray(run.mesh.vertices), meta),
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
        Path, typer.Option("--out", help="Where to write the polygon: .geojson or .json.")
    ],
    dem: Annotated[list[str], typer.Option("--dem", help=f"{DEM_HELP}.")],
    seed: Annotated[
        tuple[float, float],
        typer.Option("--seed", metavar="X Y", help="A point in the lake, in --seed-crs."),
    ],
    seed_crs: Annotated[
        str, typer.Option("--seed-crs", help="The seed's CRS, anything pyproj reads (LON LAT).")
    ] = "EPSG:4326",
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
    out_parent: Annotated[
        Path | None,
        typer.Option("--out-parent", help="Refuse any output path resolving outside this."),
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
    target = _destination(out, out_parent, out.stem)
    try:
        found = None if lakes is None else read_lake_polygons(lakes, lakes_layer, seed, seed_crs)
    except FeatureError as exc:
        raise typer.BadParameter(str(exc), param_hint="--lakes") from exc
    res = _resources(threads, memory_budget)
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
        fetch = "" if exc.box is None else " ".join(repr(v).removesuffix(".0") for v in exc.box)
        given = f" --cache {cache}" if cache else ""
        hint = f"; run: rasputin fetch {dem[0]} --bbox {fetch}{given}" if fetch else ""
        raise typer.BadParameter(f"{exc}{hint}", param_hint="--dem") from exc
    except OSError as exc:
        raise typer.BadParameter(f"cannot read {exc.filename}: {exc}", param_hint="--dem") from exc
    except ValueError as exc:
        hint = "--lakes" if isinstance(exc, LakeError) else "--dem"
        raise typer.BadParameter(_words(exc), param_hint=hint) from exc
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
    out_dir: Annotated[Path, typer.Option("--out-dir", help="Where to write the files.")],
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
    out_parent: Annotated[
        Path | None,
        typer.Option("--out-parent", help="Refuse any output path resolving outside this."),
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
    target = _destination(out_dir, out_parent, out_dir.name, "--out-dir")
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
    res = _resources(threads, memory_budget)
    repository = _dem_repository(dem, cache, res, z_unit)
    try:
        geographic = repository.footprints()[0].meta.geographic
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(
            _words(exc) if isinstance(exc, ValueError) else str(exc), param_hint="--dem"
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
            typer.echo(f"Error: {_words(exc)}", err=True)
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
    (CATCHMENTS, (catchment, station_catchments)),
    (DATA, (fetch, fetch_stations, sources)),
    (OTHER, (version,)),
):
    for command in commands:
        app.command(rich_help_panel=panel)(command)
app.command(rich_help_panel=OTHER, help=MCP_HELP)(mcp)
app.command(hidden=True)(draw)
