"""Command-line interface for the rasputin terrain engine.

This module is the **single composition root**: joining the file system to the
engine is its whole job. Paths from the command line also reach
``tin_engine.domain``, ``tin_engine.dem_input`` and ``io/repository.py``, which
read the files they name (``tin_engine.raster`` also imports ``_core``, to
build the one core raster, per ``project_structure.md``): ``viz/`` is written
against protocols and never
names a core type, while the core never sees a file, a path or a CRS. Everything
that has to know both sides lives here. (``project_structure.md``'s rule that
exactly one module constructs a core *raster* is about ``tin_engine.raster``.)

The `draw` command in particular:

1. looks a fixture up in ``viz.fixtures.GALLERY``;
2. maps the fixture's string roles onto ``ChainRole`` and calls ``build_pslg``,
   then ``node``, then ``triangulate``;
3. hands the **most-processed graph that exists** to ``build_scene`` as the
   ``PslgLike``, together with the vocabulary that graph's roles speak;
4. passes the engine's own words to ``render_svg`` as text;
5. resolves and validates the output path, and writes the bytes.

Step 3 is not a convenience. A fixture the validator *rejects* has no ``Pslg``
at all -- ``degenerate`` is refused with ``PslgError.DegenerateRing`` before
``node`` is ever called -- so the only ``PslgLike`` that exists for such a run is
the fixture, and drawing it is what keeps the failure presentation from being a
blank page. Step 3 is also why ``viz/`` names no role type: the mapping in
:data:`ROLES` is this module's to own, and so is the vocabulary that goes back
out with the source.

And the two travel as ONE value, :class:`Attempt`, because they are not
independent. A fixture's roles are strings; a ``NodedPslg``'s are ``ChainRole``
members, and ``ChainRole.Outer == "outer"`` is False. Swap the source without
swapping the vocabulary and no ring closes: one ``MASKED_EDGE_WITHOUT_CHAIN``
per ring, drawn in the alarm colour over a mesh the engine was perfectly happy
with, which every status assertion in the suite would pass. A source without its
roles is therefore made unrepresentable rather than merely discouraged.
"""

from __future__ import annotations

import asyncio
import csv
import json
import math
import os
import shlex
import sys
import tempfile
import time
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Annotated, TextIO

import numpy as np
import numpy.typing as npt
import typer
from pydantic import ValidationError

from tin_engine import edge_strip, final_check, installed_version
from tin_engine._core import (
    ChainRole,
    IndexedMesh2,
    NodedPslg,
    RefineOutcome,
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
from tin_engine.catchment_batch import BatchRequest, NoRiverLine, StationResult, run_batch, seed_for
from tin_engine.chains import start_chains
from tin_engine.crs import crs_label, reprojector, same_crs, transform_description, transform_label
from tin_engine.dem_input import (
    CachedSource,
    DemInput,
    DemRequest,
    OutCrsRequiredError,
    open_dem,
    repository_for,
)
from tin_engine.domain import DomainError, DomainPolygon, read_domain
from tin_engine.elevation import Trimmed, trim
from tin_engine.feature_input import (
    CLASS_MAPS,
    ClassMap,
    FeatureError,
    FeatureRequest,
    FeatureSet,
    FeatureSource,
    open_features,
    read_lakes,
)
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.gauge import Gauge, Placement
from tin_engine.grid_domain import default_stride, refine_start_stride, subsample
from tin_engine.hydrography import RiverSegment, Station
from tin_engine.io.cog import NotCached
from tin_engine.io.geojson import catchment_geojson
from tin_engine.io.models import DemTile, RasterMeta
from tin_engine.io.ply import write_ply
from tin_engine.io.repository import DemRepository
from tin_engine.io.rivers import read_segments
from tin_engine.io.station_set import read_references, read_stations
from tin_engine.io.vtk_legacy import write_vtk
from tin_engine.landcover import label_triangles
from tin_engine.mosaic import Bounds, Seam
from tin_engine.palettes import PALETTES, paraview_preset
from tin_engine.raster import to_core
from tin_engine.run_record import (
    SEAMS_AGREE,
    RunRecord,
    Value,
    as_json,
    file_fields,
    flat_record,
    ordinal,
    plural,
    refined_record,
    stride_record,
    summary,
)
from tin_engine.sources import SOURCES, STATION_SOURCES
from tin_engine.stats import PhaseClock, Report, Sizes, quality, render
from tin_engine.target_grid import Block, TargetGrid
from tin_engine.viz.fixtures import GALLERY, Fixture
from tin_engine.viz.protocols import PslgLike
from tin_engine.viz.scene import build_scene
from tin_engine.viz.style import PropertyStroke, SvgStyle
from tin_engine.viz.svg import render_svg

app = typer.Typer(
    name="rasputin",
    help="Parallel TIN engine for terrain meshing.",
    no_args_is_help=True,
)

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
#: Grown once, at increment 8: `bridge-over-lake`'s shoreline is an area feature
#: the vocabulary can name, and without an entry here its token is never emitted
#: -- the shoreline would draw as a plain breakline and the bridge picture would
#: lose the thing it exists to show. `coastline` sits between the two tokens
#: above, water still over infrastructure, and `svg.py` gained its CSS rule in
#: the same commit: an entry with no rule draws identically to the row above it.
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
# Increment 20, C2 (a): the start mesh's minimum angle, degrees.
DEFAULT_START_MIN_ANGLE = 25.0


@app.callback()
def main() -> None:
    """Parallel TIN engine for terrain meshing.

    The callback exists to keep Typer in sub-command mode. With a single
    registered command and no callback, Typer collapses that command into the
    root, and `rasputin version` becomes a usage error.
    """


def bounds_checks() -> str:
    """``off``, or ``on (<mode>)`` naming the checks the loaded ``_core`` has."""
    return "off" if hardening == "none" else f"on ({hardening})"


@app.command()
def version() -> None:
    """Print the installed rasputin version, then whether it is bounds-checked."""
    typer.echo(installed_version())
    typer.echo(f"bounds checks: {bounds_checks()}")


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
    field rather than a module constant because it is a property of the source.
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


def _triangulated(
    fixture: Fixture, delaunay: bool, spacing: float, clock: PhaseClock | None = None
) -> Attempt:
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
    run = _engine(fixture.vertices, chains, delaunay, spacing, clock)
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


def _destination(out: Path | None, out_parent: Path | None, name: str) -> Path:
    """Where to write, having refused every path this process should not touch.

    A hostile boundary, per the python skill. ``resolve()`` first, so that
    ``..`` cannot smuggle a write out of an explicitly permitted parent; refuse
    a symlinked ``--out`` rather than following it, because a link pointing back
    *inside* the permitted parent passes containment and still overwrites a file
    the user never named; and turn a missing directory into a message rather
    than a traceback.

    The symlink check is on the **final component only**: a symlinked parent
    directory is resolved by ``resolve()`` and then judged by containment, which
    is the right answer for ``--out-parent`` and is all this boundary claims.
    """
    if out is None:
        # "The common case is show me this now": a fresh directory of its own,
        # so a second run cannot overwrite the picture still open in a browser.
        return Path(tempfile.mkdtemp(prefix="rasputin-")) / f"{name}.svg"
    if out.is_symlink():
        raise typer.BadParameter(f"{out} is a symlink; refusing to follow it")
    target = out.resolve()
    if out_parent is not None:
        permitted = out_parent.resolve()
        if not target.is_relative_to(permitted):
            raise typer.BadParameter(f"{target} is outside the permitted parent {permitted}")
    if not target.parent.is_dir():
        raise typer.BadParameter(f"{target.parent} is not an existing directory")
    return target


@app.command()
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

    Exit code 0 means a picture was produced, which includes the two deliberate
    failure fixtures -- ``degenerate``, which the PSLG validator refuses before
    the noder is reached, and ``hole-in-hole``, which the backend refuses -- and
    a noding this spacing cannot resolve: a drawn failure is still a drawn
    picture. A caller who wants the engine's verdict rather than the renderer's
    therefore reads the header band, not the exit code; the status and the
    engine's own words are written there.

    The noder's round cap is deliberately not an option. It is a bound on a loop
    rather than a parameter of the answer -- at any value the outcome is the same
    mesh or ``NotConverged``, never a different mesh -- and on the corner graze,
    whose orbit has period two, ``Ok`` is unreachable at every setting. A knob
    whose visible effect on the failure a user is most likely to meet is "the
    same refusal, slower" is worse than not having one. The spacing is the lever,
    and it is the one the band names.
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


#: What the suffix of ``--out`` selects (increment 13, ruling 8).
MESH_SUFFIXES = (".vtk", ".ply")


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
    draws nothing (ruling 6).

    A constrained mesh edge with no entry in the noded graph gets 0, which
    ``_core.pyi`` defines as *unclassified* rather than *wrong*.
    """
    masks = _chain_masks(pslg)
    pairs = sorted(_masked_pairs(mesh))
    return (
        np.array(pairs, dtype=np.uint32).reshape(-1, 2),
        np.array([masks.get(pair, 0) for pair in pairs], dtype=np.uint32),
    )


@app.command()
def mesh(
    out: Annotated[
        Path, typer.Option("--out", help="Where to write: .vtk for ParaView, .ply for QGIS.")
    ],
    name: Annotated[
        str | None, typer.Argument(help="Gallery fixture to write; or give --dem instead.")
    ] = None,
    dem: Annotated[
        list[str] | None,
        typer.Option(
            "--dem",
            help="A GeoTIFF DEM, several (repeat --dem), one directory of tiles, or a "
            f"cached catalogue source ({', '.join(SOURCES)}): mesh its extent with z "
            "sampled from it. A file named like a source is written as a path: ./glo30.",
        ),
    ] = None,
    cache: Annotated[
        Path | None,
        typer.Option(
            "--cache",
            help="The tile cache a catalogue --dem is read from. Default: $RASPUTIN_DATA/cache.",
        ),
    ] = None,
    bbox: Annotated[
        tuple[float, float, float, float] | None,
        typer.Option(
            "--bbox",
            metavar="BOX",
            help="With --dem, mesh only XMIN YMIN XMAX YMAX in the DEM's CRS, "
            "snapped outward to its nodes.",
        ),
    ] = None,
    stride: Annotated[
        int | None,
        typer.Option("--stride", help="With --dem, every Nth DEM node. Default: <= 256 a side."),
    ] = None,
    flat: Annotated[
        bool, typer.Option("--flat", help="There is no elevation source; write z = 0.")
    ] = False,
    out_edges: Annotated[
        Path | None,
        typer.Option("--out-edges", help="With .ply, also write the constraint edges as a PLY."),
    ] = None,
    out_parent: Annotated[
        Path | None,
        typer.Option("--out-parent", help="Refuse any output path resolving outside this."),
    ] = None,
    crs: Annotated[
        str, typer.Option("--crs", help="Free text recorded as a header comment. Not validated.")
    ] = "",
    binary: Annotated[
        bool, typer.Option("--binary/--ascii", help="Packed records, or text `head` can read.")
    ] = False,
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
            help="With --dem, refine until every triangle is within this many metres "
            "of the DEM; --stride then sets the start grid. Default: no refinement.",
        ),
    ] = None,
    domain: Annotated[
        Path | None,
        typer.Option(
            "--domain",
            help="With --dem and --tolerance, mesh only this polygon (.geojson, .json or "
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
            help="With --dem, mesh in this CRS: a DEM in another CRS (required for a "
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
            "0 is off, at most 35. Default: 25.",
        ),
    ] = None,
    no_constraint_feet: Annotated[
        bool,
        typer.Option(
            "--no-constraint-feet",
            help="With --tolerance, insert a DEM node close to a constraint segment as "
            "before, instead of its foot on the segment.",
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
            help=f"The Nth --features' map from attributes to edge bits: {', '.join(CLASS_MAPS)}. "
            "Default: property.",
        ),
    ] = None,
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
    """Write a mesh as legacy VTK or as PLY, by the suffix of ``--out``.

    The mesh is a gallery fixture's, or, with ``--dem PATH``, a GeoTIFF's
    (increment 12): every ``--stride``-th DEM node inside the grid's outer ring,
    triangulated, with z sampled bilinearly from the DEM. With ``--tolerance``
    (increment 14) that grid is only the start: it is refined at DEM nodes until
    every triangle's max error is within the tolerance, and z is read at the
    nodes. With ``--domain`` (increment 16) the start is the polygon's rings
    alone, their vertices where the file puts them with bilinear z, and only
    the inside is meshed. Vertices where the DEM has no data are dropped with
    their triangles, and the count is reported.
    The file records the DEM's CRS, so ``--crs`` and ``--flat`` are refused.
    ``--dem DIR`` or several ``--dem`` files are stitched into one grid first,
    cut to ``--bbox`` if given (increment 15a, ``tin_engine.mosaic``).

    ``.vtk`` is one file for ParaView: triangles, constraint lines, their
    feature masks, one 0/1 array per feature that occurs, and the vocabulary
    (``13-bundled-mesh.md``). There is no ``--format``, because it could
    contradict the suffix. Text is the default for both formats.

    ``.ply`` is for QGIS. Two files, never one holding both element types:
    MDAL's own caveat is that a host application expects either a 1D mesh or a
    2D one, so a file carrying faces AND edges can load as nothing at all, in
    silence. ``--out`` gets the surface; ``--out-edges`` gets the constraints,
    or they are not written.
    Both repeat the identical vertex block, which is what makes the two layers
    register on each other when a person loads them side by side.

    Unlike ``draw``, a failed engine run is a non-zero exit and no file. A
    picture of a failure is still a picture and worth producing; there is no
    such thing as a picture of a failed file, so the refusal is reported in the
    engine's own words instead.

    ``--stats`` (increment 17) adds a report and changes nothing else: the
    clock always runs, the quality pass and the report only with the flag.

    A ``--features-map`` with class codes (``corine``, ``corine-water``,
    ``clc18_kode``) also gives every triangle its polygon's code, in the cell
    array ``land_cover_code`` of a ``.vtk`` or the face property of a ``.ply``
    (increment 16c); ``rasputin palette corine`` writes natural colours for it.

    What the run did is one record (increment 25, ``tin_engine.run_record``):
    the file carries the fields a user of the mesh needs, ``--stats`` and
    ``--record`` all of it, and stderr one summary.
    """
    clock = PhaseClock()
    dem_run: _DemMesh | None = None
    codes: npt.NDArray[np.int32] | None = None
    codes_text = ""
    seams: tuple[Seam, ...] = ()
    if (name is None) == (not dem):
        raise typer.BadParameter(
            "give a gallery fixture name or --dem PATH, exactly one of the two",
            param_hint="--dem",
        )
    feature_paths = features or []
    feature_crss = features_crs or []
    feature_layers = features_layer or []
    feature_maps = features_map or []
    # R3: one length rule per paired option. A list longer than --features is
    # refused naming the offending flag and the counts (M = 0 subsumes the old
    # "applies only with --features"); a shorter list takes per-source defaults.
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
    if not dem and (domain is not None or domain_crs is not None):
        raise typer.BadParameter("applies only with --dem", param_hint="--domain")
    if not dem and bbox is not None:
        raise typer.BadParameter("applies only with --dem", param_hint="--bbox")
    if domain is None and domain_crs is not None:
        raise typer.BadParameter("applies only with --domain", param_hint="--domain-crs")
    if domain is not None and bbox is not None:
        raise typer.BadParameter("--bbox and --domain exclude each other", param_hint="--bbox")
    if name is not None and name not in GALLERY:
        raise typer.BadParameter(f"unknown fixture {name}; the gallery is: {', '.join(GALLERY)}")
    if out.suffix not in MESH_SUFFIXES:
        raise typer.BadParameter(
            f"{out} has no known suffix; use {' or '.join(MESH_SUFFIXES)}", param_hint="--out"
        )
    if out.suffix == ".vtk" and out_edges is not None:
        raise typer.BadParameter(
            "a .vtk file already carries the constraint edges", param_hint="--out-edges"
        )

    if dem:
        if flat:
            raise typer.BadParameter("--dem samples z from the DEM", param_hint="--flat")
        if crs:
            raise typer.BadParameter("--dem records the DEM's own CRS", param_hint="--crs")
        if stride is not None and stride < 1:
            raise typer.BadParameter(
                f"must be a positive integer, got {stride}", param_hint="--stride"
            )
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
        if no_constraint_feet and tolerance is None:
            raise typer.BadParameter(
                "--no-constraint-feet needs --tolerance", param_hint="--no-constraint-feet"
            )
        paths, cached = _dem_sources(dem, cache_root(cache, os.environ))
        given = None
        if domain is not None:
            with clock.phase("domain read"):
                try:
                    given = read_domain(domain, domain_crs)
                except DomainError as exc:
                    raise typer.BadParameter(str(exc), param_hint="--domain") from exc
        try:
            opened = _open_dem(paths, bbox, clock, given, cached, out_crs)
        except NotCached as exc:
            area = f"--domain {domain}" if domain is not None else ""
            if bbox is not None:
                area = "--bbox " + " ".join(repr(v).removesuffix(".0") for v in bbox)
            frame = f"--out-crs {shlex.quote(out_crs)}" if out_crs else ""
            flags = " ".join(filter(None, (area, frame, f"--cache {cache}" if cache else "")))
            run = f"rasputin fetch {dem[0]} {flags}".rstrip()
            raise typer.BadParameter(f"{exc}; run: {run}", param_hint="--dem") from exc
        # 15e fix 3: what is read later comes out of `opened`, which then goes,
        # so the tile's one reference is the holder `_dem_mesh` empties.
        label, plan, grid, source_crs = opened.label, opened.plan, opened.grid, opened.source_crs
        dem_domain, checks, opened_seams = opened.domain, opened.checks, opened.seams
        held, dem_crs = _Once(opened.tile), opened.tile.meta.crs
        del opened
        found = None
        sources: tuple[FeatureSource, ...] = ()
        if feature_paths and dem_domain is not None:
            sources = _feature_sources(feature_paths, feature_crss, feature_layers, feature_maps)
            found = _open_features(sources, dem_domain, dem_crs, clock)
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
            DEFAULT_START_MIN_ANGLE if start_min_angle is None else start_min_angle,
            not no_constraint_feet,
            found,
            grid,
            checks,
        )
        surface_mesh, meta, values = dem_run.trimmed, dem_run.meta, dict(dem_run.values)
        names = [t.name for t in plan.tiles]
        size = f"{meta.cols} columns x {meta.rows} rows"
        if len(names) > 1:
            typer.echo(f"DEM: {len(names)} files, {size}", err=True)
        assumed = " (assumed: the DEM file does not say)" if meta.vertical_unit_assumed else ""
        source = cached.source if cached is not None else _ascii("; ".join(names))
        values |= {"crs": dem_crs, "dem_source": source, "dem_vertical_unit": "metres" + assumed}
        if grid is None:
            dx, dy = meta.delta_x, meta.delta_y
            apart = f"{dx:g} m" if dx == dy else f"{dx:g} m x {dy:g} m"
            values["dem_grid"] = f"{size}, {apart} apart"
        else:  # 15c-2, D7
            how = transform_description(source_crs, dem_crs)
            values |= {"dem_crs": source_crs, "dem_transform": how}
            values["resampled_grid"] = f"{grid.spacing} m square grid in {dem_crs}, {size}"
        if cached is not None:  # B16 (a): the notes the source asks to travel with it
            remote = SOURCES[cached.source]
            values["dem_credit"] = _ascii(remote.credit)
            values["licence_note"] = _ascii(remote.licence_note)
            values["cite"] = _ascii("; ".join(remote.cite)) or None
        if cached is not None or len(paths) > 1 or paths[0].is_dir():  # R11: the files used
            seams = opened_seams
            values["dem_tiles"] = _ascii("; ".join(names))
            values["dem_seams"] = _ascii("; ".join(s.entry() for s in seams)) or SEAMS_AGREE
        if given is not None:
            how = transform_label(given.crs, dem_crs)
            values |= {"domain_crs": crs_label(given.crs), "domain_transform": how}
        if found is not None and sources and dem_run.feature_counts is not None:
            # R6: one entry per source, joined. Lines and their vertices are
            # whole-run (D5, a property of the merged PSLG, not a source), so
            # every entry carries the run totals; the per-source part is the
            # name, layer, map and feature count (found.counts[i]).
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
            values |= {"features": _ascii("; ".join(texts)), "features_crs": "; ".join(crs_texts)}
            values["features_transform"] = "; ".join(transforms)
            values["features_notice"] = "; ".join(notices) or None
            # R5: labelling runs iff a source carries codes (D1's single-system
            # invariant is enforced in _feature_sources). Every coded source's
            # polygons are labelled together over the merged FeatureSet.
            coded = next((s.class_map for s in sources if s.class_map.codes), None)
            if coded is not None:
                codes, codes_text = _land_cover(dem_run.trimmed, found, coded, snap_spacing, clock)
                values["land_cover_codes"] = codes_text
        build = stride_record if tolerance is None else refined_record
        record = build(triangles=len(surface_mesh.triangles), **values)
    else:
        assert name is not None
        if stride is not None:
            raise typer.BadParameter("applies only with --dem", param_hint="--stride")
        if tolerance is not None:
            raise typer.BadParameter("applies only with --dem", param_hint="--tolerance")
        if start_min_angle is not None:
            raise typer.BadParameter("applies only with --dem", param_hint="--start-min-angle")
        if no_constraint_feet:
            raise typer.BadParameter("applies only with --dem", param_hint="--no-constraint-feet")
        if not flat:
            raise typer.BadParameter(
                "a gallery fixture has no elevation source, so z has none; pass --flat "
                "to write z = 0 and say so in the file, or mesh a DEM with --dem"
            )
        label = name
        surface_mesh = _fixture_mesh(name, delaunay, snap_spacing, clock)
        record = flat_record(triangles=len(surface_mesh.triangles), crs=crs)
        # --crs is unvalidated free text by ruling 5, so the writer's refusals
        # are refusals a person meets by typing, not internal invariants. Turn
        # the writer's ValueError into the usage error it is, in the one place
        # that knows the text came from the command line.
        try:
            write_ply(np.zeros((1, 3)), faces=np.zeros((0, 3)), comments=_comments(record))
        except ValueError as exc:
            raise typer.BadParameter(str(exc), param_hint="--crs") from exc

    vertices = surface_mesh.vertices
    surface = _destination(out, out_parent, label)
    targets = [surface]
    if out_edges is not None:
        # Resolved, because two spellings of one path are still one file. Both
        # writes succeed, the second overwrites the first, the command echoes
        # two paths and exits 0 -- the caller has lost the surface they asked
        # for and nothing said so.
        constraints = _destination(out_edges, out_parent, label)
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

    fields, comments = file_fields(record), _comments(record)
    encoders: list[Callable[[], bytes]]
    if out.suffix == ".vtk":
        encoders = [
            lambda: write_vtk(
                vertices,
                triangles=surface_mesh.triangles,
                edges=surface_mesh.edges,
                edge_masks=surface_mesh.edge_masks,
                vocabulary=DEFAULT_VOCABULARY,
                fields=fields,
                binary=binary,
                triangle_codes=codes,
                land_cover_codes=codes_text,
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
            ),
            lambda: write_ply(
                vertices,
                edges=surface_mesh.edges,
                edge_properties=surface_mesh.edge_masks,
                ascii=not binary,
                comments=comments,
                vocabulary=DEFAULT_VOCABULARY,
            ),
        ][: len(targets)]
    for target, encode in zip(targets, encoders, strict=True):
        with clock.phase("write: encode"):
            data = encode()
        with clock.phase("write: disk"):
            target.write_bytes(data)
        typer.echo(f"{target}")
    if stats is not None:
        _write_report(clock, report_target, surface_mesh, dem_run, targets, seams, record)
    if record_target is not None:  # D5: after the mesh and the report
        record_target.write_text(as_json(record, installed_version(), _command()), "ascii")
        typer.echo(f"{record_target}")


def _comments(record: RunRecord) -> list[str]:
    """The record's file fields as ``.ply`` header comments, ``name value``."""
    return [f"{name} {value}" for name, value in file_fields(record)]


def _land_cover(
    trimmed: Trimmed, found: FeatureSet, cmap: ClassMap, spacing: float, clock: PhaseClock
) -> tuple[npt.NDArray[np.int32], str]:
    """16c, R1-R3: a code per triangle from the coded polygons, timed as the
    phase ``land cover``, with its stderr line; the codes and their text."""
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


@app.command()
def palette(
    name: Annotated[str, typer.Argument(help=f"The palette: {', '.join(PALETTES)}.")],
    out: Annotated[
        Path | None, typer.Option("--out", help="Write the preset here, not to stdout.")
    ] = None,
) -> None:
    """Write a ParaView colour preset for ``land_cover_code`` (increment 16c).

    Import it in ParaView's Colour Map Editor (Choose Preset, Import), then
    colour by ``land_cover_code``; tick Interpret Values As Categories if the
    preset does not.
    """
    if name not in PALETTES:
        raise typer.BadParameter(f"unknown palette {name}; use {', '.join(PALETTES)}")
    table, title = PALETTES[name]
    text = json.dumps(paraview_preset(table, title), indent=1) + "\n"
    if out is None:
        typer.echo(text, nl=False)
        return
    out.write_text(text)
    typer.echo(f"{out}")


def _report_target(
    stats: str | None, out_parent: Path | None, label: str, meshes: list[Path]
) -> Path | None:
    """R2: where ``--stats`` writes, None for stdout (``-``) or no report,
    checked before any mesh file is written."""
    if stats is None or stats == "-":
        return None
    target = _destination(Path(stats), out_parent, label)
    if target in meshes:
        raise typer.BadParameter(
            f"resolves to {target}; the report would overwrite the mesh", param_hint="--stats"
        )
    return target


def _record_target(
    record: str | None, out_parent: Path | None, label: str, taken: list[Path | None]
) -> Path | None:
    """D5: where ``--record`` writes, checked before any file is written; never
    standard output (``--stats -``'s) nor a file this run also writes."""
    if record is None:
        return None
    if record == "-":
        raise typer.BadParameter("standard output is --stats -'s", param_hint="--record")
    target = _destination(Path(record), out_parent, label)
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
    dem_run: _DemMesh | None,
    files: list[Path],
    seams: tuple[Seam, ...],
    record: RunRecord,
) -> None:
    """Build the report from what ran and write it, or print it for ``-``.
    The total stops here; the quality pass is timed on its own line (R4)."""
    total = clock.elapsed()
    t0 = time.perf_counter_ns()
    meta = dem_run.meta if dem_run else None
    sizes = Sizes(
        output_vertices=len(trimmed.vertices),
        output_triangles=len(trimmed.triangles),
        constraint_edges=len(trimmed.edges),
        files=[(f.name, f.stat().st_size) for f in files],
        dem_nodes=(meta.rows, meta.cols) if meta else None,
        dem_spacing=(meta.delta_x, meta.delta_y) if meta else None,
        domain_vertices=dem_run.domain_vertices if dem_run else None,
        domain_holes=dem_run.domain_holes if dem_run else None,
        start_vertices=dem_run.start_vertices if dem_run else None,
        start_triangles=dem_run.start_triangles if dem_run else None,
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
            threads=os.cpu_count() if sizes.start_vertices is not None else None,
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


def _fixture_mesh(name: str, delaunay: bool, spacing: float, clock: PhaseClock) -> Trimmed:
    """A gallery fixture's mesh at z = 0, with its constraint edges."""
    attempt = _triangulated(GALLERY[name], delaunay=delaunay, spacing=spacing, clock=clock)
    if attempt.mesh is None or not isinstance(attempt.source, NodedPslg):
        raise typer.BadParameter(
            f"{name} has no mesh to write: {attempt.status}. {attempt.message}"
        )
    flat_vertices = np.asarray(attempt.mesh.vertices)
    with clock.phase("start mesh: constraint edges"):
        edges, masks = _constraint_arrays(attempt.mesh, attempt.source)
    return Trimmed(
        vertices=np.column_stack([flat_vertices, np.zeros(len(flat_vertices))]),
        triangles=np.asarray(attempt.mesh.triangles, dtype=np.uint32),
        edges=edges,
        edge_masks=masks,
        dropped=0,
    )


@app.command()
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
        Path | None, typer.Option("--cache", help="The tile cache. Default: $RASPUTIN_DATA/cache.")
    ] = None,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Read the headers, print the plan, write nothing.")
    ] = False,
    refresh: Annotated[
        bool, typer.Option("--refresh", help="Discard the source's cache and fetch anew.")
    ] = False,
    connections: Annotated[int, typer.Option("--connections", min=1)] = 8,
) -> None:
    """Copy what a mesh of the domain reads from a remote source into the tile
    cache (increment 23a-2), so ``rasputin mesh --dem SOURCE`` runs offline."""
    import asyncio

    from tin_engine.fetch.http import FetchError, RangeClient
    from tin_engine.fetch.plan import FetchRequest
    from tin_engine.fetch.run import fetch as run
    from tin_engine.io.repository import CacheError, CacheWriter

    if source not in SOURCES:
        raise typer.BadParameter(f"{source} is not a catalogue source ({', '.join(SOURCES)})")
    if (domain is None) == (bbox is None):
        raise typer.BadParameter("give exactly one of --domain and --bbox", param_hint="--bbox")
    root = cache_root(cache, os.environ)
    if root is None:
        raise typer.BadParameter(
            "no cache: set RASPUTIN_DATA (the data root; the cache is $RASPUTIN_DATA/cache) "
            "or pass --cache",
            param_hint="--cache",
        )
    try:
        given = read_domain(domain, domain_crs) if domain is not None else None
        corners = ("x_min", "y_min", "x_max", "y_max")
        box = None if bbox is None else Bounds(**dict(zip(corners, bbox, strict=True)))
        request = FetchRequest(
            source=source, domain=given, box=box, out_crs=out_crs, connections=connections,
            dry_run=dry_run, refresh=refresh,
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
        report = asyncio.run(run(request, SOURCES[source], client, writer, progress=progress))
    except (FetchError, CacheError, ValueError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(1) from exc
    for plan in report.plans:
        typer.echo(f"{plan.object_id}: {len(plan.blocks):,} blocks, {plan.bytes:,} bytes to fetch")
    head = "would fetch" if dry_run else "fetched"
    typer.echo(
        f"{source}: {report.objects} objects, {report.needed:,} blocks needed, "
        f"{report.present:,} present, {head} {report.fetched:,} ({report.empty} empty) in "
        f"{report.requests:,} requests, {report.bytes:,} bytes, {report.seconds:.1f} s"
    )
    if report.no_tile:
        typer.echo(f"no tile (sea): {', '.join(report.no_tile)}")


@app.command()
def fetch_stations(
    source: Annotated[str, typer.Argument(help=f"A station list: {', '.join(STATION_SOURCES)}.")],
    out_dir: Annotated[Path, typer.Option("--out-dir", help="Where to write the files.")],
    refresh: Annotated[
        bool, typer.Option("--refresh", help="Fetch again even if the files are there.")
    ] = False,
) -> None:
    """Fetch a list's stations, their catchment polygons and the river lines
    round them into --out-dir, with NOTICE.txt and a manifest (increment 29).
    Files already there are kept, and nothing is requested, unless --refresh."""
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


def cache_root(option: Path | None, environ: Mapping[str, str]) -> Path | None:
    """The tile cache (B7): ``--cache`` if given, else ``$RASPUTIN_DATA/cache``,
    else None. The one reader of the variable; ``mesh`` passes ``os.environ``."""
    if option is not None:
        return option
    data = environ.get("RASPUTIN_DATA")
    return Path(data) / "cache" if data else None


def _dem_sources(dem: list[str], root: Path | None) -> tuple[tuple[Path, ...], CachedSource | None]:
    """``--dem`` as given: one catalogue key alone is that cached source
    (23a-1); anything else is a path, so ``./glo30`` is a file."""
    keys = [d for d in dem if d in SOURCES]
    if not keys:
        return tuple(Path(d) for d in dem), None
    if len(dem) > 1:
        raise typer.BadParameter(
            f"{keys[0]} is a catalogue source; give it alone, not with other --dem",
            param_hint="--dem",
        )
    if root is None:
        raise typer.BadParameter(
            "no cache: set RASPUTIN_DATA (the data root; the cache is $RASPUTIN_DATA/cache) "
            "or pass --cache",
            param_hint="--dem",
        )
    return (), CachedSource(source=keys[0], cache=root)


def _open_dem(
    dem: tuple[Path, ...],
    bbox: tuple[float, float, float, float] | None,
    clock: PhaseClock,
    domain: DomainPolygon | None = None,
    cached: CachedSource | None = None,
    target_crs: str | None = None,
) -> DemInput:
    """``--dem`` and ``--bbox`` or the read ``--domain`` to one tile (increment
    15a and 15b, R11); every refusal, the reader's, the mosaic's or the
    domain's extent, is a usage error in its own words."""
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
                sources=dem, cached=cached, bounds=bounds, domain=domain, target_crs=target_crs
            )
            return open_dem(request)
    except NotCached:
        raise
    except OutCrsRequiredError as exc:  # Q11: one unwrapped line, outside the panel, to paste
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


def _feature_sources(
    paths: list[Path], crss: list[str], layers: list[str], maps: list[str]
) -> tuple[FeatureSource, ...]:
    """R1-R3: the repeatable ``--features`` flags as one source per path, paired
    by position; a shorter paired list defaults its missing entries. A flag that
    cannot apply is a usage error naming the 1-based index and path, because with
    several sources the flag alone no longer identifies which one."""
    built: list[FeatureSource] = []
    for i, path in enumerate(paths):
        name = maps[i] if i < len(maps) else "property"
        if name not in CLASS_MAPS:
            raise typer.BadParameter(
                f"unknown map {name} for --features {i + 1} ({path.name}); "
                f"use {', '.join(CLASS_MAPS)}",
                param_hint="--features-map",
            )
        layer = layers[i] if i < len(layers) else None
        if layer is not None and path.suffix.lower() != ".gpkg":
            raise typer.BadParameter(
                f"applies to a .gpkg only, not --features {i + 1} ({path.name})",
                param_hint="--features-layer",
            )
        crs = crss[i] if i < len(crss) else None
        built.append(FeatureSource(path=path, class_map=CLASS_MAPS[name], layer=layer, crs=crs))
    # R5/D1: one mesh has one land_cover_code system. Two sources with different
    # non-empty code systems are refused before any meshing. Every coded map
    # today is CORINE (same system), so this defensive refusal never fires now.
    systems = {s.class_map.codes for s in built if s.class_map.codes}
    if len(systems) > 1:
        raise typer.BadParameter(
            f"sources carry different code systems ({', '.join(sorted(systems))}); "
            "one mesh has one land_cover_code system",
            param_hint="--features-map",
        )
    return tuple(built)


def _open_features(
    sources: tuple[FeatureSource, ...], domain: DomainPolygon, dem_crs: str, clock: PhaseClock
) -> FeatureSet:
    """R5-R6 and R10: read and clip every source into one set, timed as
    ``features read`` and ``features clip``; a refusal is a usage error naming
    the feature."""
    t0 = time.perf_counter()
    try:
        found = open_features(FeatureRequest(sources=sources), domain, dem_crs)
    except FeatureError as exc:
        raise typer.BadParameter(str(exc), param_hint="--features") from exc
    clock.add("features read", time.perf_counter() - t0 - found.clip_seconds)
    clock.add("features clip", found.clip_seconds)
    typer.echo(
        f"features: {len(found.features)} kept ({found.clipped} cut at the domain outline), "
        f"{found.outside} outside the domain, {found.empty} empty",
        err=True,
    )
    for table in found.scanned:
        typer.echo(f"{table}: this layer has no spatial index, so every row was read", err=True)
    return found


def _ascii(text: str) -> str:
    """A file field's text with non-ASCII escaped, as `dem_tiles` records names."""
    return text.encode("ascii", "backslashreplace").decode("ascii")


def _words(exc: ValueError) -> str:
    """A refusal's own words: Pydantic's messages without its wrapper."""
    if isinstance(exc, ValidationError):
        return "; ".join(str(e["msg"]) for e in exc.errors())
    return str(exc)


@dataclass(frozen=True, slots=True)
class _DemMesh:
    """What ``_dem_mesh`` made: the mesh, then what ``--stats`` reports about
    the run (None where a row does not apply), and the record entries it knows
    (``run_record``'s builder arguments, increment 25)."""

    trimmed: Trimmed
    meta: RasterMeta
    domain_vertices: int | None
    domain_holes: int | None
    start_vertices: int | None
    start_triangles: int | None
    values: Mapping[str, Value]
    feature_counts: tuple[int, int] | None = None


class _Once[T]:
    """A value handed on once (15e, fix 3): `take` returns it and forgets it,
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
) -> _DemMesh:
    """Subsample, triangulate, sample or refine, and trim ``held``'s tile.

    Without ``tolerance`` this is increment 12's R6: z sampled bilinearly at
    the stride grid. With it, increment 14's R9: the stride grid is the start
    mesh, refined against the DEM's nodes; with ``domain`` (already in the
    DEM's CRS, 15b; ``domain_name`` is its file's), increment 16's R3, the
    polygon's rings are, and ``features``' lines (16b). ``min_angle`` > 0
    improves the start's angles first (increment 20); ``feet`` inserts
    constraint feet (increment 20b). With ``grid`` and its ``checks`` (15c-2),
    the refined mesh is checked against the source's nodes (D5), after the
    tile is dropped (15e, fix 3).
    Returns the mesh, the record entries it knows, and the ``--stats`` sizes;
    ``clock`` gets R5's phases.
    ``dem`` names the source in messages. Every refusal is a usage error in the
    engine's own words, and no file is written.
    """
    tile = held.take()
    meta = tile.meta
    values: dict[str, Value] = {}
    domain_vertices = domain_holes = None
    feature_counts = None
    if domain is not None:
        started = start_chains(domain, features.features if features else (), DEFAULT_VOCABULARY)
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
        t0 = time.perf_counter_ns()
        out = refine(
            to_core(tile),
            run.mesh,
            edges,
            masks,
            tolerance=tolerance,
            min_angle_deg=min_angle,
            constraint_feet=feet,
        )
        _refine_phases(clock, (time.perf_counter_ns() - t0) / 1e9, out)
        if not out.ok():
            raise typer.BadParameter(f"{dem}: {out.message}", param_hint="--dem")
        strip = edge_strip.generate(to_core(tile), out, clock)  # 15f, D6: while the tile is held
        if grid is None or checks is None:
            final = edge_strip.run(to_core(tile), strip, out, tolerance, clock)
            del tile
            # 15f, D7: refine's maximum and the strip run's make an upper bound.
            max_error = max(out.max_error, final.max_error)
            values["line_check_dem_nodes_inserted"] = final.nodes_inserted
        else:
            del tile  # 15e fix 3: phase 2 runs without the target tile
            final, n = final_check.run(out, grid, checks, tolerance, clock, strip=strip)
            max_error = final.max_error
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
            "snap_to_lines": "on" if feet else "off",
            "tolerance_m": tolerance,
            "max_error_m": max_error,
            "dem_nodes_outside_mesh": final.uncovered,
            "refinement_rounds": out.rounds,
            "points_inserted": out.inserted,
            "points_inserted_on_nodata": out.carved,
            "edge_flips": out.flips,
            "start_quality_points_inserted": out.quality_inserted,
            "start_quality_points_skipped": out.quality_skipped,
            "points_snapped_to_lines": out.feet,
            "snaps_refused": out.feet_refused,
            "start_vertices_between_dem_nodes": _off_node(np.asarray(run.mesh.vertices), meta),
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
    """R5 and R6: ``refine`` and its sub-rows; setup + output is the remainder."""
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
    col = np.round((xy[:, 0] - meta.x_min) / meta.delta_x)
    row = np.round((meta.y_max - xy[:, 1]) / meta.delta_y)
    node = (meta.x_min + col * meta.delta_x == xy[:, 0]) & (
        meta.y_max - row * meta.delta_y == xy[:, 1]
    )
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
    its reach in the file's CRS, which must be the DEM's."""
    segments, crs, _ = _segments(rivers)
    _reach_crs(crs, repository)
    ((x, y),) = reprojector(seed_crs, crs)([seed])
    try:
        placement, _, made = seed_for(Gauge(x=x, y=y), segments, crs, request)
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


@app.command()
def catchment(
    out: Annotated[
        Path, typer.Option("--out", help="Where to write the polygon: .geojson or .json.")
    ],
    dem: Annotated[
        list[Path],
        typer.Option("--dem", help="A GeoTIFF DEM, several (repeat --dem), or one directory."),
    ],
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
            help="Lake polygons (.gpkg, .geojson or .json); the one under the seed is the "
            "seed. Without it, the seed is the DEM node nearest the point.",
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
            help="River lines (GeoJSON, in the DEM's CRS): the seed is a gauge, placed on the "
            "nearest line and on the DEM's flow path along it (increment 29).",
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
) -> None:
    """Write the catchment of a lake, from the DEM, as a GeoJSON polygon in the
    DEM's CRS that ``mesh --domain`` reads (increment 22). The fine outline is
    drawn between DEM nodes, holes filled, then reduced to ``--outline-tolerance``
    keeping its area, with the seed inside. A catchment cut by the data's edge
    or by NoData is refused, and nothing is written. With ``--rivers`` the
    seed is a gauge on a mapped river (increment 29)."""
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
        found = None if lakes is None else read_lakes(lakes, lakes_layer, seed, seed_crs)
    except FeatureError as exc:
        raise typer.BadParameter(str(exc), param_hint="--lakes") from exc
    try:
        repository, _ = repository_for(tuple(dem))
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
        result = delineate(request, repository)
    except OSError as exc:
        raise typer.BadParameter(f"cannot read {exc.filename}: {exc}", param_hint="--dem") from exc
    except ValueError as exc:
        hint = "--lakes" if isinstance(exc, LakeError) else "--dem"
        raise typer.BadParameter(_words(exc), param_hint=hint) from exc
    for k, w in enumerate(result.windows, 1):
        b = w.bounds
        grown = f"window widened to the {', '.join(w.grown)}"
        typer.echo(
            f"window {k}: x {b.x_min:.0f}-{b.x_max:.0f}, y {b.y_min:.0f}-{b.y_max:.0f}, "
            f"{w.rows} x {w.cols} nodes, searched in {w.seconds:.2f} s, "
            + (grown if w.grown else "catchment inside the window"),
            err=True,
        )
    cell = result.meta.delta_x * result.meta.delta_y
    extra: dict[str, object] = {}
    if placement is not None and result.gauge is not None:
        extra = _placement_report(placement, result.gauge)
    elif result.lake_area is None:
        typer.echo(
            f"start: the outlet node at {result.seed} (an outlet must lie on the flow line; "
            "it is not moved there)",
            err=True,
        )
    else:
        typer.echo(
            f"start: lake of {result.lake_area / 1e6:.6f} km2, {result.seed_nodes} DEM nodes",
            err=True,
        )
    typer.echo(
        f"catchment: {result.nodes} DEM nodes, {result.nodes * cell / 1e6:.6f} km2",
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
    change = reduced.area - result.fine_area
    typer.echo(
        f"reduced outline: {kept} vertices, {reduced.area / 1e6:.6f} km2, difference "
        f"{change:.3g} m2 ({change / result.fine_area:.2g} relative), tolerance "
        f"{result.tolerance:g} m, {result.reduce_seconds:.2f} s",
        err=True,
    )
    properties = {"seed": list(seed), "seed_crs": seed_crs, **_outline_properties(result), **extra}
    target.write_bytes(catchment_geojson(reduced, result.crs, properties))
    typer.echo(f"{target}")


def _outline_properties(result: Catchment) -> dict[str, object]:
    """Increment 22's properties of a catchment file: the counts, the areas,
    the tolerance and the windows."""
    return {
        "nodes": result.nodes,
        "fine_vertices": len(result.fine.exterior.coords) - 1,
        "fine_area_m2": result.fine_area,
        "reduced_vertices": len(result.reduced.exterior.coords) - 1,
        "reduced_area_m2": result.reduced.area,
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
        if row.nve_in_ours is not None and row.ours_in_nve is not None:
            word += (
                f"; NVE's in ours {100 * row.nve_in_ours:.1f} %, ours in NVE's "
                f"{100 * row.ours_in_nve:.1f} %, area ratio {row.area_ratio:.3f}"
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


@app.command()
def station_catchments(
    dem: Annotated[
        list[Path],
        typer.Option("--dem", help="A GeoTIFF DEM, several (repeat --dem), or one directory."),
    ],
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
            help="Lake polygons (GeoJSON, in the river file's CRS): a gauge in a lake, or on "
            "a lake line within 30 m of it, is seeded with the whole lake.",
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
) -> None:
    """Write a catchment per station into --out-dir, as ``catchment --rivers``
    does for one, with results.csv (a row per station) and summary.json; with
    --reference each catchment is compared with its reference polygon and
    classed (match, close, miss, uncertain or refused; increment 29)."""
    for name, value in (("--map-radius", map_radius), ("--reach-up", reach_up)):
        if not (math.isfinite(value) and value > 0.0):
            raise typer.BadParameter("must be a finite number above 0", param_hint=name)
    if outline_tolerance is not None and not (
        math.isfinite(outline_tolerance) and outline_tolerance >= 0.0
    ):
        raise typer.BadParameter("must be finite and at least 0", param_hint="--outline-tolerance")
    target = _destination(out_dir, out_parent, out_dir.name)
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
        from tin_engine.io.station_set import read_lakes as read_nve_lakes  # not 22's read_lakes

        lake_list = _read_beside(read_nve_lakes, lakes, "--lakes", "lakes", crs)
    try:
        repository, _ = repository_for(tuple(dem))
    except (OSError, ValueError) as exc:
        raise typer.BadParameter(
            _words(exc) if isinstance(exc, ValueError) else str(exc), param_hint="--dem"
        ) from exc
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
