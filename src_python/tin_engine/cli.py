"""Command-line interface for the rasputin terrain engine.

This module is the **single composition root**. It is the only Python module
that has a path, and joining the file system to the engine is its whole job
(``tin_engine.raster`` also imports ``_core``, to build the one core raster, per
``project_structure.md``): ``viz/`` is written against protocols and never
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

import importlib.metadata
import math
import os
import shlex
import sys
import tempfile
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import numpy as np
import numpy.typing as npt
import typer

from tin_engine._core import (
    ChainRole,
    IndexedMesh2,
    NodedPslg,
    RefineOutcome,
    build_pslg,
    describe,
    node,
    refine,
    sample,
    triangulate,
)
from tin_engine.domain import DomainError, DomainPolygon, read_domain
from tin_engine.elevation import Trimmed, trim
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.grid_domain import default_stride, refine_start_stride, subsample
from tin_engine.io.geotiff import decode_dem
from tin_engine.io.models import GeoTiffError, RasterMeta
from tin_engine.io.ply import write_ply
from tin_engine.io.vtk_legacy import write_vtk
from tin_engine.raster import to_core
from tin_engine.stats import PhaseClock, Refinement, Report, Sizes, _exact, quality, render
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
#: draws identically to the row above it, so declaring all seven features would
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


@app.callback()
def main() -> None:
    """Parallel TIN engine for terrain meshing.

    The callback exists to keep Typer in sub-command mode. With a single
    registered command and no callback, Typer collapses that command into the
    root, and `rasputin version` becomes a usage error.
    """


@app.command()
def version() -> None:
    """Print the installed rasputin version."""
    try:
        typer.echo(importlib.metadata.version("rasputin"))
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        typer.echo("unknown (rasputin is not installed)")


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


#: Ruling 4, verbatim. Written into BOTH files, because the thing it tells a
#: reader -- that this surface is not terrain -- is equally untrue of each, and
#: a person may open either one first.
FLAT_ELEVATION = "none (z=0, --flat)"
FLAT_COMMENT = f"elevation {FLAT_ELEVATION}"

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
        Path | None,
        typer.Option("--dem", help="A GeoTIFF DEM: mesh its extent with z sampled from it."),
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
        typer.Option("--domain-crs", help="The --domain file's CRS, EPSG:n; required for .wkt."),
    ] = None,
    stats: Annotated[
        str | None,
        typer.Option(
            "--stats",
            help="Also write sizes, quality and timings as Markdown to this path (.md "
            "recommended); - prints it on stdout after the path line(s).",
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
    """
    clock = PhaseClock()
    dem_run: _DemMesh | None = None
    if (name is None) == (dem is None):
        raise typer.BadParameter(
            "give a gallery fixture name or --dem PATH, exactly one of the two",
            param_hint="--dem",
        )
    if dem is None and (domain is not None or domain_crs is not None):
        raise typer.BadParameter("applies only with --dem", param_hint="--domain")
    if domain is None and domain_crs is not None:
        raise typer.BadParameter("applies only with --domain", param_hint="--domain-crs")
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

    if dem is not None:
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
        label = dem.stem
        dem_run = _dem_mesh(
            dem, stride, delaunay, snap_spacing, tolerance, clock, domain, domain_crs
        )
        surface_mesh = dem_run.trimmed
        epsg, sentence, described = dem_run.epsg, dem_run.sentence, dem_run.described
        fields = [("crs", f"EPSG:{epsg}"), ("elevation_source", sentence)]
        comments = [f"crs EPSG:{epsg}", f"elevation {sentence}"]
        if described:
            fields.append(("domain", described))
            comments.append(f"domain {described}")
    else:
        assert name is not None
        if stride is not None:
            raise typer.BadParameter("applies only with --dem", param_hint="--stride")
        if tolerance is not None:
            raise typer.BadParameter("applies only with --dem", param_hint="--tolerance")
        if not flat:
            raise typer.BadParameter(
                "a gallery fixture has no elevation source, so z has none; pass --flat "
                "to write z = 0 and say so in the file, or mesh a DEM with --dem"
            )
        label = name
        surface_mesh = _fixture_mesh(name, delaunay, snap_spacing, clock)
        comments = [f"crs {crs}"] if crs else []
        comments.append(FLAT_COMMENT)
        fields = [("crs", crs)] if crs else []
        fields.append(("elevation_source", FLAT_ELEVATION))
        # --crs is unvalidated free text by ruling 5, so the writer's refusals
        # are refusals a person meets by typing, not internal invariants. Turn
        # the writer's ValueError into the usage error it is, in the one place
        # that knows the text came from the command line.
        try:
            write_ply(np.zeros((1, 3)), faces=np.zeros((0, 3)), comments=comments)
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
            )
        ]
    else:
        encoders = [
            lambda: write_ply(
                vertices, faces=surface_mesh.triangles, ascii=not binary, comments=comments
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
        _write_report(clock, report_target, surface_mesh, dem_run, targets)


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


def _write_report(
    clock: PhaseClock,
    target: Path | None,
    trimmed: Trimmed,
    dem_run: _DemMesh | None,
    files: list[Path],
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
        dropped=trimmed.dropped if dem_run else None,
    )
    refinement = dem_run.refinement if dem_run else None
    measured = quality(trimmed.vertices, trimmed.triangles)
    text = render(
        Report(
            command=shlex.join([Path(sys.argv[0]).name, *sys.argv[1:]]),
            sizes=sizes,
            quality=measured,
            refinement=refinement,
            phases=clock.phases(),
            total=total,
            stats_seconds=(time.perf_counter_ns() - t0) / 1e9,
            threads=os.cpu_count() if refinement else None,
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


@dataclass(frozen=True, slots=True)
class _DemMesh:
    """What ``_dem_mesh`` made: the mesh and its file fields, then what
    ``--stats`` reports about the run (None where a row does not apply)."""

    trimmed: Trimmed
    sentence: str
    epsg: int
    described: str
    meta: RasterMeta
    domain_vertices: int | None
    domain_holes: int | None
    start_vertices: int | None
    start_triangles: int | None
    refinement: Refinement | None


def _dem_mesh(
    dem: Path,
    stride: int | None,
    delaunay: bool,
    spacing: float,
    tolerance: float | None,
    clock: PhaseClock,
    domain: Path | None = None,
    domain_crs: str | None = None,
) -> _DemMesh:
    """Decode, subsample, triangulate, sample or refine, and trim.

    Without ``tolerance`` this is increment 12's R6: z sampled bilinearly at
    the stride grid. With it, increment 14's R9: the stride grid is the start
    mesh, refined against the DEM's nodes; with ``domain``, increment 16's R3,
    the polygon's rings are. Returns the mesh, the ``elevation`` sentence for
    the file, the EPSG code, the ``domain`` field (empty without one), and the
    ``--stats`` inputs; ``clock`` gets R5's phases.
    Every refusal is a usage error in the reader's or the engine's own words,
    and no file is written.
    """
    try:
        with clock.phase("decode"), dem.open("rb") as stream:
            tile = decode_dem(stream)
    except OSError as exc:
        raise typer.BadParameter(
            f"cannot read {dem}: {exc.strerror or exc}", param_hint="--dem"
        ) from exc
    except GeoTiffError as exc:
        raise typer.BadParameter(str(exc), param_hint="--dem") from exc

    meta = tile.meta
    described = ""
    domain_vertices = domain_holes = None
    if domain is not None:
        with clock.phase("domain read"):
            try:
                polygon = read_domain(domain, meta, domain_crs)
            except DomainError as exc:
                raise typer.BadParameter(str(exc), param_hint="--domain") from exc
            xy, chains, described = _domain_chains(polygon, domain.name)
        domain_vertices, domain_holes = len(xy), len(polygon.polygon.interiors)
        start = "start domain boundary, boundary z bilinear"
        run = _engine(xy, chains, delaunay, spacing, clock)
    else:
        if stride is not None:
            step = stride
        else:
            step = default_stride(meta) if tolerance is None else refine_start_stride(meta)
        xy, ring = subsample(meta, step)
        start = f"start stride {step}"
        run = _engine(xy, [(ring, ChainRole.Outer, 0)], delaunay, spacing, clock)
    if run.mesh is None or run.noded is None:
        raise typer.BadParameter(f"{dem} has no mesh to write: {run.status}. {run.message}")

    with clock.phase("start mesh: constraint edges"):
        edges, masks = _constraint_arrays(run.mesh, run.noded)
    refinement = None
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
        sentence = f"bilinear from DEM, stride {step}"  # no domain without tolerance
        report = ""
    else:
        t0 = time.perf_counter_ns()
        out = refine(to_core(tile), run.mesh, edges, masks, tolerance=tolerance)
        _refine_phases(clock, (time.perf_counter_ns() - t0) / 1e9, out)
        if not out.ok():
            raise typer.BadParameter(f"{dem}: {out.message}", param_hint="--dem")
        with clock.phase("trim"):
            trimmed = trim(
                vertices=out.vertices,
                triangles=out.triangles,
                edges=out.edges,
                edge_masks=out.masks,
                z=out.z,
                valid=out.valid,
            )
        refinement = Refinement(
            tolerance, out.max_error, out.rounds, out.inserted, out.flips, out.uncovered, out.carved
        )
        sentence = (
            f"refined from DEM nodes, constrained Delaunay, tolerance {_exact(tolerance)} m, "
            f"achieved max error {_exact(out.max_error)} m, {start}, "
            f"{out.uncovered} valid DEM nodes not covered"
        )
        report = (
            f"{out.rounds} rounds, {out.inserted} points inserted, {out.flips} flips, "
            f"{len(trimmed.triangles)} triangles, achieved max error "
            f"{_exact(out.max_error)} m, {out.uncovered} valid DEM nodes not covered, "
            f"{len(run.mesh.triangles)} start triangles, "
            f"{_off_node(np.asarray(run.mesh.vertices), meta)} start vertices off-node, "
        )
    if len(trimmed.triangles) == 0:
        raise typer.BadParameter(
            f"{dem} has no data under any triangle; nothing to write", param_hint="--dem"
        )
    sentence += f", {trimmed.dropped} vertices without data dropped"
    if meta.vertical_unit_assumed:
        sentence += ", vertical unit assumed metres"
    typer.echo(f"{report}{trimmed.dropped} vertices without data dropped", err=True)
    started = refinement is not None
    return _DemMesh(
        trimmed=trimmed,
        sentence=sentence,
        epsg=meta.epsg,
        described=described,
        meta=meta,
        domain_vertices=domain_vertices,
        domain_holes=domain_holes,
        start_vertices=len(run.mesh.vertices) if started else None,
        start_triangles=len(run.mesh.triangles) if started else None,
        refinement=refinement,
    )


def _refine_phases(clock: PhaseClock, seconds: float, out: RefineOutcome) -> None:
    """R5 and R6: ``refine`` and its sub-rows; setup + output is the remainder."""
    inner = (
        ("refine: legalise start", out.legalise_seconds),
        ("refine: scan (parallel)", out.scan_seconds),
        ("refine: split + flip (serial)", out.split_seconds),
    )
    clock.add("refine", seconds)
    for name, part in inner:
        clock.add(name, part)
    clock.add("refine: setup + output", max(0.0, seconds - sum(p for _, p in inner)))


def _domain_chains(
    domain: DomainPolygon, name: str
) -> tuple[npt.NDArray[np.float64], list[tuple[list[int], ChainRole, int]], str]:
    """R3: the outer ring as ``Outer`` and each hole as ``Hole``, mask 0 (U5),
    and the ``domain`` field."""
    rings = [domain.polygon.exterior, *domain.polygon.interiors]
    points: list[tuple[float, float]] = []
    chains: list[tuple[list[int], ChainRole, int]] = []
    for k, ring in enumerate(rings):
        coords = list(ring.coords)[:-1]
        chains.append(
            (
                list(range(len(points), len(points) + len(coords))),
                ChainRole.Outer if k == 0 else ChainRole.Hole,
                0,
            )
        )
        points += coords
    holes = len(rings) - 1
    described = f"{name}, 1 ring {holes} hole{'' if holes == 1 else 's'}, {len(points)} vertices"
    return np.array(points, dtype=np.float64), chains, described


def _off_node(xy: npt.NDArray[np.float64], meta: RasterMeta) -> int:
    """How many of ``xy`` are not a DEM node bit for bit, as ``refine`` classifies."""
    col = np.round((xy[:, 0] - meta.x_min) / meta.delta_x)
    row = np.round((meta.y_max - xy[:, 1]) / meta.delta_y)
    node = (meta.x_min + col * meta.delta_x == xy[:, 0]) & (
        meta.y_max - row * meta.delta_y == xy[:, 1]
    )
    return int(np.count_nonzero(~node))
