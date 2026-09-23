"""Command-line interface for the rasputin terrain engine.

This module is the **single composition root**. It is the only Python module
that imports ``tin_engine._core`` and the only one that has a path, and joining
those two is its whole job: ``viz/`` is written against protocols and never
names a core type, while the core never sees a file, a path or a CRS. Everything
that has to know both sides lives here, per ``project_structure.md``'s rule that
exactly one module constructs a core object.

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
import tempfile
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
    build_pslg,
    describe,
    node,
    triangulate,
)
from tin_engine.features import DEFAULT_VOCABULARY
from tin_engine.io.ply import write_ply
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

PROPERTY_STROKES = tuple(
    PropertyStroke(bit=_BIT_OF[name], token=name) for name in _PRECEDENCE
)

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
    result = build_pslg(np.asarray(fixture.vertices), chains)
    if not result.ok or result.pslg is None:
        return Attempt(
            source=fixture,
            closed_roles=FIXTURE_CLOSED_ROLES,
            mesh=None,
            ok=False,
            status=", ".join(d.error.name for d in result.diagnostics),
            message="; ".join(d.message for d in result.diagnostics),
        )

    noded = node(result.pslg, spacing)
    if not noded.ok() or noded.pslg is None:
        return Attempt(
            source=fixture,
            closed_roles=FIXTURE_CLOSED_ROLES,
            mesh=None,
            ok=False,
            status=noded.status.name,
            message=_band(describe(noded.status), noded.message),
        )

    outcome = triangulate(noded.pslg, delaunay)
    return Attempt(
        source=noded.pslg,
        closed_roles=CORE_CLOSED_ROLES,
        mesh=outcome.mesh if outcome.ok() else None,
        ok=outcome.ok(),
        status=outcome.status.name,
        message=_band(describe(outcome.status), outcome.message),
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
FLAT_COMMENT = "elevation none (z=0, --flat)"


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
    name: Annotated[str, typer.Argument(help="Gallery fixture to write.")],
    out: Annotated[Path, typer.Option("--out", help="Where to write the surface PLY.")],
    flat: Annotated[
        bool, typer.Option("--flat", help="There is no elevation source; write z = 0.")
    ] = False,
    out_edges: Annotated[
        Path | None,
        typer.Option("--out-edges", help="Also write the constraint edges, as a second PLY."),
    ] = None,
    out_parent: Annotated[
        Path | None,
        typer.Option("--out-parent", help="Refuse any output path resolving outside this."),
    ] = None,
    crs: Annotated[
        str, typer.Option("--crs", help="Free text recorded as a header comment. Not validated.")
    ] = "",
    ascii_: Annotated[
        bool, typer.Option("--ascii", help="Write the bodies as text, so `head` can read them.")
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
) -> None:
    """Write a gallery fixture's mesh as PLY.

    Two files, never one holding both element types: MDAL's own caveat is that
    a host application expects either a 1D mesh or a 2D one, so a file carrying
    faces AND edges can load as nothing at all, in silence. ``--out`` gets the
    surface; ``--out-edges`` gets the constraints, or they are not written.
    Both repeat the identical vertex block, which is what makes the two layers
    register on each other when a person loads them side by side.

    Unlike ``draw``, a failed engine run is a non-zero exit and no file. A
    picture of a failure is still a picture and worth producing; there is no
    such thing as a picture of a failed file, so the refusal is reported in the
    engine's own words instead.
    """
    if name not in GALLERY:
        raise typer.BadParameter(f"unknown fixture {name}; the gallery is: {', '.join(GALLERY)}")
    if not flat:
        raise typer.BadParameter(
            "nothing in this tree samples elevation yet, so z has no source; pass --flat "
            "to write z = 0 and say so in the file"
        )

    attempt = _triangulated(GALLERY[name], delaunay=delaunay, spacing=snap_spacing)
    if attempt.mesh is None or not isinstance(attempt.source, NodedPslg):
        raise typer.BadParameter(
            f"{name} has no mesh to write: {attempt.status}. {attempt.message}"
        )

    flat_vertices = np.asarray(attempt.mesh.vertices)
    vertices = np.column_stack([flat_vertices, np.zeros(len(flat_vertices))])
    comments = [f"crs {crs}"] if crs else []
    comments.append(FLAT_COMMENT)

    # --crs is unvalidated free text by ruling 5, so the writer's refusals are
    # refusals a person meets by typing, not internal invariants. A degree sign
    # in a projection string is ordinary and used to exit 1 with a 23-line
    # traceback. Turn the writer's ValueError into the usage error it is, in
    # the one place that knows the text came from the command line.
    try:
        write_ply(np.zeros((1, 3)), faces=np.zeros((0, 3)), comments=comments)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--crs") from exc

    surface = _destination(out, out_parent, name)
    if out_edges is not None:
        # Resolved, because two spellings of one path are still one file. Both
        # writes succeed, the second overwrites the first, the command echoes
        # two paths and exits 0 -- the caller has lost the surface they asked
        # for and nothing said so.
        constraints_target = _destination(out_edges, out_parent, name)
        if constraints_target == surface:
            raise typer.BadParameter(
                f"--out and --out-edges both resolve to {surface}; "
                "the second would overwrite the first",
                param_hint="--out-edges",
            )
    surface.write_bytes(
        write_ply(
            vertices,
            faces=np.asarray(attempt.mesh.triangles),
            ascii=ascii_,
            comments=comments,
        )
    )
    typer.echo(f"{surface}")

    if out_edges is not None:
        edges, masks = _constraint_arrays(attempt.mesh, attempt.source)
        constraints = _destination(out_edges, out_parent, name)
        constraints.write_bytes(
            write_ply(
                vertices,
                edges=edges,
                edge_properties=masks,
                ascii=ascii_,
                comments=comments,
            )
        )
        typer.echo(f"{constraints}")
