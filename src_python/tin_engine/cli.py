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
   then ``triangulate``;
3. hands the **fixture itself** to ``build_scene`` as the ``PslgLike``, with
   ``closed_roles=("outer", "hole")``;
4. passes the engine's own words to ``render_svg`` as text;
5. resolves and validates the output path, and writes the bytes.

Step 3 is not a convenience. A fixture the validator *rejects* has no ``Pslg``
at all -- ``degenerate`` is refused with ``PslgError.DegenerateRing`` before
``triangulate`` is ever called -- so the only ``PslgLike`` that exists for such
a run is the fixture, and drawing it is what keeps the failure presentation from
being a blank page. Step 3 is also why the roles are strings: ``viz/`` cannot
name the enum, so the mapping in :data:`ROLES` is this module's to own.
"""

from __future__ import annotations

import importlib.metadata
import tempfile
from pathlib import Path
from typing import Annotated

import numpy as np
import typer

from tin_engine._core import ChainRole, IndexedMesh2, build_pslg, describe, triangulate
from tin_engine.viz.fixtures import GALLERY, Fixture
from tin_engine.viz.scene import build_scene
from tin_engine.viz.style import SvgStyle
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

#: Which roles name chains that are rings, and therefore carry a closing edge
#: no ``Pslg`` stores. ``scene.py`` cannot work this out -- it may not name the
#: enum -- and without it every ring closure is drawn as a false alarm.
CLOSED_ROLES = ("outer", "hole")

DEFAULT_LABEL_LIMIT = 500


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


def _triangulated(fixture: Fixture, delaunay: bool) -> tuple[IndexedMesh2 | None, bool, str, str]:
    """Run the engine on a fixture and reduce what it said to four values.

    Returns ``(mesh, ok, status, message)``, where the last two are plain text
    for the header band. Two distinct refusals reach this point and neither is
    an exception: the validator can reject the constraint set, in which case
    there is no ``Pslg`` and the words come from every ``PslgDiagnostic`` rather
    than only the first; or the backend can refuse to triangulate a valid one.
    """
    chains = [
        ([int(i) for i in fixture.indices_of(c)], ROLES[chain.role], chain.is_river)
        for c, chain in enumerate(fixture.chains)
    ]
    result = build_pslg(np.asarray(fixture.vertices), chains)
    if not result.ok or result.pslg is None:
        return (
            None,
            False,
            ", ".join(d.error.name for d in result.diagnostics),
            "; ".join(d.message for d in result.diagnostics),
        )
    outcome = triangulate(result.pslg, delaunay)
    mesh = outcome.mesh if outcome.ok() else None
    message = " ".join(part for part in (describe(outcome.status), outcome.message) if part)
    return mesh, outcome.ok(), outcome.status.name, message


def _destination(out: Path | None, out_parent: Path | None, name: str) -> Path:
    """Where to write, having refused every path this process should not touch.

    A hostile boundary, per the python skill. ``resolve()`` first, so that
    ``..`` cannot smuggle a write out of an explicitly permitted parent; refuse
    a symlink outright rather than following it, because a link pointing back
    *inside* the permitted parent passes containment and still overwrites a file
    the user never named; and turn a missing directory into a message rather
    than a traceback.
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
) -> None:
    """Draw a gallery fixture as an SVG.

    Exit code 0 means a picture was produced, which includes the two deliberate
    failure fixtures: a drawn failure is still a drawn picture, and the engine's
    refusal is in the header band where a person can read it.
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
    mesh, ok, status, message = _triangulated(fixture, delaunay)
    scene = build_scene(fixture, mesh, ok=ok, closed_roles=CLOSED_ROLES)
    if labels and len(scene.vertices) > label_limit:
        raise typer.BadParameter(
            f"{len(scene.vertices)} vertices exceeds the label limit of {label_limit}; "
            f"raise --label-limit or drop --labels"
        )

    target = _destination(out, out_parent, name)
    document = render_svg(
        scene,
        SvgStyle(show_vertices=show_vertices),
        title=title,
        status=status,
        message=message,
        labels=labels,
    )
    target.write_text(document, encoding="utf-8")
    typer.echo(f"{target}")
