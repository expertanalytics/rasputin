"""The command-line helpers the core's commands share with a plugin's.

The ``--threads`` option, the output-path boundary, the resolved threads and
memory budget, and a refusal's own words. Each is what the ``rasputin``
commands use, so a plugin's commands behave as the core's do.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Annotated

import typer
from pydantic import ValidationError

from tin_engine.resources import VARIABLE, Resources, resources

#: The ``--threads`` option every parallel command takes.
THREADS = Annotated[
    int | None,
    typer.Option(
        "--threads",
        min=0,
        metavar="N",
        help="Threads for every parallel step (refinement, decoding, resampling, "
        "reprojection); 0 means every core, a container's CPU quota counted. Default: "
        "$RASPUTIN_THREADS, else every core.",
    ),
]


def destination(out: Path | None, out_parent: Path | None, name: str, hint: str = "--out") -> Path:
    """Where to write, having refused every path this process should not touch; a
    refusal names the option, ``hint``.

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
        raise typer.BadParameter(f"{out} is a symlink; refusing to follow it", param_hint=hint)
    target = out.resolve()
    if out_parent is not None:
        permitted = out_parent.resolve()
        if not target.is_relative_to(permitted):
            raise typer.BadParameter(
                f"{target} is outside the permitted parent {permitted}", param_hint=hint
            )
    if not target.parent.is_dir():
        raise typer.BadParameter(f"{target.parent} is not an existing directory", param_hint=hint)
    return target


def cli_resources(threads: int | None, memory_budget: str | None) -> Resources:
    """``--threads`` and ``--memory-budget`` resolved against the environment,
    built once per command; a refusal is a usage error naming the option or
    the variable."""
    try:
        return resources(threads, memory_budget, os.environ)
    except ValueError as exc:
        hint = VARIABLE if str(exc).startswith(VARIABLE) else "--memory-budget"
        raise typer.BadParameter(str(exc), param_hint=hint) from exc


def words(exc: ValueError) -> str:
    """A refusal's own words: Pydantic's messages without its wrapper."""
    if isinstance(exc, ValidationError):
        return "; ".join(str(e["msg"]) for e in exc.errors())
    return str(exc)
