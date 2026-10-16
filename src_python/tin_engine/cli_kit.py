"""The command-line helpers the core's commands share with a plugin's.

The ``--threads`` option, the output-path boundary, the resolved threads and
memory budget, and a refusal's own words. Each is what the ``rasputin``
commands use, so a plugin's commands behave as the core's do.
"""

from __future__ import annotations

import os
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


def destination(out: Path, hint: str = "--out") -> Path:
    """Where to write ``out``, its missing folders made; a refusal names the
    option, ``hint``.

    A symlinked final part is refused rather than followed: a link that points
    elsewhere inside an allowed folder still overwrites a file the caller never
    named. The path is then resolved, so ``..`` makes the folder it reaches, not
    the one it spells. It fences nothing: a caller that must keep writes inside
    a folder checks that before calling.
    """
    if out.is_symlink():
        raise typer.BadParameter(f"{out} is a symlink; refusing to follow it", param_hint=hint)
    target = out.resolve()
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise typer.BadParameter(
            f"cannot make the folder {target.parent}: {exc.strerror}", param_hint=hint
        ) from exc
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
