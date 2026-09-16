"""Command-line interface for the rasputin terrain engine."""

from __future__ import annotations

import importlib.metadata

import typer

app = typer.Typer(
    name="rasputin",
    help="Parallel TIN engine for terrain meshing.",
    no_args_is_help=True,
)


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
