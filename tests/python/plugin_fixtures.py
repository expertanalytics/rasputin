"""Installed-looking plugins for the `rasputin.commands` hook. Test support only.

`docs/increments/51-open-core.md`, *9.1 Public (P1 unless marked)*.
`write_plugin` lays out what an installed distribution leaves on ``sys.path``:
one module and one ``<name>-<version>.dist-info`` directory holding
``METADATA`` and ``entry_points.txt``. Put that directory on ``sys.path`` (or
``PYTHONPATH`` for a child process) and ``importlib.metadata`` finds it the
way it finds a real install: discovery is not replaced by anything.

Each variant changes one thing: the module's body, the entry point's
target, or the distribution's name.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

#: The plugin of the design: one Typer command, ``hello``, that prints ``hello``.
HELLO = """\
import typer

app = typer.Typer()


@app.command()
def hello() -> None:
    typer.echo("hello")
"""

#: A plugin whose one command has a core command's name.
MESH = """\
import typer

app = typer.Typer()


@app.command("mesh")
def mesh() -> None:
    \"\"\"The plugin's mesh, which must never replace the core's.\"\"\"
    typer.echo("plugin mesh ran")
"""

#: A plugin whose one command writes a line to a ``writable=True`` ``--out``, by
#: itself, without `cli_kit.destination`: the MCP server's output root is its
#: only fence.
SCRIBBLE = """\
from pathlib import Path
from typing import Annotated

import typer

app = typer.Typer()


@app.command()
def scribble(out: Annotated[Path, typer.Option("--out", writable=True)]) -> None:
    out.write_text("a line\\n")
"""


def scribble_defaulting_to(default: Path) -> str:
    """`SCRIBBLE` with ``--out`` optional, defaulting to ``default``: a call that
    leaves ``out`` unset writes there."""
    return SCRIBBLE.replace(
        'typer.Option("--out", writable=True)])',
        f'typer.Option("--out", writable=True)] = Path({str(default)!r}))',
    )


#: A plugin whose module cannot be imported.
IMPORT_ERROR = 'raise ImportError("fakeplug is broken on purpose")\n'

#: A plugin whose module raises something other than ImportError on import.
RUNTIME_ERROR = 'raise RuntimeError("fakeplug fails at import on purpose")\n'

#: A plugin whose entry point names something that is not a ``typer.Typer``.
NOT_TYPER = 'app = "not a Typer app"\n'


def hello_from(dist: str) -> str:
    """A ``hello`` command that says which distribution it came from."""
    return HELLO.replace('typer.echo("hello")', f'typer.echo("hello from {dist}")')


def write_plugin(
    root: Path,
    *,
    dist: str = "fakeplug",
    version: str = "0.1",
    module: str = "fakeplug",
    body: str = HELLO,
    entry: str = "hello",
    target: str = "app",
) -> Path:
    """One distribution under ``root``, with one ``rasputin.commands`` entry
    point ``entry = module:target``; returns ``root``, the directory to put
    on the path."""
    root.mkdir(parents=True, exist_ok=True)
    (root / f"{module}.py").write_text(body)
    info = root / f"{dist}-{version}.dist-info"
    info.mkdir()
    (info / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {dist}\nVersion: {version}\n")
    (info / "entry_points.txt").write_text(f"[rasputin.commands]\n{entry} = {module}:{target}\n")
    return root


#: Variables that make Typer force rich, coloured output whatever ``NO_COLOR``
#: says (``typer/rich_utils.py``); GitHub Actions sets the first.
FORCES_COLOUR = ("GITHUB_ACTIONS", "FORCE_COLOR", "PY_COLORS")

#: An ANSI CSI escape sequence, such as the bold ``\x1b[1m`` rich writes.
ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def plain(text: str) -> str:
    """``text`` without ANSI escapes, so a parser never depends on the terminal."""
    return ANSI.sub("", text)


def child_env(*path: Path) -> dict[str, str]:
    """This process's environment for a child: ``path`` ahead of the inherited
    ``PYTHONPATH``, wide columns so help never wraps a name, and no colour,
    whether or not the tests run under CI."""
    env = dict(os.environ, COLUMNS="200", NO_COLOR="1")
    for name in FORCES_COLOUR:
        env.pop(name, None)
    entries = [str(p) for p in path]
    if env.get("PYTHONPATH"):
        entries.append(env["PYTHONPATH"])
    if entries:
        env["PYTHONPATH"] = os.pathsep.join(entries)
    else:
        env.pop("PYTHONPATH", None)
    return env


def run_rasputin(
    *args: str, path: list[Path] | None = None, timeout: float = 60.0
) -> subprocess.CompletedProcess[str]:
    """``python -m tin_engine *args`` in a child, in ``child_env(*path)``.

    The timeout is a hang guard only: ``--help`` measured about 0.3 s.
    """
    env = child_env(*(path or []))
    return subprocess.run(
        [sys.executable, "-m", "tin_engine", *args],
        capture_output=True,
        text=True,
        env=env,
        timeout=timeout,
        check=False,
    )


def lines_naming(text: str, word: str) -> list[str]:
    """The lines of ``text``, escapes stripped, that contain ``word``."""
    return [line for line in plain(text).splitlines() if word in line]


def write_source_plugin(
    root: Path,
    *,
    dist: str = "fakeplug",
    package: str = "fakeplug_data",
    sources: str | None = None,
    value: str | None = None,
    init: str = "",
) -> Path:
    """One distribution under ``root`` with one ``rasputin.sources`` entry
    point ``<dist> = <value>`` (default: ``package``), and the package
    directory holding ``__init__.py`` (body ``init``) and, unless ``sources``
    is None, ``sources.toml``. Returns ``root``, the directory to put on the path."""
    root.mkdir(parents=True, exist_ok=True)
    directory = root / package
    directory.mkdir()
    (directory / "__init__.py").write_text(init)
    if sources is not None:
        (directory / "sources.toml").write_text(sources, encoding="utf-8")
    info = root / f"{dist}-0.1.dist-info"
    info.mkdir()
    (info / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {dist}\nVersion: 0.1\n")
    entry = package if value is None else value
    (info / "entry_points.txt").write_text(f"[rasputin.sources]\n{dist} = {entry}\n")
    return root
