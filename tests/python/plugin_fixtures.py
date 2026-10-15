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


def run_rasputin(
    *args: str, path: list[Path] | None = None, timeout: float = 60.0
) -> subprocess.CompletedProcess[str]:
    """``python -m tin_engine *args`` in a child, with ``path`` ahead of the
    inherited ``PYTHONPATH``. Wide columns, so help never wraps a name.

    The timeout is a hang guard only: ``--help`` measured about 0.3 s.
    """
    env = dict(os.environ, COLUMNS="200", NO_COLOR="1")
    entries = [str(p) for p in path or []]
    if env.get("PYTHONPATH"):
        entries.append(env["PYTHONPATH"])
    if entries:
        env["PYTHONPATH"] = os.pathsep.join(entries)
    else:
        env.pop("PYTHONPATH", None)
    return subprocess.run(
        [sys.executable, "-m", "tin_engine", *args],
        capture_output=True,
        text=True,
        env=env,
        timeout=timeout,
        check=False,
    )


def lines_naming(text: str, word: str) -> list[str]:
    """The lines of ``text`` that contain ``word``."""
    return [line for line in text.splitlines() if word in line]
