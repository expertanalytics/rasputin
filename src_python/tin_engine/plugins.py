"""The ``rasputin.commands`` plugin hook: installed packages that add commands.

A distribution declares an entry point in the group ``rasputin.commands`` whose
value is a ``typer.Typer``; each of that app's commands becomes a top-level
``rasputin`` command. A plugin never replaces a core command, and a broken
plugin never stops the core: each refusal is one line on stderr naming the
distribution, and the core runs on.
"""

from __future__ import annotations

import copy
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from importlib.metadata import entry_points

import typer
from typer.main import get_command_name
from typer.models import CommandInfo

from tin_engine import cli

GROUP = "rasputin.commands"


@dataclass(frozen=True, slots=True)
class Plugin:
    """One entry point of the group: who declared it, and its app or why there is none."""

    dist: str
    version: str
    entry: str
    app: typer.Typer | None
    error: str | None = None


def discover() -> tuple[Plugin, ...]:
    """One record per entry point in the group, in distribution-name order.

    Never raises: a plugin that fails to import, or whose value is not a
    ``typer.Typer``, is recorded with the reason.
    """
    found = []
    for ep in entry_points(group=GROUP):
        dist = ep.dist.name if ep.dist else ep.module
        version = ep.dist.version if ep.dist else "unknown"
        try:
            loaded = ep.load()
        except Exception as exc:  # a broken plugin is reported, never raised
            found.append(Plugin(dist, version, ep.name, None, f"{type(exc).__name__}: {exc}"))
            continue
        if isinstance(loaded, typer.Typer):
            found.append(Plugin(dist, version, ep.name, loaded))
        else:
            reason = f"{ep.value} is a {type(loaded).__name__}, not a typer.Typer"
            found.append(Plugin(dist, version, ep.name, None, reason))
    return tuple(sorted(found, key=lambda p: (p.dist, p.entry)))


def _name(info: CommandInfo) -> str:
    """The command's name as Typer gives it on the command line."""
    if info.name:
        return info.name
    return get_command_name(info.callback.__name__) if info.callback else ""


def application(plugins: Sequence[Plugin] | None = None) -> typer.Typer:
    """A new app: the core's commands, then each plugin's (``None`` means
    `discover`). ``cli.app`` is never changed. Each refusal is one line on stderr."""
    app = copy.copy(cli.app)
    app.registered_commands = list(cli.app.registered_commands)
    app.registered_groups = list(cli.app.registered_groups)
    owner = {_name(info): "" for info in app.registered_commands}
    for plugin in discover() if plugins is None else plugins:
        if plugin.app is None:
            _refuse(f"plugin {plugin.dist} skipped: {plugin.error}")
            continue
        for info in plugin.app.registered_commands:
            name = _name(info)
            if name in owner:
                by = f"plugin {owner[name]}" if owner[name] else "the core"
                _refuse(f"plugin {plugin.dist}: command {name} skipped, {by} has one")
                continue
            owner[name] = plugin.dist
            app.registered_commands.append(info)
    return app


def _refuse(line: str) -> None:
    print(f"rasputin: {line}", file=sys.stderr)


def loaded() -> list[str]:
    """``<distribution> <version>`` for each plugin that loaded, for ``rasputin version``."""
    return [f"{p.dist} {p.version}" for p in discover() if p.app is not None]


#: The app `main` built and runs, so ``mcp`` serves it rather than building a
#: second one (and printing each refusal twice); ``None`` outside `main`.
served: typer.Typer | None = None


def main() -> None:
    """The ``rasputin`` console script: the core and every installed plugin."""
    global served
    served = application()
    served()
