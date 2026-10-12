"""The rasputin MCP server over stdio with a test's child runner and heartbeat.

Increment 37 (`docs/increments/37-mcp-server.md`, sections 6.5, 6.6 and 9).
``rasputin mcp`` runs the real CLI as its child with a 5 s heartbeat; T8 to
T10 over stdio need the server with `mcp_fake_child` as the child and, for
T9, a heartbeat of 0.1 s. This file is the server process the SDK's client
spawns for those tests: ``python mcp_stdio_launcher.py CONFIG.json``, where
the config, written by `mcp_driver.stdio_server`, holds

``child``
    the argv prefix the server runs (`mcp_driver.fake_child`'s), or absent
    for the production default;
``jobs``
    the concurrent-jobs cap (default 1);
``heartbeat``
    seconds between heartbeats, or absent for the module's constant.

Everything else is the production server: `build_server` on the real Typer
app, run by the SDK's own ``stdio_server``. Test support only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anyio
from mcp.server.stdio import stdio_server

from tin_engine.cli import app
from tin_engine.mcp import server as server_module


async def main(config: dict[str, object]) -> None:
    if "heartbeat" in config:
        server_module.HEARTBEAT_SECONDS = float(config["heartbeat"])
    options: dict[str, object] = {"jobs": config.get("jobs", 1)}
    if "child" in config:
        options["child"] = config["child"]
    server = server_module.build_server(app, **options)
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())


if __name__ == "__main__":
    anyio.run(main, json.loads(Path(sys.argv[1]).read_text()))
