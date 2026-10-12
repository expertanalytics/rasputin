"""Increment 37, T11: ``rasputin mcp`` is a command whether or not the extra is installed.

`docs/increments/37-mcp-server.md`, sections 6.1, 6.8 and 9. The SDK is the
optional extra ``rasputin[mcp]``; without it, ``rasputin mcp`` refuses with
exit 2 and the install line, and every other command is unchanged. This file
does not skip when ``mcp`` is absent: that is the case it tests. The import is
made to fail by planting ``None`` for the SDK's packages in ``sys.modules``,
and the server module an earlier test may have imported is dropped, so the
command's lazy import cannot be satisfied from the cache.
"""

from __future__ import annotations

import re
import sys

import pytest

from cli_driver import USAGE, invoke, plain, runner
from tin_engine.cli import app

#: The SDK's two top-level packages (`mcp` and its types), both absent without the extra.
SDK_PACKAGES = ("mcp", "mcp_types")
SERVER_MODULES = ("tin_engine.mcp.server", "tin_engine.mcp.schema", "tin_engine.mcp")


def test_mcp_refuses_with_the_install_line_when_the_extra_is_absent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    for name in SDK_PACKAGES:
        monkeypatch.setitem(sys.modules, name, None)
    for name in SERVER_MODULES:
        monkeypatch.delitem(sys.modules, name, raising=False)
    code, output = invoke("mcp")
    assert code == USAGE, output
    assert "Traceback" not in output
    assert "the mcp extra is not installed" in output, output
    assert 'pip install "rasputin[mcp]"' in output, output


def test_help_lists_mcp_beside_the_other_commands() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0, result.output
    assert re.search(r"\bmcp\b", plain(result.output)), plain(result.output)
