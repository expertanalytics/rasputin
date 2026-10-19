"""Probe: does the MCP SDK's HTTP auth tell a tool handler which user called it?

Two static tokens, a `TokenVerifier` comparing in constant time, the low-level
`Server` served by uvicorn on a free loopback port; two clients call one tool
at once, then a request with no token and one with a wrong token.
Run with a venv that has the `mcp` extra: `python -I http_auth_probe.py`.
"""

import asyncio
import hmac
import socket
import threading
import time
from typing import Any

import httpx2
import mcp_types as types
import uvicorn
from mcp.client import Client
from mcp.client.streamable_http import streamable_http_client
from mcp.server.auth.middleware.auth_context import get_access_token
from mcp.server.auth.provider import AccessToken
from mcp.server.auth.settings import AuthSettings
from mcp.server.lowlevel import Server

TOKENS = {"ola": "a" * 40, "kari": "b" * 40}
ACCEPT = {"Accept": "application/json, text/event-stream"}
LIST = {"jsonrpc": "2.0", "id": 1, "method": "tools/list"}


class Verifier:
    async def verify_token(self, token: str) -> AccessToken | None:
        hit = None
        for user, known in TOKENS.items():
            if hmac.compare_digest(token.encode(), known.encode()):
                hit = user
        return AccessToken(token=token, client_id=hit, scopes=[]) if hit else None


async def list_tools(ctx: Any, params: Any) -> types.ListToolsResult:
    return types.ListToolsResult(tools=[types.Tool(name="who", input_schema={"type": "object"})])


async def call_tool(ctx: Any, params: Any) -> types.CallToolResult:
    request = ctx.request
    user = request.scope.get("user") if request is not None else None
    var = get_access_token()
    await asyncio.sleep(0.3)  # both calls in flight at once
    text = (
        f"req={type(request).__name__} user={getattr(user, 'username', None)} "
        f"cv={var.client_id if var else None}"
    )
    return types.CallToolResult(content=[types.TextContent(type="text", text=text)])


def start() -> str:
    server = Server("probe", on_list_tools=list_tools, on_call_tool=call_tool)
    auth = AuthSettings(issuer_url="http://127.0.0.1", resource_server_url=None)
    app = server.streamable_http_app(auth=auth, token_verifier=Verifier())
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    options = {"app": app, "host": "127.0.0.1", "port": port, "log_level": "warning"}
    threading.Thread(target=uvicorn.run, kwargs=options, daemon=True).start()
    time.sleep(1.5)
    return f"http://127.0.0.1:{port}/mcp"


async def one(url: str, user: str) -> tuple[str, str]:
    http = httpx2.AsyncClient(headers={"Authorization": f"Bearer {TOKENS[user]}"})
    async with Client(streamable_http_client(url, http_client=http)) as client:
        result = await client.call_tool("who", {})
        block = result.content[0]
        return user, block.text if isinstance(block, types.TextContent) else ""


async def main() -> None:
    url = start()
    print(await asyncio.gather(one(url, "ola"), one(url, "kari")))
    bare = httpx2.post(url, json=LIST, headers=ACCEPT)
    print("no token:", bare.status_code, bare.headers.get("www-authenticate"))
    wrong = httpx2.post(url, json=LIST, headers={**ACCEPT, "Authorization": "Bearer " + "c" * 40})
    print("bad token:", wrong.status_code)


asyncio.run(main())
