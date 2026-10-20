"""Which paths a tool call may write and read: one call's scope, and the refusal in words.

No child, no state: a `CallScope` is a value, and `refusal` judges a call's
arguments against it. Every path is resolved as the child would open it, from
the scope's working folder, ``~`` expanded, ``..`` and every link followed.
Over HTTP a call reads only inside its read roots; a value with ``/`` or a
leading ``~`` is judged whether or not it names a file, so a refusal never
tells a caller which files exist, and it must be a plain path, since a path
inside another syntax (a PROJ init file named in a CRS) is opened by the
child too.

`docs/increments/69-mcp-http.md`, *5.2 The read limit (HTTP only)*
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anyio

__all__ = ["HTTP_REFUSED", "CallScope", "overlap", "refusal"]

#: The characters of a plain path: a value with ``/`` or a leading ``~`` must be one.
PLAIN = re.compile(r"[A-Za-z0-9._~/-]+")
#: The (tool, flag) pairs refused over HTTP when true, with the reason.
HTTP_REFUSED = {
    ("fetch", "refresh"): "--refresh discards the shared tile cache for every user, so it is "
    "refused over HTTP; give --cache inside your folder for a fresh copy",
}


@dataclass(frozen=True, kw_only=True)
class CallScope:
    """Where one call runs: ``cwd`` None inherits the server's; ``read_roots`` None is no
    read limit; ``refused`` the (tool, flag) pairs refused when set to true."""

    cwd: Path | None
    write_root: Path
    read_roots: tuple[Path, ...] | None
    user: str | None
    refused: frozenset[tuple[str, str]] = frozenset()


def _place(scope: CallScope, value: str) -> Path:
    path = Path(value).expanduser()
    return (path if scope.cwd is None else scope.cwd / path).resolve()


def _judge(
    scope: CallScope,
    tool_name: str,
    arguments: Mapping[str, Any],
    writes: Sequence[tuple[str, str]],
    reads: Sequence[tuple[str, str]],
) -> str | None:
    for tool, flag in sorted(scope.refused):
        if tool == tool_name and arguments.get(flag) is True:
            return f"{tool}: {HTTP_REFUSED.get((tool, flag), f'--{flag} is refused here')}"
    for argument, given in writes:
        if not (resolved := _place(scope, given)).is_relative_to(scope.write_root):
            return (
                f"{argument}: {given} resolves to {resolved}, outside the output root "
                f"{scope.write_root}; write inside it (rasputin mcp --out-root sets the root)"
            )
    if scope.read_roots is None:
        return None
    for argument, given in reads:
        named = "/" in given or given.startswith("~")
        if named and not PLAIN.fullmatch(given):
            return (
                f"{argument}: {given} is not a plain file path (letters, digits and . _ ~ / - "
                "only); give a file by such a path, a CRS by an EPSG code or PROJ string"
            )
        resolved = _place(scope, given)
        if not any(resolved.is_relative_to(root) for root in scope.read_roots) and (
            named or resolved.exists()
        ):
            return (
                f"{argument}: {given} resolves to {resolved}, outside your folder and the "
                "server's data sources; read only from those"
            )
    return None


async def refusal(
    scope: CallScope,
    tool_name: str,
    arguments: Mapping[str, Any],
    writes: Sequence[tuple[str, str]],
    reads: Sequence[tuple[str, str]],
) -> str | None:
    """The first refusal of a call in plain words, None when it may run: a refused flag,
    then a write outside ``write_root``, then a read outside every read root."""
    return await anyio.to_thread.run_sync(_judge, scope, tool_name, arguments, writes, reads)


def overlap(out_root: Path, shared: Sequence[Path]) -> str | None:
    """The refusal of a shared read root that lies inside the output root or contains it,
    since it would open other users' folders; None when none does."""
    for root in shared:
        if root.is_relative_to(out_root) or out_root.is_relative_to(root):
            return (
                f"the shared read root {root} and the output root {out_root} overlap; "
                "give an --out-root apart from the tile cache and the local data sources"
            )
    return None
