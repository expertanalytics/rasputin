"""One user's folder: outputs stored under their sha256, the bytes it holds, and the sweep.

``<user>/work`` is where the user's calls write; ``<user>/store`` holds
``<sha256><suffix>`` files that only the server writes, so no call can change
the bytes under a link already handed out. Walks follow no link and remove
only regular files. The hashing and copying run in a worker thread.
"""

from __future__ import annotations

import hashlib
import os
import re
import secrets
import shutil
import stat
from dataclasses import dataclass
from pathlib import Path

import anyio

__all__ = ["Store", "Stored"]

#: A suffix kept on a stored name; any other is dropped.
SUFFIX = re.compile(r"\.[a-z0-9]{1,8}")


@dataclass(frozen=True)
class Stored:
    sha256: str
    size: int
    path: Path


def _regular(top: Path) -> list[tuple[Path, os.stat_result]]:
    found = []
    for folder, _, names in os.walk(top, followlinks=False):
        for name in names:
            path = Path(folder) / name
            try:
                info = path.lstat()
            except OSError:  # removed while the walk ran
                continue
            if stat.S_ISREG(info.st_mode):
                found.append((path, info))
    return found


class Store:
    """The store of the user whose folder is ``user_dir``."""

    def __init__(self, user_dir: Path) -> None:
        self.user_dir = Path(user_dir)

    def _put(self, path: Path) -> Stored:
        folder = self.user_dir / "store"
        folder.mkdir(parents=True, exist_ok=True)
        partial = folder / f".partial-{secrets.token_hex(8)}"
        try:  # hash the copy, so the stored bytes are the hashed ones
            shutil.copyfile(path, partial)
            with partial.open("rb") as f:
                digest = hashlib.file_digest(f, "sha256").hexdigest()
            suffix = path.suffix.lower() if SUFFIX.fullmatch(path.suffix.lower()) else ""
            target = folder / f"{digest}{suffix}"
            os.replace(partial, target)  # the same bytes again: replaced whole, mtime now
        finally:
            partial.unlink(missing_ok=True)
        return Stored(digest, target.stat().st_size, target)

    async def put(self, path: Path) -> Stored:
        """``path`` copied to ``store/<sha256><suffix>``; the work copy stays."""
        return await anyio.to_thread.run_sync(self._put, Path(path))

    def _receive(self, name: str, data: bytes) -> Path:
        folder = self.user_dir / "work" / "uploads"
        folder.mkdir(parents=True, exist_ok=True)
        if folder.resolve() != folder:
            raise ValueError(f"{folder} resolves to {folder.resolve()}, not to itself")
        partial, target = folder / f".part-{secrets.token_hex(8)}", folder / name
        try:
            partial.write_bytes(data)
            os.replace(partial, target)
        finally:
            partial.unlink(missing_ok=True)
        return target

    async def receive(self, name: str, data: bytes) -> Path:
        """``data`` written whole to ``work/uploads/<name>``, replacing an earlier file of
        that name; a `ValueError` when that folder resolves anywhere but itself."""
        return await anyio.to_thread.run_sync(self._receive, name, data)

    async def usage(self) -> int:
        """The bytes of the regular files under the user's folder."""
        found = await anyio.to_thread.run_sync(_regular, self.user_dir)
        return sum(info.st_size for _, info in found)

    def sweep_now(self, now: float, keep_seconds: float) -> None:
        for top in (self.user_dir / "work", self.user_dir / "store"):
            for path, info in _regular(top):
                if now - info.st_mtime > keep_seconds:
                    path.unlink(missing_ok=True)

    async def sweep(self, now: float, keep_seconds: float) -> None:
        """Remove the regular files in ``work`` and ``store`` older than ``keep_seconds``."""
        await anyio.to_thread.run_sync(self.sweep_now, now, keep_seconds)
