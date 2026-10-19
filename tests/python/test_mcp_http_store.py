"""Increment 69, T7: the per-user store, its usage, the sweep, and the quota.

`docs/increments/69-mcp-http.md`, sections 4 (``store.py``), 7 and 9. `Store`
on a temporary folder, ``now`` passed in; then the quota over a loopback
server whose quota is 1 MiB.

Whether `Store`'s methods are coroutines (the hashing and copying may run in a
worker thread inside them) is left to the implementation: `settle` awaits a
result only when it is awaitable.
"""

from __future__ import annotations

import hashlib
import inspect
import os
import time
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("mcp")

from mcp_http_driver import HttpServer, call, error_text, payload, rasputin_server

DAY = 86_400
KEEP = 7 * DAY
NOW = 2_000_000_000.0


async def settle(value: Any) -> Any:
    return await value if inspect.isawaitable(value) else value


@pytest.fixture
def user_dir(tmp_path: Path) -> Path:
    folder = tmp_path / "out" / "ola"
    (folder / "work").mkdir(parents=True)
    (folder / "store").mkdir()
    return folder


def store_of(user_dir: Path) -> Any:
    from tin_engine.mcp.store import Store

    return Store(user_dir)


def read(path: Path) -> bytes:
    return Path(path).read_bytes()


def mtime(path: Path) -> float:
    return Path(path).stat().st_mtime


def aged(path: Path, seconds: float) -> Path:
    """``path`` with its mtime ``seconds`` before `NOW`, the link itself for a link."""
    stamp = NOW - seconds
    os.utime(path, (stamp, stamp), follow_symlinks=False)
    return path


class TestPut:
    async def test_put_stores_under_the_sha256_with_the_suffix(self, user_dir: Path) -> None:
        data = payload(3000, seed=1)
        source = user_dir / "work" / "a.gpkg"
        source.write_bytes(data)
        stored = await settle(store_of(user_dir).put(source))
        digest = hashlib.sha256(data).hexdigest()
        assert stored.sha256 == digest
        assert stored.size == len(data)
        assert Path(stored.path) == user_dir / "store" / f"{digest}.gpkg"
        assert read(stored.path) == data
        assert source.read_bytes() == data, "the work copy stays for later calls"

    async def test_the_same_bytes_again_add_no_file_and_touch_the_stored_one(
        self, user_dir: Path
    ) -> None:
        source = user_dir / "work" / "a.gpkg"
        source.write_bytes(payload(3000, seed=2))
        store = store_of(user_dir)
        first = await settle(store.put(source))
        three_days_ago = time.time() - 3 * DAY  # put touches with the real clock
        os.utime(first.path, (three_days_ago, three_days_ago))
        again = await settle(store.put(source))
        assert again.path == first.path
        assert sorted(p.name for p in (user_dir / "store").iterdir()) == [Path(first.path).name]
        assert mtime(first.path) > three_days_ago + DAY, "not touched"


class TestUsage:
    async def test_usage_counts_work_and_store_and_follows_no_link(
        self, user_dir: Path, tmp_path: Path
    ) -> None:
        (user_dir / "work" / "a.bin").write_bytes(b"x" * 1000)
        (user_dir / "work" / "sub").mkdir()
        (user_dir / "work" / "sub" / "b.bin").write_bytes(b"x" * 200)
        (user_dir / "store" / "c.bin").write_bytes(b"x" * 30)
        big = tmp_path / "big.bin"
        big.write_bytes(b"x" * 100_000)
        (user_dir / "work" / "link.bin").symlink_to(big)
        (user_dir / "work" / "dirlink").symlink_to(tmp_path, target_is_directory=True)
        assert await settle(store_of(user_dir).usage()) == 1230


class TestSweep:
    async def test_files_older_than_seven_days_go_younger_ones_stay(self, user_dir: Path) -> None:
        old_work = user_dir / "work" / "old.vtk"
        old_work.write_bytes(b"o")
        old_store = user_dir / "store" / "old.gpkg"
        old_store.write_bytes(b"o")
        young = user_dir / "work" / "young.vtk"
        young.write_bytes(b"y")
        aged(old_work, KEEP + 1)
        aged(old_store, KEEP + 1)
        aged(young, KEEP - 1)
        await settle(store_of(user_dir).sweep(NOW, KEEP))
        assert not old_work.exists() and not old_store.exists()
        assert young.exists()

    async def test_no_directory_and_no_link_target_is_removed(
        self, user_dir: Path, tmp_path: Path
    ) -> None:
        folder = user_dir / "work" / "sub"
        folder.mkdir()
        target = tmp_path / "target.gpkg"
        target.write_bytes(b"t")
        aged(target, KEEP + DAY)
        link = user_dir / "work" / "link.gpkg"
        link.symlink_to(target)
        aged(link, KEEP + DAY)
        outside_dir = tmp_path / "outside"
        outside_dir.mkdir()
        (outside_dir / "old.gpkg").write_bytes(b"o")
        aged(outside_dir / "old.gpkg", KEEP + DAY)
        (user_dir / "work" / "dirlink").symlink_to(outside_dir, target_is_directory=True)
        aged(folder, KEEP + DAY)
        await settle(store_of(user_dir).sweep(NOW, KEEP))
        assert folder.is_dir()
        assert target.exists()
        assert (outside_dir / "old.gpkg").exists(), "the sweep walked through a link"


# ---------------------------------------------------------------- over HTTP

MIB = 1024 * 1024


@pytest.fixture(scope="module")
def server(tmp_path_factory: pytest.TempPathFactory) -> Iterator[HttpServer]:
    with rasputin_server(tmp_path_factory.mktemp("quota"), quota=MIB) as running:
        yield running


class TestQuota:
    async def test_a_call_at_the_quota_is_refused_naming_it_before_the_child(
        self, server: HttpServer
    ) -> None:
        work = server.work("ola")
        work.mkdir(parents=True, exist_ok=True)
        (work / "full.bin").write_bytes(b"x" * MIB)
        count = server.out_root.parent / "quota-count"
        server.set_spec(count_file=str(count))
        try:
            words = error_text(await call(server, "ola", "version", {}))
            assert "quota" in words, words
            assert not count.exists(), "the child ran over the quota"
            refused = await call(server, "ola", "download_link", {"path": "full.bin"})
            assert "quota" in error_text(refused)
        finally:
            (work / "full.bin").unlink()

    async def test_another_user_under_the_quota_is_not_refused(self, server: HttpServer) -> None:
        work = server.work("ola")
        work.mkdir(parents=True, exist_ok=True)
        (work / "full.bin").write_bytes(b"x" * MIB)
        count = server.out_root.parent / "kari-count"
        server.set_spec(count_file=str(count))
        try:
            result = await call(server, "kari", "version", {})
            assert not result.is_error, result
            assert count.exists()
        finally:
            (work / "full.bin").unlink()
