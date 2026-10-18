"""Increment 69, T4: signed download links, pure and over the download route.

`docs/increments/69-mcp-http.md`, sections 5.5, 6.4 and 9. Invariant critical.
`links` is pure, with ``now`` passed in. Over HTTP, a link comes from a real
call (the fake child writes the file), and the signing key is the suite's
`KEY`, so a test can also sign a link the server accepts: one already expired,
or one to a stored file it has removed.
"""

from __future__ import annotations

import hmac
import os
import time
from collections.abc import Iterator
from typing import Any

import httpx2
import pytest

pytest.importorskip("mcp")

from mcp_http_driver import (
    KEY,
    HttpServer,
    Link,
    call,
    downloads,
    payload,
    rasputin_server,
    raw_get,
    sha256,
)

USER, NAME, FILENAME, EXPIRES = "ola", "ab" * 32 + ".gpkg", "mesh.gpkg", 1_900_000_000


def links_module() -> Any:
    from tin_engine.mcp import links

    return links


def signed(**changes: Any) -> str:
    parts = {"key": KEY, "user": USER, "name": NAME, "filename": FILENAME, "expires": EXPIRES}
    parts.update(changes)
    return str(links_module().sign(**parts))


def valid(sig: str, now: int, **changes: Any) -> bool:
    parts = {"key": KEY, "user": USER, "name": NAME, "filename": FILENAME, "expires": EXPIRES}
    parts.update(changes)
    return bool(links_module().valid(**parts, sig=sig, now=now))


# ---------------------------------------------------------------- pure


class TestSignAndCheck:
    def test_a_signed_link_checks_out_one_second_before_expiry(self) -> None:
        assert valid(signed(), now=EXPIRES - 1)

    def test_at_the_expiry_second_it_is_expired(self) -> None:
        """The check is ``now < expires``."""
        assert not valid(signed(), now=EXPIRES)

    @pytest.mark.parametrize(
        "change",
        [
            {"user": "kari"},
            {"name": "cd" * 32 + ".gpkg"},
            {"filename": "other.gpkg"},
            {"expires": EXPIRES + 1},
            {"key": bytes(32)},
        ],
        ids=["user", "name", "filename", "expiry", "key"],
    )
    def test_any_one_part_changed_fails(self, change: dict[str, Any]) -> None:
        assert not valid(signed(), now=EXPIRES - 3600, **change)

    def test_the_signature_is_hex(self) -> None:
        sig = signed()
        assert sig == sig.lower() and int(sig, 16) >= 0 and len(sig) == 64

    def test_the_signature_is_compared_once_through_hmac_compare_digest(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        links = links_module()
        calls: list[tuple[Any, Any]] = []
        original = hmac.compare_digest

        def recording(a: Any, b: Any) -> bool:
            calls.append((a, b))
            return original(a, b)

        monkeypatch.setattr(links.hmac, "compare_digest", recording)
        sig = signed()
        assert valid(sig, now=EXPIRES - 1)
        assert len(calls) == 1, calls
        assert sig in calls[0] or sig.encode() in calls[0], calls


# ---------------------------------------------------------------- over HTTP

DATA = payload(64 * 1024 + 7, seed=4)


@pytest.fixture(scope="module")
def server(tmp_path_factory: pytest.TempPathFactory) -> Iterator[HttpServer]:
    root = tmp_path_factory.mktemp("links")
    with rasputin_server(root) as running:
        yield running


@pytest.fixture
async def link(server: HttpServer, tmp_path: Any) -> Link:
    """A fresh link to a stored ``t4-<n>.gpkg`` holding `DATA`, from a real call."""
    source = tmp_path / "payload.gpkg"
    source.write_bytes(DATA)
    out = f"t4-{tmp_path.name}.gpkg"
    server.set_spec(copy=[[str(source), out]], stdout=[out])
    result = await call(server, "ola", "version", {})
    assert not result.is_error, result
    (found,) = downloads(result)
    return Link.parse(str(found["url"]))


def get(server: HttpServer, target: str, **headers: str) -> httpx2.Response:
    return httpx2.get(f"{server.base}{target}", headers=headers, timeout=30)


class TestTheRoute:
    def test_a_stored_file_downloads_with_its_name_and_type(
        self, server: HttpServer, link: Link
    ) -> None:
        response = get(server, link.target())
        assert response.status_code == 200
        assert response.content == DATA
        assert link.name == f"{sha256(DATA)}.gpkg"
        disposition = response.headers["content-disposition"]
        assert disposition.startswith("attachment") and link.fn in disposition, disposition
        assert response.headers["content-type"].startswith("application/geopackage+sqlite3")
        assert response.headers["cache-control"] == "private, no-store"
        assert response.headers["referrer-policy"] == "no-referrer"

    def test_the_link_is_signed_with_the_servers_key(self, link: Link) -> None:
        assert valid(
            link.sig,
            now=int(link.exp) - 1,
            user=link.user,
            name=link.name,
            filename=link.fn,
            expires=int(link.exp),
        )

    def test_a_range_request_gets_206_and_ten_bytes(self, server: HttpServer, link: Link) -> None:
        response = get(server, link.target(), Range="bytes=0-9")
        assert response.status_code == 206
        assert response.content == DATA[:10]

    def test_a_changed_signature_is_403(self, server: HttpServer, link: Link) -> None:
        sig = ("0" if link.sig[0] != "0" else "1") + link.sig[1:]
        assert get(server, link.target(sig=sig)).status_code == 403

    def test_a_past_expiry_correctly_signed_is_410(self, server: HttpServer, link: Link) -> None:
        past = int(time.time()) - 10
        sig = signed(user=link.user, name=link.name, filename=link.fn, expires=past)
        assert get(server, link.target(exp=str(past), sig=sig)).status_code == 410

    @pytest.mark.parametrize("exp", ["1e9", "-1", "1" * 13], ids=["float", "negative", "13 digits"])
    def test_a_malformed_expiry_is_404(self, server: HttpServer, link: Link, exp: str) -> None:
        assert get(server, link.target(exp=exp)).status_code == 404

    def test_an_unknown_user_is_404(self, server: HttpServer, link: Link) -> None:
        sig = signed(user="nobody", name=link.name, filename=link.fn, expires=int(link.exp))
        assert get(server, link.target(user="nobody", sig=sig)).status_code == 404

    @pytest.mark.parametrize("name", ["../x", "%2e%2e"], ids=["dot-dot raw", "dot-dot encoded"])
    def test_a_traversing_name_sent_raw_is_404(
        self, server: HttpServer, link: Link, name: str
    ) -> None:
        """Sent with `http.client`: httpx2 would normalise ``/files/ola/../x`` to
        ``/files/x`` before sending."""
        query = link.target().split("?", 1)[1]
        assert raw_get(server, f"/files/{link.user}/{name}?{query}") == 404

    def test_a_stored_file_removed_by_hand_is_404(self, server: HttpServer, link: Link) -> None:
        (server.store(link.user) / link.name).unlink()
        assert get(server, link.target()).status_code == 404

    def test_a_stored_file_replaced_by_a_link_is_404(self, server: HttpServer, link: Link) -> None:
        stored = server.store(link.user) / link.name
        stored.unlink()
        work_copy = next(server.work(link.user).glob("t4-*.gpkg"))
        os.symlink(work_copy, stored)
        try:
            assert get(server, link.target()).status_code == 404
        finally:
            stored.unlink()  # the next test's call stores the file again

    def test_a_bad_filename_is_404(self, server: HttpServer, link: Link) -> None:
        fn = "a b.gpkg"
        sig = signed(user=link.user, name=link.name, filename=fn, expires=int(link.exp))
        assert get(server, link.target(fn=fn, sig=sig)).status_code == 404
