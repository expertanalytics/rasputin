"""Increment 78, T1 to T7: the server's own links, and uploads, as tool inputs over HTTP.

`docs/increments/78-mcp-inputs.md`, sections 5 and 6. Invariant critical:
T1 and T3. T1 and T2 are pure (`workspace.linked`, `DerivedCommand.substitute`);
T3 to T6 run a loopback server with the fake child (`mcp_fake_child`), which
records the argv it was given (``argv_file``), counts its runs
(``count_file``) and copies files (``copy``); T7 runs the real
``python -m tin_engine``.

Bounds: the upload limit is 2,097,152 bytes after decoding, checked at the
limit and one byte over, as text and as base64 (both base64 forms are
2,796,204 characters, so only the check after decoding tells them apart).
"""

from __future__ import annotations

import base64
import json
import logging
import os
import sys
import time
import zlib
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import httpx2
import pytest

pytest.importorskip("mcp")

import mcp_http_driver
from cli_helpers import UTM33
from geotiff_fixtures import COLS, DELTA_X, DELTA_Y, ROWS, TIE_X, TIE_Y, micro_tiff
from mcp_driver import call as stdio_call
from mcp_driver import fake_child, option_value
from mcp_http_driver import (
    TOKENS,
    HttpServer,
    Link,
    call,
    client,
    downloads,
    error_text,
    payload,
    rasputin_server,
    sha256,
)
from tin_engine.cli import app as cli_app
from tin_engine.mcp import server as mcp_server
from tin_engine.mcp.schema import derive

TOOLS = {d.tool_name: d for d in derive(cli_app)}
MIB = 1024 * 1024
#: The upload limit after decoding (section 5.2).
LIMIT = 2 * MIB
#: base64 of `LIMIT` bytes and of one byte more: the same number of characters.
LIMIT_BASE64 = 2_796_204

PUBLIC = "https://mcp.example.org"
PREFIX = f"{PUBLIC}/files/"
SHA = "0123456789abcdef" * 4
OLAS = f"{PREFIX}ola/{SHA}.gpkg"
QUERY = "?fn=mesh.gpkg&exp=1900000000&sig=" + "ab" * 32

# The refusals of section 5.3, word for word.
NOT_OURS = (
    "{arg}: {value} is not a link this server made; the server opens no other address. "
    "Send the file with the upload tool and pass the link it returns."
)
OTHER_USER = "{arg}: this link names another user's file; a call reads only your own."
GONE = (
    "{arg}: the file this link names is no longer kept (files are kept 7 days after their "
    "last use); upload it again, or rerun the call that made it."
)
URI_OUT = "{arg}: a link names a stored file, which no call writes; give a path in your folder."
BAD_NAME = (
    "name: give a file name of letters, digits, '.', '_' and '-', not starting with '.', "
    "up to 128 characters."
)
ONE_OF = "give the file as text or as base64, one of the two."
TOO_LARGE = (
    "the file is {n} bytes, over the upload limit of 2097152 (2 MiB); a larger input must be "
    "one of the server's data sources."
)
BAD_BASE64 = "base64: not valid base64 (A-Z, a-z, 0-9, + and /, padded with =)."
INSTRUCTION_SENTENCES = (
    "An input may also be a download link this server gave you, of any age, while its file "
    "is kept.",
    "upload sends a file of your own, up to 2 MiB, and returns its link.",
)


def workspace() -> Any:
    from tin_engine.mcp import workspace as module

    return module


# ---------------------------------------------------------------- T1 linked, pure


def _from_workspace() -> bool:
    """Whether the stack holds a frame of ``workspace.py`` and none of coverage's tracer
    (which stats files of its own while any frame runs)."""
    names = []
    frame = sys._getframe(2)
    while frame is not None:
        names.append(frame.f_code.co_filename)
        frame = frame.f_back
    in_workspace = any(n.endswith(f"{os.sep}mcp{os.sep}workspace.py") for n in names)
    return in_workspace and not any(f"{os.sep}coverage{os.sep}" in n for n in names)


@pytest.fixture
def no_disk(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """A stat or listing of the file system reached from `workspace` fails the test:
    `linked` is pure. Calls from anywhere else (pytest, coverage) pass through."""
    for name in ("stat", "lstat", "listdir", "scandir"):
        original = getattr(os, name)

        def guarded(*args: Any, _original: Any = original, **kwargs: Any) -> Any:
            if _from_workspace():
                raise AssertionError(f"linked touched the disk: {args}")
            return _original(*args, **kwargs)

        monkeypatch.setattr(os, name, guarded)
    yield


@pytest.fixture
def out_root(tmp_path: Path) -> Path:
    """An output root that does not exist, so a value judged on disk could not pass."""
    return tmp_path / "absent" / "out"


def http_scope(out_root: Path, user: str = "ola") -> Any:
    work = out_root / user / "work"
    return workspace().CallScope(
        cwd=work,
        write_root=work,
        read_roots=(work, out_root / user / "store"),
        user=user,
        link_prefix=PREFIX,
    )


def stdio_scope(out_root: Path) -> Any:
    return workspace().CallScope(cwd=None, write_root=out_root, read_roots=None, user=None)


def judged(scope: Any, value: str) -> Any:
    return workspace().linked(scope, "domain", value)


ACCEPTED = {
    "no query": OLAS,
    "its own query": OLAS + QUERY,
    "expiry in the past": f"{OLAS}?fn=mesh.gpkg&exp=1&sig=" + "ab" * 32,
    "changed signature": f"{OLAS}?fn=mesh.gpkg&exp=1900000000&sig=" + "cd" * 32,
    "a fragment": OLAS + "#x",
}
NOT_THIS_SERVERS = {
    "host of the same length": f"https://mcp.example.com/files/ola/{SHA}.gpkg",
    "look-alike host": f"https://mcp.example.org.evil/files/ola/{SHA}.gpkg",
    "filesx": f"https://mcp.example.org/filesx/ola/{SHA}.gpkg",
    "http for https": f"http://mcp.example.org/files/ola/{SHA}.gpkg",
    "file URI": "file:///etc/hosts",
    "data URI": "data:text/plain,x",
    "data URI upper case": "DATA:text/plain,x",
    "dot-dot name": f"{PREFIX}ola/../x",
    "encoded dot-dot": f"{PREFIX}ola/%2e%2e",
    "upper-case hex": f"{PREFIX}ola/{SHA.upper()}.gpkg",
    "a path after the name": f"{PREFIX}ola/{SHA}.gpkg/../x",
    "65 hex digits": f"{PREFIX}ola/{SHA}0",
    "encoded slashes after the name": f"{PREFIX}ola/{SHA}.gpkg%2f..%2fx",
    "no user": f"{PREFIX}{SHA}.gpkg",
    "a user out of pattern": f"{PREFIX}Ola/{SHA}.gpkg",
}
PLAIN_VALUES = ("glo30", "EPSG:3006", "in.geojson")


class TestT1Linked:
    @pytest.mark.parametrize("value", list(ACCEPTED.values()), ids=list(ACCEPTED))
    def test_olas_own_link_is_the_stored_file_by_a_relative_path(
        self, out_root: Path, no_disk: None, value: str
    ) -> None:
        ws = workspace()
        found = judged(http_scope(out_root), value)
        assert found == ws.Linked(
            text=f"../store/{SHA}.gpkg", place=out_root / "ola" / "store" / f"{SHA}.gpkg"
        )

    def test_the_suffix_is_optional(self, out_root: Path, no_disk: None) -> None:
        ws = workspace()
        found = judged(http_scope(out_root), f"{PREFIX}ola/{SHA}")
        assert found == ws.Linked(text=f"../store/{SHA}", place=out_root / "ola" / "store" / SHA)

    @pytest.mark.parametrize("value", list(NOT_THIS_SERVERS.values()), ids=list(NOT_THIS_SERVERS))
    def test_a_uri_this_server_did_not_make_is_refused_without_the_disk(
        self, out_root: Path, no_disk: None, value: str
    ) -> None:
        assert judged(http_scope(out_root), value) == NOT_OURS.format(arg="domain", value=value)

    def test_karis_link_is_another_users_file(self, out_root: Path, no_disk: None) -> None:
        karis = f"{PREFIX}kari/{SHA}.gpkg{QUERY}"
        assert judged(http_scope(out_root), karis) == OTHER_USER.format(arg="domain")

    def test_olas_link_in_karis_call_is_another_users_file(
        self, out_root: Path, no_disk: None
    ) -> None:
        assert judged(http_scope(out_root, "kari"), OLAS) == OTHER_USER.format(arg="domain")

    @pytest.mark.parametrize("value", PLAIN_VALUES)
    def test_a_value_that_is_not_a_uri_is_none(
        self, out_root: Path, no_disk: None, value: str
    ) -> None:
        assert judged(http_scope(out_root), value) is None

    @pytest.mark.parametrize(
        "value",
        [*ACCEPTED.values(), *NOT_THIS_SERVERS.values(), *PLAIN_VALUES],
        ids=[*ACCEPTED, *NOT_THIS_SERVERS, *PLAIN_VALUES],
    )
    def test_on_stdio_every_value_is_none(self, out_root: Path, no_disk: None, value: str) -> None:
        assert judged(stdio_scope(out_root), value) is None

    def test_a_scope_has_no_link_prefix_by_default(self, out_root: Path) -> None:
        assert stdio_scope(out_root).link_prefix is None


# ---------------------------------------------------------------- T2 substitute, pure

LINK = f"{PREFIX}ola/{SHA}.geojson"
STORED = f"../store/{SHA}.geojson"


def validated(tool: str, arguments: dict[str, Any]) -> Any:
    return TOOLS[tool].model.model_validate(arguments)


@pytest.mark.parametrize(
    ("tool", "given", "expected"),
    [
        pytest.param(
            "mesh",
            {"out": "m.vtk", "dem": ["glo30"], "domain": LINK, "tolerance": 1.0},
            {"out": "m.vtk", "dem": ["glo30"], "domain": STORED, "tolerance": 1.0},
            id="plain option",
        ),
        pytest.param(
            "mesh",
            {"out": "m.vtk", "dem": ["glo30"], "features": ["a.geojson", LINK, "b.geojson"]},
            {"out": "m.vtk", "dem": ["glo30"], "features": ["a.geojson", STORED, "b.geojson"]},
            id="repeatable item",
        ),
        pytest.param(
            "mesh",
            {"out": "m.vtk", "dem": ["glo30"], "tolerance_near": [LINK, 0.5]},
            {"out": "m.vtk", "dem": ["glo30"], "tolerance_near": [STORED, 0.5]},
            id="path of a tuple",
        ),
        pytest.param(
            "gml_to_gpkg",
            {"source": LINK, "out": "x.gpkg"},
            {"source": STORED, "out": "x.gpkg"},
            id="positional",
        ),
    ],
)
def test_t2_substitute_replaces_the_value_and_keeps_the_rest(
    tool: str, given: dict[str, Any], expected: dict[str, Any]
) -> None:
    command = TOOLS[tool]
    arguments = validated(tool, given)
    before = arguments.model_dump()
    replaced = command.substitute(arguments, {LINK: STORED})
    want = validated(tool, expected)
    assert replaced.model_dump() == want.model_dump()
    assert replaced.model_fields_set == arguments.model_fields_set
    assert command.argv(replaced) == command.argv(want)
    assert arguments.model_dump() == before, "substitute changed the arguments it was given"


def test_t2_a_value_not_in_the_mapping_is_unchanged() -> None:
    command = TOOLS["mesh"]
    arguments = validated("mesh", {"out": "m.vtk", "dem": ["glo30"], "domain": "d.geojson"})
    replaced = command.substitute(arguments, {LINK: STORED})
    assert command.argv(replaced) == command.argv(arguments)
    assert replaced.model_fields_set == arguments.model_fields_set


# ---------------------------------------------------------------- the HTTP server


@pytest.fixture(scope="module")
def server(tmp_path_factory: pytest.TempPathFactory) -> Iterator[HttpServer]:
    with rasputin_server(tmp_path_factory.mktemp("inputs")) as running:
        for user in TOKENS:
            running.work(user).mkdir(parents=True, exist_ok=True)
        yield running


@dataclass(frozen=True)
class Made:
    """A file a call wrote and the server stored: its bytes, sha256 and link."""

    data: bytes
    sha: str
    url: str

    def stored(self, server: HttpServer, user: str = "ola") -> Path:
        return server.store(user) / f"{self.sha}.geojson"


async def made_by_a_call(server: HttpServer, tmp_path: Path, seed: int) -> Made:
    """A 1,024-byte ``a-<seed>.geojson`` written by a fake-child call of Ola's, and its link."""
    data = payload(1024, seed=seed)
    source = tmp_path / "a.geojson"
    source.write_bytes(data)
    out = f"a-{seed}.geojson"
    server.set_spec(copy=[[str(source), out]], stdout=[out])
    result = await call(server, "ola", "version", {})
    assert not result.is_error, result
    (found,) = downloads(result)
    assert found["sha256"] == sha256(data)
    return Made(data, sha256(data), str(found["url"]))


@pytest.fixture
async def made(server: HttpServer, tmp_path: Path, request: pytest.FixtureRequest) -> Made:
    return await made_by_a_call(server, tmp_path, zlib.crc32(request.node.name.encode()))


def mesh_with(domain: str, **more: Any) -> dict[str, Any]:
    return {"out": "m.vtk", "dem": ["glo30"], "domain": domain, "tolerance": 1.0, **more}


def counting(server: HttpServer, tmp_path: Path) -> Path:
    """The next call's child counts its runs in the returned file, and does nothing else."""
    count = tmp_path / "count.txt"
    server.set_spec(count_file=str(count))
    return count


def folder(server: HttpServer, user: str = "ola") -> dict[str, bytes]:
    """Every regular file under the user's folder, by relative path, with its bytes."""
    top = server.out_root / user
    return {
        str(p.relative_to(top)): p.read_bytes()
        for p in sorted(top.rglob("*"))
        if p.is_file() and not p.is_symlink()
    }


def fetch(url: str) -> bytes:
    response = httpx2.get(url, timeout=60)
    assert response.status_code == 200, response
    return response.content


# ---------------------------------------------------------------- T3 a link round trip


class TestT3LinkRoundTrip:
    async def test_the_link_reaches_the_child_as_the_stored_file(
        self, server: HttpServer, made: Made, tmp_path: Path
    ) -> None:
        stored = made.stored(server)
        old = time.time() - 1000.0
        os.utime(stored, (old, old))
        argv_file = tmp_path / "argv.json"
        server.set_spec(
            argv_file=str(argv_file), copy=[[f"../store/{made.sha}.geojson", "b.geojson"]]
        )
        result = await call(server, "ola", "mesh", mesh_with(made.url))
        assert not result.is_error, result
        argv = json.loads(argv_file.read_text())
        assert option_value(argv, "--domain") == f"../store/{made.sha}.geojson"
        assert made.url not in argv
        assert (server.work("ola") / "b.geojson").read_bytes() == made.data
        assert stored.stat().st_mtime > old + 500.0, "the stored file was not touched"

    async def test_a_stored_file_removed_by_hand_is_no_longer_kept(
        self, server: HttpServer, made: Made, tmp_path: Path
    ) -> None:
        made.stored(server).unlink()
        count = counting(server, tmp_path)
        result = await call(server, "ola", "mesh", mesh_with(made.url))
        assert error_text(result) == GONE.format(arg="domain")
        assert not count.exists(), "the child ran"

    async def test_a_stored_file_replaced_by_a_symlink_is_no_longer_kept(
        self, server: HttpServer, made: Made, tmp_path: Path
    ) -> None:
        stored = made.stored(server)
        stored.unlink()
        work_copy = next(server.work("ola").glob("a-*.geojson"))
        os.symlink(work_copy, stored)
        try:
            count = counting(server, tmp_path)
            result = await call(server, "ola", "mesh", mesh_with(made.url))
            assert error_text(result) == GONE.format(arg="domain")
            assert not count.exists(), "the child ran"
        finally:
            stored.unlink()

    async def test_karis_call_with_olas_link_is_refused(
        self, server: HttpServer, made: Made, tmp_path: Path
    ) -> None:
        count = counting(server, tmp_path)
        result = await call(server, "kari", "mesh", mesh_with(made.url))
        assert error_text(result) == OTHER_USER.format(arg="domain")
        assert not count.exists(), "the child ran"

    @pytest.mark.parametrize(
        "value",
        ["https://evil.example/x.geojson", "file:///etc/hosts"],
        ids=["another host", "file URI"],
    )
    async def test_a_uri_not_this_servers_is_refused(
        self, server: HttpServer, tmp_path: Path, value: str
    ) -> None:
        count = counting(server, tmp_path)
        result = await call(server, "ola", "mesh", mesh_with(value))
        assert error_text(result) == NOT_OURS.format(arg="domain", value=value)
        assert not count.exists(), "the child ran"

    @pytest.mark.parametrize("own", [True, False], ids=["own link", "another address"])
    async def test_out_set_to_a_link_is_refused(
        self, server: HttpServer, made: Made, tmp_path: Path, own: bool
    ) -> None:
        count = counting(server, tmp_path)
        out = made.url if own else "https://evil.example/m.vtk"
        result = await call(server, "ola", "mesh", {"out": out, "dem": ["glo30"]})
        assert error_text(result) == URI_OUT.format(arg="out")
        assert not count.exists(), "the child ran"


async def test_t3_a_restarted_server_with_a_new_key_accepts_an_old_link(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "root"
    with rasputin_server(root) as first:
        first.work("ola").mkdir(parents=True, exist_ok=True)
        old = await made_by_a_call(first, tmp_path, seed=78)
    monkeypatch.setattr(mcp_http_driver, "KEY", bytes(range(1, 33)))
    with rasputin_server(root) as second:
        link = Link.parse(old.url)
        renamed = old.url.replace(link.origin, second.base)
        assert httpx2.get(renamed, timeout=30).status_code == 403, "the key did not change"
        argv_file = tmp_path / "argv.json"
        second.set_spec(argv_file=str(argv_file))
        result = await call(second, "ola", "mesh", mesh_with(renamed))
        assert not result.is_error, result
        argv = json.loads(argv_file.read_text())
        assert option_value(argv, "--domain") == f"../store/{old.sha}.geojson"


# ---------------------------------------------------------------- T4 upload

#: 1,024 bytes of UTF-8 GeoJSON text, letters outside ASCII included.
GEOJSON_TEXT = '{"type": "Polygon", "coordinates": [], "name": "æøå"}'
GEOJSON_TEXT += " " * (1024 - len(GEOJSON_TEXT.encode()))


async def upload(server: HttpServer, user: str = "ola", **arguments: Any) -> Any:
    return await call(server, user, "upload", arguments)


def uploaded(result: Any) -> dict[str, Any]:
    """The one `Download` of an upload, checked against its one resource link."""
    assert not result.is_error, result
    body = result.structured_content
    assert isinstance(body, dict), result
    assert sorted(body) == ["bytes", "expires", "path", "sha256", "url"], body
    linked = [b for b in result.content if b.type == "resource_link"]
    assert [b.uri for b in linked] == [body["url"]], result.content
    return body


class TestT4Upload:
    async def test_text_is_stored_and_linked(self, server: HttpServer) -> None:
        data = GEOJSON_TEXT.encode()
        assert len(data) == 1024
        body = uploaded(await upload(server, name="x.geojson", text=GEOJSON_TEXT))
        assert body["path"] == "uploads/x.geojson"
        assert body["sha256"] == sha256(data) and body["bytes"] == 1024
        assert (server.work("ola") / "uploads" / "x.geojson").read_bytes() == data
        assert (server.store("ola") / f"{sha256(data)}.geojson").read_bytes() == data
        assert fetch(body["url"]) == data

    async def test_base64_of_every_byte_value_round_trips(self, server: HttpServer) -> None:
        data = bytes(range(256))
        encoded = base64.b64encode(data).decode()
        body = uploaded(await upload(server, name="bytes.tif", base64=encoded))
        assert body["path"] == "uploads/bytes.tif"
        assert body["sha256"] == sha256(data)
        assert (server.work("ola") / "uploads" / "bytes.tif").read_bytes() == data
        assert fetch(body["url"]) == data

    async def test_text_of_exactly_the_limit_is_accepted(self, server: HttpServer) -> None:
        body = uploaded(await upload(server, name="limit.txt", text="x" * LIMIT))
        assert body["bytes"] == LIMIT

    async def test_base64_of_exactly_the_limit_is_accepted(self, server: HttpServer) -> None:
        encoded = base64.b64encode(payload(LIMIT, seed=1)).decode()
        assert len(encoded) == LIMIT_BASE64
        body = uploaded(await upload(server, name="limit.bin", base64=encoded))
        assert body["bytes"] == LIMIT

    @pytest.mark.parametrize(
        ("arguments", "size"),
        [
            pytest.param({"text": "x" * (LIMIT + 1)}, LIMIT + 1, id="text one byte over"),
            pytest.param(
                {"text": "é" * (LIMIT // 2 + 1)}, LIMIT + 2, id="text over in UTF-8 bytes"
            ),
            pytest.param(
                {"base64": base64.b64encode(payload(LIMIT + 1, seed=2)).decode()},
                LIMIT + 1,
                id="base64 one byte over",
            ),
        ],
    )
    async def test_over_the_limit_is_refused_and_nothing_written(
        self, server: HttpServer, arguments: dict[str, Any], size: int
    ) -> None:
        if "base64" in arguments:
            assert len(arguments["base64"]) == LIMIT_BASE64
        before = folder(server)
        result = await upload(server, name="over.bin", **arguments)
        assert error_text(result) == TOO_LARGE.format(n=size)
        assert folder(server) == before

    @pytest.mark.parametrize(
        ("arguments", "words"),
        [
            pytest.param({"name": "n.txt"}, ONE_OF, id="neither"),
            pytest.param({"name": "n.txt", "text": "a", "base64": "YQ=="}, ONE_OF, id="both"),
            pytest.param({"name": "..", "text": "a"}, BAD_NAME, id="dot-dot"),
            pytest.param({"name": ".x", "text": "a"}, BAD_NAME, id="leading dot"),
            pytest.param({"name": "a/b", "text": "a"}, BAD_NAME, id="slash"),
            pytest.param({"name": "a b", "text": "a"}, BAD_NAME, id="space"),
            pytest.param({"name": "a" * 129, "text": "a"}, BAD_NAME, id="129 characters"),
            pytest.param({"name": "n.bin", "base64": "@@@"}, BAD_BASE64, id="not base64"),
            pytest.param({"name": "n.bin", "base64": "QQ"}, BAD_BASE64, id="unpadded"),
        ],
    )
    async def test_a_bad_upload_is_refused_and_nothing_written(
        self, server: HttpServer, arguments: dict[str, Any], words: str
    ) -> None:
        before = folder(server)
        result = await call(server, "ola", "upload", arguments)
        assert error_text(result) == words
        assert folder(server) == before

    async def test_a_name_of_128_characters_is_accepted(self, server: HttpServer) -> None:
        name = "a" * 124 + ".txt"
        body = uploaded(await upload(server, name=name, text="a"))
        assert body["path"] == f"uploads/{name}"

    async def test_the_same_name_twice_keeps_the_second_and_both_links(
        self, server: HttpServer
    ) -> None:
        first = uploaded(await upload(server, name="twice.geojson", text="first"))
        second = uploaded(await upload(server, name="twice.geojson", text="second"))
        assert (server.work("ola") / "uploads" / "twice.geojson").read_bytes() == b"second"
        assert fetch(first["url"]) == b"first"
        assert fetch(second["url"]) == b"second"

    async def test_an_uploads_link_is_an_input(self, server: HttpServer, tmp_path: Path) -> None:
        body = uploaded(await upload(server, name="in.geojson", text=GEOJSON_TEXT))
        argv_file = tmp_path / "argv.json"
        server.set_spec(argv_file=str(argv_file))
        result = await call(server, "ola", "mesh", mesh_with(body["url"]))
        assert not result.is_error, result
        argv = json.loads(argv_file.read_text())
        assert option_value(argv, "--domain") == f"../store/{body['sha256']}.geojson"

    async def test_upload_is_listed_beside_download_link(self, server: HttpServer) -> None:
        async with client(server, "ola") as c:
            listed = {t.name: t for t in (await c.list_tools()).tools}
        assert "download_link" in listed
        notes = listed["upload"].annotations
        assert notes is not None
        assert (
            notes.read_only_hint,
            notes.destructive_hint,
            notes.idempotent_hint,
            notes.open_world_hint,
        ) == (False, False, True, False)


async def test_t4_a_folder_at_the_quota_takes_no_upload(tmp_path: Path) -> None:
    with rasputin_server(tmp_path, quota=MIB) as small:
        small.work("ola").mkdir(parents=True, exist_ok=True)
        (small.work("ola") / "full.bin").write_bytes(payload(MIB, seed=3))
        before = folder(small)
        result = await upload(small, name="more.txt", text="a")
        words = error_text(result)
        assert words.startswith(f"your folder holds {MIB} bytes, at or over the quota"), words
        assert folder(small) == before


# ---------------------------------------------------------------- T5 stdio unchanged


async def test_t5_stdio_lists_no_upload(tmp_path: Path) -> None:
    from mcp.client import Client

    built = mcp_server.build_server(cli_app, child=fake_child(tmp_path), out_root=tmp_path)
    async with Client(built) as c:
        names = [t.name for t in (await c.list_tools()).tools]
    assert "upload" not in names


async def test_t5_on_stdio_a_link_reaches_the_child_unchanged(tmp_path: Path) -> None:
    argv_file = tmp_path / "argv.json"
    child = fake_child(tmp_path, argv_file=str(argv_file))
    built = mcp_server.build_server(cli_app, child=child, out_root=tmp_path)
    link = f"https://mcp.example.org/files/ola/{SHA}.geojson{QUERY}"
    out = str(tmp_path / "m.vtk")
    result = await stdio_call(built, "mesh", {"out": out, "dem": ["glo30"], "domain": link})
    assert not result.is_error, result
    assert option_value(json.loads(argv_file.read_text()), "--domain") == link


async def test_t5_on_stdio_out_set_to_a_uri_is_not_refused_by_the_link_check(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With the server's working folder inside the output root, the write guard passes
    ``https://x/m.vtk`` as a relative path; only the HTTP check would refuse it."""
    monkeypatch.chdir(tmp_path)
    argv_file = tmp_path / "argv.json"
    child = fake_child(tmp_path, argv_file=str(argv_file))
    built = mcp_server.build_server(cli_app, child=child, out_root=tmp_path)
    result = await stdio_call(built, "mesh", {"out": "https://x/m.vtk", "dem": ["glo30"]})
    assert not result.is_error, result
    assert option_value(json.loads(argv_file.read_text()), "--out") == "https://x/m.vtk"


# ---------------------------------------------------------------- T6 start and log


def test_t6_a_command_named_upload_stops_start(tmp_path: Path) -> None:
    import typer

    from tin_engine.mcp import http as mcp_http

    toy = typer.Typer()

    @toy.command("upload")
    def upload_command(name: str) -> None:
        """A plugin command whose tool name the server needs for itself."""

    @toy.command("other")
    def other() -> None:
        """Another command, so the app keeps its subcommands."""

    settings = mcp_http.HttpSettings(
        host="127.0.0.1",
        port=8765,
        public_url="http://127.0.0.1:8765",
        out_root=tmp_path / "out",
        quota=1024,
    )
    tokens = mcp_http.tokens_from_env({"RASPUTIN_MCP_TOKENS": mcp_http_driver.TOKENS_VALUE})
    with pytest.raises(ValueError, match="upload"):
        mcp_http.build_app(toy, settings, tokens)


#: The server's loggers: `tin_engine`'s, uvicorn's and the MCP SDK's server side, never
#: the test's own client, whose request lines and debug records the server cannot govern.
SERVER_LOGGERS = ("tin_engine", "uvicorn", "mcp.server", "mcp.shared")


def is_the_servers(record: logging.LogRecord) -> bool:
    return any(record.name == n or record.name.startswith(f"{n}.") for n in SERVER_LOGGERS)


@pytest.mark.parametrize("sentence", INSTRUCTION_SENTENCES)
def test_t6_the_http_instructions_name_links_and_upload(sentence: str) -> None:
    from tin_engine.mcp import http as mcp_http

    assert sentence in mcp_http.HTTP_INSTRUCTIONS


async def test_t6_an_upload_logs_one_line_without_its_text_or_signature(
    server: HttpServer, caplog: pytest.LogCaptureFixture
) -> None:
    marker = "MARKER-the-uploaded-text-7f3a"
    text = json.dumps({"type": "Polygon", "coordinates": [], "name": marker})
    caplog.set_level(logging.DEBUG)
    body = uploaded(await upload(server, name="logged.geojson", text=text))
    sig = Link.parse(body["url"]).sig
    ours = [r.getMessage() for r in caplog.records if r.name == "tin_engine.mcp.http"]
    lines = [t for t in ours if t.startswith("upload ")]
    assert lines == [
        f"upload user=ola file={sha256(text.encode())[:12]} bytes={len(text.encode())}"
    ], ours
    servers = [r.getMessage() for r in caplog.records if is_the_servers(r)]
    assert servers, "no record from the server's loggers was captured"
    for message in servers:
        assert marker not in message, message
        assert sig not in message, message


# ---------------------------------------------------------------- T7 the whole exchange

#: The micro tile's node extent (3 x 4 nodes, pixel-is-point) as a GeoJSON polygon.
X_MIN, X_MAX = TIE_X, TIE_X + DELTA_X * (COLS - 1)
Y_MIN, Y_MAX = TIE_Y - DELTA_Y * (ROWS - 1), TIE_Y
POLYGON = json.dumps(
    {
        "type": "Polygon",
        "coordinates": [
            [[X_MIN, Y_MIN], [X_MAX, Y_MIN], [X_MAX, Y_MAX], [X_MIN, Y_MAX], [X_MIN, Y_MIN]]
        ],
        "crs": {"type": "name", "properties": {"name": UTM33}},
    }
)
#: `POLYGON`'s size in bytes, recorded so a change to the fixture shows here.
POLYGON_BYTES = 233


async def test_t7_upload_a_polygon_then_mesh_with_its_link(tmp_path: Path) -> None:
    from vtkread import read_vtk

    assert len(POLYGON.encode()) == POLYGON_BYTES
    with rasputin_server(tmp_path, child=list(mcp_server.CHILD)) as real:
        real.work("ola").mkdir(parents=True, exist_ok=True)
        (real.work("ola") / "micro.tif").write_bytes(micro_tiff().getvalue())
        polygon = uploaded(await upload(real, name="domain.geojson", text=POLYGON))
        arguments = {"dem": ["micro.tif"], "domain": polygon["url"], "out": "m.vtk"}
        result = await call(real, "ola", "mesh", {**arguments, "tolerance": 1.0})
        assert not result.is_error, result
        assert result.structured_content is not None
        assert result.structured_content["exit_code"] == 0
        vtk = [d for d in downloads(result) if str(d["path"]).endswith(".vtk")]
        assert [d["path"] for d in vtk] == ["m.vtk"]
        assert read_vtk(fetch(str(vtk[0]["url"]))) is not None
