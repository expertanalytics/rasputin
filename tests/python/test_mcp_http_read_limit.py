"""Increment 69, T3: the read limit, the write guard per user, and the refused flags.

`docs/increments/69-mcp-http.md`, sections 5.2, 5.3, 6.2, 7 and 9. Invariant
critical. Most cases run the server's own pipeline without a server: the real
CLI's tool, its validated arguments, ``DerivedCommand.writes`` and ``.reads``,
then `workspace.refusal` on a temporary tree laid out as the server lays out a
user's folder. The rest go over a loopback server with the fake child: a
refused call never starts the child, and the check runs again after the job
slot is taken.

The tree: ``out/ola/work`` (the caller's cwd and write root), ``out/kari/work``
(another user), ``cache`` and ``datasets/local.tif`` (shared read roots), and
``outside`` (a folder no root contains).
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anyio
import pytest

pytest.importorskip("mcp")

from dataset_fixtures import local_entry, write_sources
from mcp_http_driver import HttpServer, call, error_text, rasputin_server, until_exists
from tin_engine.cli import app
from tin_engine.mcp import server as mcp_server
from tin_engine.mcp.schema import DerivedCommand, derive

TOOLS = {d.tool_name: d for d in derive(app)}
#: A PROJJSON transformation whose grid file parameter names a file by absolute path.
PROJJSON_GRID = json.dumps(
    {
        "$schema": "https://proj.org/schemas/v0.7/projjson.schema.json",
        "type": "Transformation",
        "name": "grid shift",
        "method": {"name": "Geographic3D to GravityRelatedHeight (gtx)"},
        "parameters": [{"name": "Geoid (height correction) model file", "value": "/etc/hosts"}],
    }
)


@dataclass(frozen=True)
class Tree:
    root: Path
    work: Path
    kari: Path
    cache: Path
    local: Path
    outside: Path


@pytest.fixture
def tree(tmp_path: Path) -> Tree:
    root = tmp_path.resolve()
    work, kari = root / "out" / "ola" / "work", root / "out" / "kari" / "work"
    cache, outside = root / "cache", root / "outside"
    for folder in (work, kari, cache, outside):
        folder.mkdir(parents=True)
    local = root / "datasets" / "local.tif"
    local.parent.mkdir()
    local.write_bytes(b"II*\x00")
    (work / "in.geojson").write_text("{}")
    (outside / "secret.geojson").write_text("{}")
    return Tree(root, work, kari, cache, local, outside)


def http_scope(tree: Tree) -> Any:
    from tin_engine.mcp.workspace import CallScope

    return CallScope(
        cwd=tree.work,
        write_root=tree.work,
        read_roots=(tree.work, tree.cache, tree.local),
        user="ola",
        refused=frozenset({("fetch", "refresh")}),
    )


def stdio_scope(tree: Tree) -> Any:
    from tin_engine.mcp.workspace import CallScope

    return CallScope(
        cwd=None, write_root=tree.root, read_roots=None, user=None, refused=frozenset()
    )


async def judge(scope: Any, tool_name: str, arguments: dict[str, Any]) -> str | None:
    """The server's checks on one call, in its order: validate, then `refusal`."""
    from tin_engine.mcp.workspace import refusal

    tool: DerivedCommand = TOOLS[tool_name]
    model = tool.model.model_validate(arguments)
    return await refusal(
        scope, tool_name, model.model_dump(), tool.writes(model), tool.reads(model)
    )


def mesh(**more: Any) -> dict[str, Any]:
    """A mesh call that writes inside ``work`` and reads a bare source name."""
    return {"out": "mesh.vtk", "dem": ["glo30"], **more}


# ------------------------------------------------- reads() itself


class TestReadsFindsEveryValue:
    def test_tuple_elements_items_and_defaults_are_read_values(self) -> None:
        tool = TOOLS["mesh"]
        model = tool.model.model_validate(
            mesh(tolerance_near=["/etc/hosts", 5], features=["a.geojson", "b.geojson"])
        )
        reads = tool.reads(model)
        assert ("tolerance_near", "/etc/hosts") in reads
        assert ("features", "a.geojson") in reads and ("features", "b.geojson") in reads
        assert ("dem", "glo30") in reads
        assert all(value != "5" and value != "5.0" for _, value in reads), reads

    def test_a_writable_path_is_not_a_read(self) -> None:
        tool = TOOLS["mesh"]
        reads = tool.reads(tool.model.model_validate(mesh(stats="s.md")))
        assert all(argument not in {"out", "stats"} for argument, _ in reads), reads


# ------------------------------------------------- refused


class TestRefused:
    @pytest.mark.parametrize(
        ("argument", "arguments"),
        [
            pytest.param("domain", {"domain": "/etc/hosts"}, id="absolute outside"),
            pytest.param("domain", {"domain": "../../kari/work/x.geojson"}, id="other user"),
            pytest.param("domain", {"domain": "/nonexistent/x.tif"}, id="missing absolute"),
            pytest.param("domain", {"domain": ".."}, id="bare dot-dot"),
            pytest.param(
                "tolerance_near", {"tolerance_near": ["/etc/hosts", 5]}, id="tuple element"
            ),
            pytest.param(
                "features", {"features": ["in.geojson", "/etc/hosts"]}, id="repeatable item"
            ),
        ],
    )
    async def test_a_read_outside_every_root_names_the_argument(
        self, tree: Tree, argument: str, arguments: dict[str, Any]
    ) -> None:
        words = await judge(http_scope(tree), "mesh", mesh(**arguments))
        assert words is not None, f"{arguments} passed the read limit"
        assert argument in words, words

    async def test_a_directory_outside_given_as_dem(self, tree: Tree) -> None:
        words = await judge(http_scope(tree), "mesh", mesh(dem=[str(tree.outside)]))
        assert words is not None and "dem" in words, words

    async def test_a_home_path_is_expanded_before_it_is_judged(
        self, tree: Tree, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        home = tree.root / "home"
        home.mkdir()
        (home / "x").write_text("{}")
        monkeypatch.setenv("HOME", str(home))
        words = await judge(http_scope(tree), "mesh", mesh(domain="~/x"))
        assert words is not None and "domain" in words, words

    async def test_a_link_in_work_pointing_outside_is_followed(self, tree: Tree) -> None:
        (tree.work / "link.geojson").symlink_to(tree.outside / "secret.geojson")
        words = await judge(http_scope(tree), "mesh", mesh(domain="link.geojson"))
        assert words is not None and "domain" in words, words
        assert str(tree.outside / "secret.geojson") in words, "the place it resolves to"

    async def test_a_path_into_another_users_folder_is_refused_before_the_file_exists(
        self, tree: Tree
    ) -> None:
        """The same words whether or not Kari's file exists: a refusal never tells a
        caller which files exist."""
        given = "../../kari/work/out.geojson"
        before = await judge(http_scope(tree), "mesh", mesh(domain=given))
        (tree.kari / "out.geojson").write_text("{}")
        after = await judge(http_scope(tree), "mesh", mesh(domain=given))
        assert before is not None and "domain" in before, before
        assert before == after

    @pytest.mark.parametrize(
        "named", ["/etc/hosts", "/nonexistent/rasputin-init"], ids=["present", "missing"]
    )
    async def test_a_proj_init_with_a_path_is_not_a_plain_path(
        self, tree: Tree, named: str
    ) -> None:
        words = await judge(http_scope(tree), "mesh", mesh(out_crs=f"+init={named}:x"))
        assert words is not None and "not a plain file path" in words, words
        assert "out_crs" in words, words

    async def test_the_proj_init_refusal_does_not_depend_on_the_file(self, tree: Tree) -> None:
        present = await judge(http_scope(tree), "mesh", mesh(out_crs="+init=/etc/hosts:x"))
        missing = await judge(http_scope(tree), "mesh", mesh(out_crs="+init=/nonexistent/y:x"))
        assert present is not None and missing is not None
        assert present.replace("/etc/hosts", "<v>") == missing.replace("/nonexistent/y", "<v>")

    async def test_a_projjson_with_a_grid_path_is_not_a_plain_path(self, tree: Tree) -> None:
        words = await judge(http_scope(tree), "mesh", mesh(out_crs=PROJJSON_GRID))
        assert words is not None and "not a plain file path" in words, words

    async def test_a_file_name_with_a_space_is_not_a_plain_path(self, tree: Tree) -> None:
        (tree.work / "in file.geojson").write_text("{}")
        words = await judge(http_scope(tree), "mesh", mesh(domain="./in file.geojson"))
        assert words is not None and "not a plain file path" in words, words
        assert "domain" in words, words

    async def test_fetch_refresh_is_refused_over_http_naming_the_flag(self, tree: Tree) -> None:
        words = await judge(http_scope(tree), "fetch", {"source": "glo30", "refresh": True})
        assert words is not None and "--refresh" in words, words

    async def test_a_write_outside_work_is_refused(self, tree: Tree) -> None:
        words = await judge(http_scope(tree), "mesh", mesh(out="../../kari/work/m.vtk"))
        assert words is not None and "out" in words, words

    async def test_a_relative_write_is_judged_against_the_callers_cwd(
        self, tree: Tree, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """``mesh.vtk`` lands in ``work``, whatever the server's own cwd is."""
        monkeypatch.chdir(tree.outside)
        assert await judge(http_scope(tree), "mesh", mesh(out="mesh.vtk")) is None

    async def test_an_explicit_fetch_cache_outside_work_is_refused(self, tree: Tree) -> None:
        words = await judge(
            http_scope(tree), "fetch", {"source": "glo30", "cache": str(tree.cache)}
        )
        assert words is not None and "cache" in words, words


# ------------------------------------------------- allowed


class TestAllowed:
    @pytest.mark.parametrize(
        "arguments",
        [
            pytest.param({"domain": "in.geojson"}, id="relative file in work"),
            pytest.param({"domain": "sub/new.geojson"}, id="relative path naming no file"),
            pytest.param({"dem": ["glo30"]}, id="bare source name"),
            pytest.param({"out_crs": "EPSG:3006"}, id="bare CRS code"),
            pytest.param({"tolerance_near": ["in.geojson", 5]}, id="tuple element in work"),
        ],
    )
    async def test_inside_work_or_naming_no_file(
        self, tree: Tree, arguments: dict[str, Any]
    ) -> None:
        assert await judge(http_scope(tree), "mesh", mesh(**arguments)) is None

    async def test_a_shared_read_root(self, tree: Tree) -> None:
        arguments = mesh(dem=[str(tree.local)], cache=str(tree.cache))
        assert await judge(http_scope(tree), "mesh", arguments) is None

    async def test_fetch_without_refresh_is_not_refused_over_http(self, tree: Tree) -> None:
        assert await judge(http_scope(tree), "fetch", {"source": "glo30", "refresh": False}) is None

    async def test_stdio_keeps_no_read_limit_and_refuses_no_flag(self, tree: Tree) -> None:
        scope = stdio_scope(tree)
        out = str(tree.root / "m.vtk")
        assert await judge(scope, "mesh", mesh(out=out, domain="/etc/hosts")) is None
        assert await judge(scope, "fetch", {"source": "glo30", "refresh": True}) is None


# ------------------------------------------------- over HTTP


@pytest.fixture(scope="module")
def roots(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A data root whose sources file names one local dataset, outside the output root."""
    root = tmp_path_factory.mktemp("read-limit").resolve()
    dataset = root / "datasets" / "local.tif"
    dataset.parent.mkdir()
    dataset.write_bytes(b"II*\x00")
    write_sources(root / "data" / "sources.toml", {"localdem": local_entry(dataset)})
    return root


@pytest.fixture(scope="module")
def server(roots: Path) -> Iterator[HttpServer]:
    with rasputin_server(roots / "server", data=roots / "data", jobs=1) as running:
        running.work("ola").mkdir(parents=True, exist_ok=True)
        yield running


class TestOverHttp:
    async def test_a_refused_call_never_starts_the_child(self, server: HttpServer) -> None:
        count = server.out_root.parent / "refused-count"
        server.set_spec(count_file=str(count))
        result = await call(server, "ola", "mesh", mesh(domain="/etc/hosts"))
        assert "domain" in error_text(result)
        assert not count.exists(), "the child ran on a refused call"

    async def test_a_local_catalogue_path_is_a_read_root(
        self, server: HttpServer, roots: Path
    ) -> None:
        count = server.out_root.parent / "catalogue-count"
        server.set_spec(count_file=str(count))
        dataset = str(roots / "datasets" / "local.tif")
        result = await call(server, "ola", "mesh", mesh(dem=[dataset]))
        assert not result.is_error, result
        assert count.read_text().count("\n") == 1

    async def test_fetch_refresh_is_refused_before_the_child(self, server: HttpServer) -> None:
        count = server.out_root.parent / "refresh-count"
        server.set_spec(count_file=str(count))
        result = await call(server, "ola", "fetch", {"source": "glo30", "refresh": True})
        assert "--refresh" in error_text(result)
        assert not count.exists()

    async def test_a_link_swapped_in_while_waiting_for_the_slot_is_refused(
        self, server: HttpServer, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The first call holds the one job slot until the test releases it; the second
        passes the first check (its heartbeat says it waits), then ``in.geojson`` becomes
        a link to ``/etc/hosts``; the check after the slot refuses it."""
        monkeypatch.setattr(mcp_server, "HEARTBEAT_SECONDS", 0.05)
        scratch = server.out_root.parent
        count, release, started = scratch / "slot-count", scratch / "release", scratch / "pid"
        server.set_spec(count_file=str(count), wait_for=str(release), pid_file=str(started))
        target = server.work("ola") / "slot-in.geojson"
        target.write_text("{}")
        waiting = anyio.Event()

        async def progress(value: float, total: float | None, message: str | None) -> None:
            if message and message.startswith("waiting for a running job"):
                waiting.set()

        results: dict[str, Any] = {}

        async def first() -> None:
            results["first"] = await call(server, "ola", "version", {})

        async def second() -> None:
            results["second"] = await call(
                server, "ola", "mesh", mesh(domain="slot-in.geojson"), progress_callback=progress
            )

        with anyio.fail_after(30):
            async with anyio.create_task_group() as group:
                group.start_soon(first)
                await until_exists(started)
                group.start_soon(second)
                await waiting.wait()
                target.unlink()
                target.symlink_to("/etc/hosts")
                release.write_text("")
        assert not results["first"].is_error, results["first"]
        assert "domain" in error_text(results["second"])
        assert count.read_text().count("\n") == 1, "the second call's child ran"
