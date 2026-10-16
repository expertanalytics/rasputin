"""`rasputin sources`, and a source name as `--dem` of every command, on the CLI.

`docs/increments/55-dataset-discovery.md`, sections 4.4 and 4.5, tests T6 (a
cached key in `catchment` and `station-catchments`), T7 (the listing's own
refusals), and T5's `fetch` refusal of a local key, with sockets closed.

Every test here finds the user's sources file the way an agent's server does
with no `RASPUTIN_DATA`: in `$XDG_CONFIG_HOME/rasputin/sources.toml`, where
`conftest.py` gives each test its own empty `XDG_CONFIG_HOME`.
"""

from __future__ import annotations

import json
import os
import re
import socket
from pathlib import Path
from typing import Any

import pytest

import batch_fixtures as bf
from catchment_fixtures import EPSG
from cli_helpers import USAGE, invoke, plain, runner, squashed
from cog_fixtures import write_cache
from dataset_fixtures import (
    CACHED,
    LOCAL,
    NODES_NO_LAKE,
    SEED,
    SEED_BLOCK_ROW,
    SHIPPED,
    block_row,
    by_id,
    cached_bytes,
    cached_entry,
    local_entry,
    write_local,
    write_sources,
)
from test_cli_mesh_mosaic import write_tiles
from tin_engine.cli import app

STATIONS = "bowl-stations"


@pytest.fixture
def config(monkeypatch: pytest.MonkeyPatch) -> Path:
    """The user's sources file, found as an agent's server finds it: no `RASPUTIN_DATA`."""
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    return Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"


@pytest.fixture
def local(tmp_path: Path) -> Path:
    return write_local(tmp_path / "tiles")


def listing(*args: str) -> dict[str, Any]:
    """`rasputin sources ARGS`: exit 0 and one JSON document on one stdout line."""
    result = runner.invoke(app, ["sources", *args])
    assert result.exit_code == 0, result.output
    lines = result.stdout.splitlines()
    assert len(lines) == 1, result.stdout
    document = json.loads(lines[0])
    assert isinstance(document, dict), document
    return document


def run(command: str, *args: str) -> tuple[int, str]:
    """As `cli_helpers.invoke`, with an exception that escaped the command
    named in the output: a key read as a file escapes as `FileNotFoundError`."""
    result = runner.invoke(app, [command, *args])
    escaped = "" if result.exception is None else f" [escaped: {result.exception!r}]"
    return result.exit_code, plain(result.output) + escaped


def no_network(*args: Any, **kwargs: Any) -> Any:
    raise AssertionError("the command opened a socket")


def seed_args() -> tuple[str, ...]:
    return ("--seed", repr(SEED[0]), repr(SEED[1]), "--seed-crs", EPSG)


def nodes_of(out: Path) -> int:
    (feature,) = json.loads(out.read_text())["features"]
    return int(feature["properties"]["nodes"])


# ------------------------------------------------------------ T6


class TestT6CachedKeyInCatchment:
    @pytest.fixture
    def cache(self, tmp_path: Path, config: Path) -> Path:
        write_sources(config, {CACHED: cached_entry()})
        return tmp_path / "cache"

    def catchment(self, cache: Path, out: Path) -> tuple[int, str]:
        return invoke(
            "catchment", "--dem", CACHED, "--cache", str(cache), *seed_args(), "--out", str(out)
        )

    def test_the_full_cache(self, cache: Path, tmp_path: Path) -> None:
        write_cache(cache, CACHED, {"dem": cached_bytes()})
        out = tmp_path / "c.geojson"
        code, output = self.catchment(cache, out)
        assert code == 0, output
        assert nodes_of(out) == NODES_NO_LAKE

    def test_a_block_row_the_catchment_never_reaches_need_not_be_cached(
        self, cache: Path, tmp_path: Path
    ) -> None:
        """The windowed read: without the northern block row (rows 0-15, far
        from the catchment) the whole-tile read refused with "13 of the 325
        blocks needed are not in the cache"; the windows never meet it."""
        write_cache(cache, CACHED, {"dem": cached_bytes()}, skip={"dem": block_row(0)})
        out = tmp_path / "c.geojson"
        code, output = self.catchment(cache, out)
        assert code == 0, output
        assert nodes_of(out) == NODES_NO_LAKE

    def test_missing_blocks_under_the_seed_say_what_to_fetch(
        self, cache: Path, tmp_path: Path
    ) -> None:
        write_cache(cache, CACHED, {"dem": cached_bytes()}, skip={"dem": block_row(SEED_BLOCK_ROW)})
        out = tmp_path / "c.geojson"
        code, output = self.catchment(cache, out)
        assert code == USAGE, output
        assert not out.exists()
        assert re.search(r"\d+ of the \d+ blocks needed are not in the cache", output), output
        assert "in the source's CRS" in output, output
        hint = re.search(
            rf"run: rasputin fetch {CACHED} --bbox (\S+) (\S+) (\S+) (\S+?)(?: --cache .*)?\.?$",
            output.strip(),
        )
        assert hint, output
        x0, y0, x1, y1 = map(float, hint.groups())
        assert x0 <= SEED[0] <= x1 and y0 <= SEED[1] <= y1, (hint.groups(), SEED)

    def test_the_fetch_hint_names_the_cache_it_was_given(self, cache: Path, tmp_path: Path) -> None:
        """An agent that follows the hint must fetch into the cache the run
        read, or the rerun fails again: the hint carries `--cache PATH`, last."""
        write_cache(cache, CACHED, {"dem": cached_bytes()}, skip={"dem": block_row(SEED_BLOCK_ROW)})
        code, output = self.catchment(cache, tmp_path / "c.geojson")
        assert code == USAGE, output
        hint = squashed(output).partition("run:rasputinfetch")[2]
        assert hint, output
        assert hint.removesuffix(".").endswith(squashed(f"--cache {cache}")), output

    def test_without_cache_the_fetch_hint_names_none(self, config: Path, tmp_path: Path) -> None:
        """No `--cache` given: the run read the default cache, which a fetch
        in the same environment also writes, so the hint adds no `--cache`."""
        write_sources(config, {CACHED: cached_entry()})
        default = Path(os.environ["XDG_CACHE_HOME"]) / "rasputin"
        write_cache(
            default, CACHED, {"dem": cached_bytes()}, skip={"dem": block_row(SEED_BLOCK_ROW)}
        )
        code, output = invoke(
            "catchment", "--dem", CACHED, *seed_args(), "--out", str(tmp_path / "c.geojson")
        )
        assert code == USAGE, output
        hint = squashed(output).partition("run:rasputinfetch")[2]
        assert hint.startswith(f"{CACHED}--bbox"), output
        assert "--cache" not in hint, output

    def test_station_catchments_resolves_the_cached_key(self, cache: Path, tmp_path: Path) -> None:
        """The key reaches the cache: whatever the run says of these stations
        (they lie off the bowl), it is not a file named after the key."""
        write_cache(cache, CACHED, {"dem": cached_bytes()})
        _code, output = run(
            "station-catchments", "--dem", CACHED, "--cache", str(cache),
            "--stations", str(bf.write_stations(tmp_path / "s.geojson", [bf.TREFF])),
            "--rivers", str(bf.write_rivers(tmp_path / "r.geojson")),
            "--out-dir", str(tmp_path / "out"),
        )  # fmt: skip
        for wrong in ("cannot read", "No such option", "No such file", "escaped"):
            assert wrong not in output, output


class TestT6LocalKeyInStationCatchments:
    def test_one_station_from_a_local_key(
        self, tmp_path: Path, config: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        dem = tmp_path / "basin"
        write_tiles(dem, bf.basin_tiles())
        write_sources(config, {STATIONS: local_entry(dem)})
        monkeypatch.chdir(tmp_path)  # so a file named after the key would be looked for here
        out = tmp_path / "out"
        code, output = run(
            "station-catchments", "--dem", STATIONS,
            "--stations", str(bf.write_stations(tmp_path / "s.geojson", [bf.TREFF])),
            "--rivers", str(bf.write_rivers(tmp_path / "r.geojson")),
            "--out-dir", str(out),
        )  # fmt: skip
        assert code == 0, output
        assert (out / "results.csv").is_file(), sorted(p.name for p in out.iterdir())


# ------------------------------------------------------------ fetch of a local key


def test_fetch_refuses_a_local_key_before_any_network_use(
    config: Path, local: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    write_sources(config, {LOCAL: local_entry(local)})
    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    code, output = invoke("fetch", LOCAL, "--bbox", "15.15", "59.33", "15.20", "59.38")
    assert code != 0, output
    assert f"{LOCAL} is a local source, read in place; it needs no fetch" in output, output
    assert "Traceback" not in output, output


# ------------------------------------------------------------ T7


class TestT7ListingRefusals:
    def test_a_refused_user_file_still_lists_the_shipped_sources(
        self, config: Path, local: Path
    ) -> None:
        write_sources(config, {LOCAL: local_entry(local, colour="red")})
        document = listing()
        problem = document["user_file_problem"]
        assert isinstance(problem, str) and str(config) in plain(problem), problem
        assert "colour" in problem, problem
        found = by_id(document)
        assert set(SHIPPED) <= set(found), sorted(found)
        assert LOCAL not in found, sorted(found)

    def test_a_missing_path_is_that_sources_problem_alone(
        self, config: Path, local: Path, tmp_path: Path
    ) -> None:
        gone = tmp_path / "nowhere"
        write_sources(config, {"gone": local_entry(gone), LOCAL: local_entry(local)})
        document = listing()
        assert document["user_file_problem"] is None
        found = by_id(document)
        assert found["gone"]["present"] is False
        assert isinstance(found["gone"]["problem"], str) and found["gone"]["problem"].strip()
        assert found["gone"]["coverage"] is None
        assert (found[LOCAL]["present"], found[LOCAL]["tiles"]) == (True, 4)
        assert found[LOCAL]["problem"] is None

    def test_no_user_file_says_where_to_register_one(
        self, config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        document = listing()
        assert document["user_file"] is None
        assert document["looked_at"] == [str(config)]  # XDG set: the home config is not looked at
        assert document["register_at"] == str(config)
        note = document["note"]
        assert isinstance(note, str), document
        assert str(config) in note, note
        assert 'kind = "local"' in note, note
        assert set(SHIPPED) <= set(by_id(document))

    def test_the_user_file_read_is_named(self, config: Path, local: Path) -> None:
        write_sources(config, {LOCAL: local_entry(local)})
        document = listing()
        assert document["user_file"] == str(config)
        assert document["register_at"] == str(config)
        assert document["user_file_problem"] is None
