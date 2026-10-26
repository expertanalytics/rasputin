"""`network-catchment --network`, `--lakes` and `--outlets` take a catalogue
source of features or a file, through the helper `mesh --reservoirs` uses.

`docs/increments/106-wednesday-cells.md`, section 3. The network, its outlet
points and the lake are `test_cli_network_catchment.py`'s (1 km squares in
EPSG:3006); the sources file is `$RASPUTIN_DATA/sources.toml`, written by the
test, each entry a local file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from shapely.geometry import box

from cli_helpers import squashed
from dataset_fixtures import write_sources
from gpkg_fixtures import HAS_RTREE, Layer, Row, write_gpkg
from test_cli_network_catchment import (
    ABOVE_SE100,
    AREA,
    ATTRIBUTES,
    DOWN,
    ID,
    NAME,
    NETWORK,
    at,
    lakes_geojson,
    network_geojson,
    outlets_geojson,
    refuse,
    run,
    utm,
)
from vector_fixtures import vector_entry

KEYS = {"--network": "test-units", "--lakes": "test-lakes", "--outlets": "test-outlets"}
#: The seed is 250 m from SE100's outlet, the lake under it 85 m (the lake test's).
SEED = utm(*at(500, 250))


@pytest.fixture
def files(tmp_path: Path) -> dict[str, Path]:
    """The three files, each also registered under its key in `KEYS`."""
    own = tmp_path / "own"
    own.mkdir()
    return {
        "--network": network_geojson(own / "units.geojson"),
        "--lakes": lakes_geojson(own / "lakes.geojson", box(*at(400, 85), *at(600, 300))),
        "--outlets": outlets_geojson(own / "outlets.geojson"),
    }


def install(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, entries: dict[str, dict[str, Any]]
) -> None:
    data = tmp_path / "data"
    write_sources(data / "sources.toml", entries)
    monkeypatch.setenv("RASPUTIN_DATA", str(data))


@pytest.mark.parametrize("option", list(KEYS))
def test_a_key_writes_the_bytes_its_file_does(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, files: dict[str, Path], option: str
) -> None:
    install(
        monkeypatch,
        tmp_path,
        {KEYS[o]: vector_entry(p, f"Test {o}") for o, p in files.items()},
    )
    by_path = [a for o, p in files.items() for a in (o, str(p))]
    by_key = [a for o, p in files.items() for a in (o, KEYS[o] if o == option else str(p))]
    props, _ = run(tmp_path, *by_path, *ATTRIBUTES, *SEED, out="path.geojson")
    assert props["units"] == ABOVE_SE100
    run(tmp_path, *by_key, *ATTRIBUTES, *SEED, out="key.geojson")
    assert (tmp_path / "key.geojson").read_bytes() == (tmp_path / "path.geojson").read_bytes()


def two_tables(path: Path) -> Path:
    """The network as table `units`, beside a table `decoy` whose ids are all
    other, so `--unit SE100` is found only in `units`."""
    columns = (ID, DOWN, AREA, NAME)
    rows = [
        Row(k, g, {ID: o, DOWN: d, AREA: a, NAME: n})
        for k, (o, d, a, n, g) in enumerate(NETWORK, 1)
    ]
    decoy = [Row(r.pk, r.geometry, {**r.values, ID: "X" + str(r.values[ID])}) for r in rows]
    return write_gpkg(
        path,
        [Layer("decoy", 3006, decoy, columns=columns), Layer("units", 3006, rows, columns=columns)],
    )


@pytest.mark.skipif(not HAS_RTREE, reason="this SQLite has no R-tree module")
@pytest.mark.parametrize(
    ("entry_layer", "given"), [("units", ()), ("decoy", ("--network-layer", "units"))]
)
def test_the_entrys_layer_unless_one_is_given(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, entry_layer: str, given: tuple[str, ...]
) -> None:
    gpkg = two_tables(tmp_path / "two.gpkg")
    install(monkeypatch, tmp_path, {"test-units": vector_entry(gpkg, "Units", layer=entry_layer)})
    props, _ = run(tmp_path, "--network", "test-units", *given, *ATTRIBUTES, "--unit", "SE100")
    assert props["units"] == ABOVE_SE100


@pytest.mark.parametrize("option", list(KEYS))
def test_neither_a_source_nor_a_file_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, files: dict[str, Path], option: str
) -> None:
    install(monkeypatch, tmp_path, {})
    args = [a for o, p in files.items() for a in (o, "no-such-thing" if o == option else str(p))]
    output = refuse(tmp_path, *args, *ATTRIBUTES, *SEED, says=(f"Invalid value for {option}",))
    wording = "no-such-thing is neither a catalogue source of features nor a .gpkg, .geojson, "
    assert squashed(wording + ".json file") in squashed(output), output


def test_an_uncached_source_names_the_fetch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    entry = {
        "kind": "arcgis-features", "url": "https://maps.example.org/arcgis/rest/services/U/MapServer/1",
        "crs": "EPSG:3006", "fields": [ID, DOWN],
        "credit": "Test units, {url} (accessed on {date})", "licence_note": "Test licence",
    }  # fmt: skip
    install(monkeypatch, tmp_path, {"test-units": entry})
    output = refuse(
        tmp_path, "--network", "test-units", "--cache", str(tmp_path / "cache"), *ATTRIBUTES,
        "--unit", "SE100", says=("Invalid value for --network",),
    )  # fmt: skip
    flat = squashed(output)
    assert "notinthecache" in flat and "rasputinfetchtest-units" in flat, output


@pytest.mark.parametrize("option", list(KEYS))
def test_a_key_that_is_also_a_file_here_is_refused_naming_both(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, files: dict[str, Path], option: str
) -> None:
    """A catalogue id with a vector suffix, and a file of that name in the
    working directory: neither is taken silently; a usage error names both."""
    key = f"{KEYS[option]}.geojson"
    quoted = f'"{key}"'  # a dotted TOML key needs quotes
    install(monkeypatch, tmp_path, {quoted: vector_entry(files[option], f"Test {option}")})
    here = tmp_path / "here"
    here.mkdir()
    (here / key).write_bytes(files[option].read_bytes())
    monkeypatch.chdir(here)
    args = [a for o, p in files.items() for a in (o, key if o == option else str(p))]
    output = refuse(tmp_path, *args, *ATTRIBUTES, *SEED, says=(f"Invalid value for {option}",))
    flat = squashed(output)
    assert squashed(key) in flat and "catalogue" in flat and "file" in flat, output
    assert "neither" not in flat, output
