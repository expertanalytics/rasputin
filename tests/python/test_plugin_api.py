"""`tin_engine.plugin_api`: the names a plugin may import, each the core's own object.

`docs/increments/51-open-core.md`, *4.2 The Python API a plugin may import
(`tin_engine.plugin_api`)* and *9.1 Public (P1 unless marked)*, T7. The list
is the design's table, written out here: a name added to or dropped from
``__all__`` is a visible edit of this file. Identity, not equality: a
re-export that wrapped or copied an object would let the two drift apart.
The modules under test are imported inside each test, so each fails on its own.
"""

from __future__ import annotations

from importlib import import_module

import pytest

#: Section 4.2's table: home module, then the names re-exported from it.
PROMISED: dict[str, tuple[str, ...]] = {
    "crs": ("frame_factors", "parse_crs", "reprojector", "same_crs", "site_of"),
    "dem_input": ("DemRequest", "open_dem", "repository_for"),
    "io.models": ("Bounds", "DemTile", "RasterMeta"),
    "io.vtk_read": ("read_mesh",),
    "io.geotiff_write": ("encode_geotiff",),
    "raster": ("core_args",),
    "resources": ("Resources",),
    "run_record": ("plural",),
    "cli_kit": ("destination", "cli_resources", "words", "THREADS"),
    "catchment_core": ("accumulate",),
}

NAMES = [(home, name) for home, names in PROMISED.items() for name in names]

#: The helpers that moved from ``cli`` to ``cli_kit``, under the names ``cli``
#: still calls them by.
MOVED = {"_destination": "destination", "_resources": "cli_resources", "_words": "words"}


def test_all_is_exactly_the_promised_list() -> None:
    api = import_module("tin_engine.plugin_api")
    assert sorted(api.__all__) == sorted(name for _, name in NAMES)
    assert len(api.__all__) == len(set(api.__all__)), "a name listed twice"


@pytest.mark.parametrize(("home", "name"), NAMES, ids=[f"{h}.{n}" for h, n in NAMES])
def test_each_name_is_the_object_in_its_home_module(home: str, name: str) -> None:
    api = import_module("tin_engine.plugin_api")
    assert getattr(api, name) is getattr(import_module(f"tin_engine.{home}"), name)


@pytest.mark.parametrize(("old", "new"), sorted(MOVED.items()))
def test_cli_calls_the_moved_helpers_under_their_old_names(old: str, new: str) -> None:
    cli, kit = import_module("tin_engine.cli"), import_module("tin_engine.cli_kit")
    assert getattr(cli, old) is getattr(kit, new)


def test_cli_shares_the_threads_option_with_cli_kit() -> None:
    cli, kit = import_module("tin_engine.cli"), import_module("tin_engine.cli_kit")
    assert cli.THREADS is kit.THREADS
