"""A local features entry may keep only the features another source names:
`named_by = { source, field, key, prefix }` on a sources-file entry.

`docs/increments/106-wednesday-cells.md`, section 5, tests 1 and 2. A feature
is kept when `prefix + str(its key)` equals the `field` of some feature of
`source`; the four fields are required text, with no other match rule, and
the entry is a local `.geojson` source of features.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The model's field is `RemoteSource.named_by`, a frozen model whose
  attributes are `source`, `field`, `key` and `prefix`; an unknown key in it
  is refused, as on the entry.
- Each refusal's message names `named_by`.
- The pure filter is `tin_engine.sources.kept_named(keys, names, prefix)`:
  the indices, ascending, of the keys some name equals with the prefix
  before the key's `str`; a name that is None matches nothing.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from tin_engine.sources import RemoteSource

NAMED_BY = {"source": "test-dams", "field": "vy_eucd", "key": "SJOID", "prefix": "SE"}


def entry(path: str, /, **changes: Any) -> dict[str, Any]:
    """A local entry for `path` with `named_by`; a change to None drops that key."""
    table: dict[str, Any] = {
        "id": "test-reservoirs",
        "kind": "local",
        "path": path,
        "named_by": dict(NAMED_BY),
        "credit": "Test lakes",
        "licence_note": "Test licence",
    }
    table.update(changes)
    return {k: v for k, v in table.items() if v is not None}


def test_named_by_is_accepted_on_a_local_geojson_entry(tmp_path: Path) -> None:
    source = RemoteSource.model_validate(entry(str(tmp_path / "lakes.geojson")))
    assert source.values == "features"
    assert source.named_by is not None
    got = {name: getattr(source.named_by, name) for name in NAMED_BY}
    assert got == NAMED_BY


def test_an_entry_without_named_by_has_none(tmp_path: Path) -> None:
    source = RemoteSource.model_validate(entry(str(tmp_path / "lakes.geojson"), named_by=None))
    assert source.named_by is None


def layer() -> dict[str, Any]:
    return {
        "kind": "arcgis-features", "url": "https://maps.example.org/arcgis/rest/services/L/MapServer/1",
        "crs": "EPSG:3006", "fields": ["SJOID"], "path": None,
        "credit": "Test lakes, {url} (accessed on {date})",
    }  # fmt: skip


def no_prefix() -> dict[str, Any]:
    return {"named_by": {k: v for k, v in NAMED_BY.items() if k != "prefix"}}


def other_rule() -> dict[str, Any]:
    return {"named_by": {**NAMED_BY, "match": "contains"}}


@pytest.mark.parametrize(
    ("path", "changes"),
    [
        ("dem.tif", {}),
        ("dem.tif", layer()),
        ("lakes.gpkg", {}),
        ("lakes.gpkg", {"layer": "lakes"}),
        ("lakes.geojson", no_prefix()),
        ("lakes.geojson", other_rule()),
    ],
    ids=["heights", "arcgis-features", "gpkg", "gpkg-layer", "prefix-left-out", "other-key"],
)
def test_named_by_is_refused_off_a_local_geojson_entry_or_incomplete(
    tmp_path: Path, path: str, changes: dict[str, Any]
) -> None:
    with pytest.raises(ValidationError) as caught:
        RemoteSource.model_validate(entry(str(tmp_path / path), **changes))
    assert "named_by" in str(caught.value), caught.value


# ------------------------------------------------------------ the filter


def kept_named(keys: list[Any], names: list[Any], prefix: str) -> list[int]:
    found: list[int] = importlib.import_module("tin_engine.sources").kept_named(keys, names, prefix)
    return list(found)


def test_the_filter_keeps_the_keys_named_with_the_prefix() -> None:
    """Keys 1-1, 2-2, 3-3; names SE1-1, SE3-3, none, NO9; prefix SE."""
    assert kept_named(["1-1", "2-2", "3-3"], ["SE1-1", "SE3-3", None, "NO9"], "SE") == [0, 2]


def test_a_name_without_the_prefix_names_nothing() -> None:
    """SMHI's dams write `SE` before the lake's SJOID: the bare key is no match."""
    assert kept_named(["1-1", "2-2"], ["1-1", "2-2"], "SE") == []


def test_a_key_that_is_a_number_is_compared_as_its_str() -> None:
    assert kept_named([7, 8], ["SE8"], "SE") == [1]


def test_no_names_keep_nothing() -> None:
    assert kept_named(["1-1"], [], "SE") == []
