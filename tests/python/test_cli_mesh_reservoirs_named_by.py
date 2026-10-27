"""`mesh --reservoirs <entry>` where the entry keeps only the lakes a dam
names, and every other option refusing such an entry.

`docs/increments/106-wednesday-cells.md`, section 5, tests 3 and 4 and the
refusals. The scene is `test_cli_mesh_reservoirs.py`'s (the polygon route with
`corine`). Two lake squares lie wholly inside its water polygon and its
domain: lake A (`SJOID` "1-1") is named by a dam (`vy_eucd` "SE1-1"); lake B
("2-2") is named only by a dam that writes "2-2", without the "SE" prefix
SMHI's dams carry, so a join that drops the prefix would keep it.

THE ORACLE: the same run with lake A alone given as a file. The reservoir
column must be the same, element for element; the run with both lakes as a
file must differ from it, so lake B is seen to add reservoir when kept.
Lake A lies wholly in water, so every triangle's reservoir share is its share
under A, and their area-weighted sum is A's area (bound 1e-6 relative: the
shares are float32, 6e-8).

PINNED HERE, where the design leaves it open (listed for `@architect`):

- The dams are read whole, not with the run's box: a dam outside the
  domain's box (as at a lake's outlet on the domain's edge) still names its lake.
- A refusal of `named_by.source` (not in the catalogue, or not a local
  `.geojson` of features) is a usage error before any phase, naming
  `--reservoirs` and the source's id.
- `mesh --catchments` refuses such an entry with the same wording as
  `network-catchment`'s three options.
"""

from __future__ import annotations

import os
import re
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import shapely
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import Point, box

from cli_helpers import USAGE, squashed
from dataset_fixtures import write_sources
from feature_fixtures import Feat, write_geojson
from fraction_fixtures import MESH_CRS, dem_origin
from test_cli_land_cover_polygons import cover_polygons, write_cover_geojson
from test_cli_mesh_cells import column, corners, mesh, scene, written
from test_cli_mesh_reservoirs import phases, record
from test_cli_network_catchment import ATTRIBUTES, refuse
from test_cli_network_catchment_sources import KEYS, SEED, files, install
from vector_fixtures import vector_entry

__all__ = ["files", "phases", "scene"]  # fixtures, used by name

KEY = "test-reservoirs"
DAMS = "test-dams"
NAMED_BY = {"source": DAMS, "field": "vy_eucd", "key": "SJOID", "prefix": "SE"}
REFUSAL = f"{KEY} keeps only the features another source names; only --reservoirs reads it"
#: On an area, relative: see the module docstring.
AREA = 1e-6


def lake_a() -> shapely.Polygon:
    x0, y0 = dem_origin()
    return box(x0 + 1500.2, y0 - 1600.3, x0 + 1900.7, y0 - 1200.6)


def lake_b() -> shapely.Polygon:
    x0, y0 = dem_origin()
    return box(x0 + 1500.4, y0 - 1100.2, x0 + 1900.9, y0 - 700.8)


def write_lakes(path: Path, *which: str) -> Path:
    shapes = {"1-1": lake_a(), "2-2": lake_b()}
    feats = [Feat(k + 1, shapes[s], {"SJOID": s}) for k, s in enumerate(which)]
    return write_geojson(path, feats, crs=MESH_CRS)


def write_dams(path: Path, outside: bool) -> Path:
    """A's dam at its south edge, or 5 km beyond the domain's box; B's dam
    without the prefix; one with no name, one naming a Norwegian lake."""
    x0, y0 = dem_origin()
    a = Point(x0 + 7000.3, y0 - 7000.1) if outside else Point(x0 + 1700.5, y0 - 1600.3)
    feats = [
        Feat(1, a, {"vy_eucd": "SE1-1"}),
        Feat(2, Point(x0 + 1700.6, y0 - 700.8), {"vy_eucd": "2-2"}),
        Feat(3, Point(x0 + 500.1, y0 - 500.2), {"vy_eucd": None}),
        Feat(4, Point(x0 + 600.1, y0 - 600.2), {"vy_eucd": "NO9"}),
    ]
    return write_geojson(path, feats, crs=MESH_CRS)


def install_named(monkeypatch: pytest.MonkeyPatch, entries: dict[str, dict[str, Any]]) -> None:
    """`entries` as the user's sources file (`test_cli_mesh_reservoirs`' place)."""
    monkeypatch.delenv("RASPUTIN_DATA", raising=False)
    config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
    write_sources(config, entries)


def named_entry(lakes: Path, **named_by: str) -> dict[str, Any]:
    return vector_entry(lakes, "Test reservoirs", named_by={**NAMED_BY, **named_by})


@pytest.fixture
def cover(tmp_path: Path) -> Path:
    return write_cover_geojson(tmp_path / "cover.geojson", cover_polygons())


def args(cover: Path) -> tuple[str, ...]:
    """The land cover and its legend, then `--reservoirs` awaiting its value."""
    return ("--land-cover", str(cover), "--land-cover-rules", "corine", "--reservoirs")


def area_weighted(conn: sqlite3.Connection, name: str) -> float:
    triangles = shapely.area(shapely.polygons(corners(conn)[:, :, :2]))
    return float((column(conn, "cells", name) * triangles).sum())


@pytest.mark.parametrize("outside", [False, True], ids=["dam-inside", "dam-outside-the-box"])
def test_only_the_lake_a_dam_names_is_a_reservoir(
    scene: dict[str, Path],
    tmp_path: Path,
    cover: Path,
    monkeypatch: pytest.MonkeyPatch,
    outside: bool,
) -> None:
    own = tmp_path / "own"
    own.mkdir()
    lakes = write_lakes(own / "lakes.geojson", "1-1", "2-2")
    dams = write_dams(own / "dams.geojson", outside)
    install_named(monkeypatch, {KEY: named_entry(lakes), DAMS: vector_entry(dams, "Test dams")})
    only_a = write_lakes(tmp_path / "a.geojson", "1-1")
    with (
        closing(written(scene, tmp_path / "key.gpkg", *args(cover), KEY)) as got,
        closing(written(scene, tmp_path / "a.gpkg", *args(cover), str(only_a))) as want,
        closing(written(scene, tmp_path / "both.gpkg", *args(cover), str(lakes))) as both,
    ):
        reservoir = column(got, "cells", "reservoir")
        assert_array_equal(reservoir, column(want, "cells", "reservoir"))
        assert_array_equal(column(got, "cells", "lake"), column(want, "cells", "lake"))
        assert not np.array_equal(reservoir, column(both, "cells", "reservoir")), (
            "lake B adds no reservoir even when kept: the oracle cannot fail"
        )
        assert_allclose(area_weighted(got, "reservoir"), lake_a().area, rtol=AREA)
        value = record(got, "reservoirs")
    assert value is not None and not re.search(r"\bnone\b", value), value
    assert KEY in value, value


# ------------------------------------------------------------ refusals


@pytest.mark.parametrize("missing", ["absent", "gpkg"], ids=["not-in-catalogue", "a-gpkg"])
def test_a_named_by_source_that_is_not_a_local_geojson_is_refused_naming_it(
    scene: dict[str, Path],
    tmp_path: Path,
    cover: Path,
    monkeypatch: pytest.MonkeyPatch,
    phases: list[str],
    missing: str,
) -> None:
    own = tmp_path / "own"
    own.mkdir()
    lakes = write_lakes(own / "lakes.geojson", "1-1", "2-2")
    entries = {KEY: named_entry(lakes)}
    if missing == "gpkg":
        dams = own / "dams.gpkg"
        dams.write_bytes(b"")  # refused by its suffix, never read
        entries[DAMS] = vector_entry(dams, "Test dams")
    install_named(monkeypatch, entries)
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, *args(cover), KEY)
    assert code == USAGE, output
    flat = squashed(output)
    assert "--reservoirs" in flat and DAMS in flat, output
    assert phases == [], phases
    assert not target.exists()


@pytest.mark.parametrize("option", list(KEYS))
def test_network_catchment_refuses_an_entry_with_named_by(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, files: dict[str, Path], option: str
) -> None:
    install(
        monkeypatch,
        tmp_path,
        {KEY: named_entry(files[option]), DAMS: vector_entry(files["--outlets"], "Test dams")},
    )
    given = [a for o, p in files.items() for a in (o, KEY if o == option else str(p))]
    output = refuse(tmp_path, *given, *ATTRIBUTES, *SEED, says=(f"Invalid value for {option}",))
    assert squashed(REFUSAL) in squashed(output), output


def test_mesh_catchments_refuses_an_entry_with_named_by(
    scene: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, phases: list[str]
) -> None:
    own = tmp_path / "own"
    own.mkdir()
    lakes = write_lakes(own / "lakes.geojson", "1-1", "2-2")
    dams = write_dams(own / "dams.geojson", outside=False)
    install_named(monkeypatch, {KEY: named_entry(lakes), DAMS: vector_entry(dams, "Test dams")})
    target = tmp_path / "x.gpkg"
    code, output = mesh(scene, target, "--catchments", KEY)
    assert code == USAGE, output
    flat = squashed(output)
    assert "Invalidvaluefor--catchments" in flat and squashed(REFUSAL) in flat, output
    assert phases == [], phases
    assert not target.exists()
