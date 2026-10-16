"""Vector entries in the catalogue, and plugins that register sources.

`docs/increments/57-vector-sources.md`, sections 4.2 and 4.7, tests T1 (the
catalogue's rules for a vector entry, `~` in a sources file, a features key
refused by `--dem` and `--land-cover`) and T2 (plugin sources files: found,
listed with their origin, skipped whole on a clash or a fault, never raising).

A plugin here is an installed-looking distribution
(`plugin_fixtures.write_source_plugin`): a `.dist-info` with a
`[rasputin.sources]` entry point naming a package directory that holds
`sources.toml`. Its root is put on `sys.path` (or `PYTHONPATH` for a child),
so `importlib.metadata` finds it as it finds a real install.
"""

from __future__ import annotations

import importlib
import json
import os
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from cli_helpers import USAGE, invoke, plain, runner
from dataset_fixtures import (
    CREDIT,
    LICENCE,
    SHIPPED,
    by_id,
    lake_file,
    write_coarse,
    write_sources,
)
from plugin_fixtures import run_rasputin, write_source_plugin
from tin_engine import sources
from tin_engine.cli import app
from vector_fixtures import LAKE_DESCRIPTION, LAKES_ID, landcover_gpkg, vector_entry

#: Package names the fake plugins import; dropped from `sys.modules` after each test.
PACKAGES = ("fakeplug_data", "aplug_data", "bplug_data")
#: A package that cannot be imported.
BROKEN = 'raise ImportError("fakeplug_data is broken on purpose")\n'


def entry(path: str, **fields: Any) -> sources.RemoteSource:
    return sources.RemoteSource(
        id="v", kind="local", path=Path(path), credit=CREDIT, licence_note=LICENCE, **fields
    )


def refused_with(words: str, **fields: Any) -> str:
    """The message of the refusal of `RemoteSource(**fields)`, which must contain `words`."""
    with pytest.raises(ValueError) as caught:
        sources.RemoteSource(credit=CREDIT, licence_note=LICENCE, id="v", **fields)
    message = str(caught.value)
    assert words in message, message
    return message


# ------------------------------------------------------------ T1


class TestT1VectorEntries:
    @pytest.mark.parametrize("name", ["l.gpkg", "l.geojson", "l.json", "L.GPKG", "l.GeoJSON"])
    def test_a_vector_path_with_no_values_holds_features(self, name: str) -> None:
        assert entry(f"/data/{name}").values == "features"

    def test_a_raster_path_with_no_values_still_holds_heights(self) -> None:
        assert entry("/data/dem.tif").values == "heights"
        assert entry("/data/tiles").values == "heights"

    def test_features_need_a_local_entry(self) -> None:
        message = refused_with(
            "local",
            kind="one-cog",
            url="https://example.invalid/x.gpkg",
            crs="EPSG:25833",
            values="features",
        )
        assert "kind" in message or "values" in message, message

    @pytest.mark.parametrize("name", ["dem.tif", "lines.gml", "tiles"])
    def test_features_need_a_vector_path(self, name: str) -> None:
        refused_with(
            "values: features need a .gpkg, .geojson or .json path",
            kind="local", path=Path(f"/data/{name}"), values="features",
        )  # fmt: skip

    def test_features_take_no_height_unit(self) -> None:
        message = refused_with("features hold no heights", kind="local",
                               path=Path("/data/l.gpkg"), z_unit="metre")  # fmt: skip
        assert "z_unit" in message, message

    @pytest.mark.parametrize("values", ["heights", "classes"])
    def test_a_vector_path_cannot_hold_cells(self, values: str) -> None:
        refused_with("values: a .gpkg path holds features",
                     kind="local", path=Path("/data/l.gpkg"), values=values)  # fmt: skip

    def test_a_layer_needs_a_geopackage(self) -> None:
        assert entry("/data/l.gpkg", layer="landcover").layer == "landcover"
        assert entry("/data/l.gpkg").layer is None
        for name in ("l.geojson", "dem.tif"):
            message = refused_with("layer", kind="local", path=Path(f"/data/{name}"),
                                   layer="landcover")  # fmt: skip
            assert ".gpkg" in message, message

    def test_tilde_in_path_and_class_map_is_the_home(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        home = tmp_path / "home"
        monkeypatch.setenv("HOME", str(home))
        user = write_sources(
            tmp_path / "conf" / "sources.toml",
            {"mine": vector_entry(Path("~/x.gpkg"), LAKE_DESCRIPTION, class_map="~/m.toml")},
        )
        got = sources.catalogue(user).sources["mine"]
        assert got.path == home / "x.gpkg"
        assert Path(got.class_map) == home / "m.toml"

    def test_a_relative_path_is_still_beside_the_file(self, tmp_path: Path) -> None:
        user = write_sources(
            tmp_path / "conf" / "sources.toml",
            {"mine": vector_entry(Path("data/x.gpkg"), LAKE_DESCRIPTION, class_map="m.toml")},
        )
        got = sources.catalogue(user).sources["mine"]
        assert got.path == tmp_path / "conf" / "data" / "x.gpkg"
        assert Path(got.class_map) == tmp_path / "conf" / "m.toml"

    def test_an_absolute_path_wins_the_join(self, tmp_path: Path) -> None:
        user = write_sources(
            tmp_path / "conf" / "sources.toml",
            {"mine": vector_entry(tmp_path / "elsewhere" / "x.gpkg", LAKE_DESCRIPTION)},
        )
        assert sources.catalogue(user).sources["mine"].path == tmp_path / "elsewhere" / "x.gpkg"

    @pytest.mark.parametrize(
        ("option", "wanted"), [("--dem", "heights"), ("--land-cover", "classes")]
    )
    def test_a_features_key_is_refused_by_the_cell_options(
        self, option: str, wanted: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Refused before the path is used: the key's file does not exist."""
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        config = Path(os.environ["XDG_CONFIG_HOME"]) / "rasputin" / "sources.toml"
        write_sources(config, {LAKES_ID: vector_entry(tmp_path / "gone.geojson", LAKE_DESCRIPTION)})
        dem = str(write_coarse(tmp_path / "dem" / "bowl.tif"))
        args = ["--dem", LAKES_ID] if option == "--dem" else ["--dem", dem, option, LAKES_ID]
        code, output = invoke("mesh", *args, "--out", str(tmp_path / "m.vtk"))
        assert code == USAGE, output
        words = f"{option} takes a source of {wanted}; {LAKES_ID} holds features"
        assert words in " ".join(output.split()), output
        assert not (tmp_path / "m.vtk").exists()


# ------------------------------------------------------------ T2


@pytest.fixture
def plugins(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """A directory on `sys.path` for fake distributions; their packages are
    forgotten afterwards, so the next test imports its own."""
    root = tmp_path / "site"
    root.mkdir()
    monkeypatch.syspath_prepend(str(root))
    importlib.invalidate_caches()
    yield root
    for name in PACKAGES:
        sys.modules.pop(name, None)
    importlib.invalidate_caches()


def plugin(root: Path, entries: dict[str, dict[str, Any]], **kwargs: Any) -> Path:
    """A source plugin under `root` whose `sources.toml` holds `entries`."""
    package = kwargs.get("package", "fakeplug_data")
    write_source_plugin(root, **kwargs)
    write_sources(root / package / "sources.toml", entries)
    return root / package / "sources.toml"


def mkdir(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def lakes_table(directory: Path) -> dict[str, Any]:
    return vector_entry(lake_file(directory / "lakes.geojson"), LAKE_DESCRIPTION)


def the_plugin_file(dist: str) -> Any:
    found = [p for p in sources.plugin_source_files() if p.dist == dist]
    assert len(found) == 1, sources.plugin_source_files()
    return found[0]


class TestT2PluginFiles:
    def test_the_fake_distributions_file_is_found(self, plugins: Path, tmp_path: Path) -> None:
        toml = plugin(plugins, {"plug-lakes": lakes_table(tmp_path)})
        found = the_plugin_file("fakeplug")
        assert found.error is None, found
        assert found.path is not None and Path(found.path) == toml

    def test_its_entries_are_catalogued_with_their_origin(
        self, plugins: Path, tmp_path: Path
    ) -> None:
        plugin(plugins, {"plug-lakes": lakes_table(tmp_path)})
        got = sources.catalogue(None, sources.plugin_source_files())
        assert got.sources["plug-lakes"].values == "features"
        assert got.origins["plug-lakes"] == "plugin fakeplug"
        for key in SHIPPED:
            assert got.origins[key] == "package", key
        assert got.refused == ()

    def test_a_relative_path_is_beside_the_plugins_file(self, plugins: Path) -> None:
        toml = plugin(plugins, {"plug-lakes": vector_entry(Path("lakes.geojson"), "x")})
        got = sources.catalogue(None, sources.plugin_source_files())
        assert got.sources["plug-lakes"].path == toml.parent / "lakes.geojson"

    def test_the_listing_names_the_plugin(
        self, plugins: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        plugin(plugins, {"plug-lakes": lakes_table(tmp_path)})
        result = runner.invoke(app, ["sources"])
        assert result.exit_code == 0, result.output
        document = json.loads(result.stdout)
        found = by_id(document)
        assert found["plug-lakes"]["origin"] == "plugin fakeplug"
        assert found["glo30"]["origin"] == "package"
        assert document["plugin_problems"] == []

    @pytest.mark.parametrize("clash", ["package", "user", "earlier plugin"])
    def test_a_clash_skips_the_plugins_whole_file(
        self, clash: str, plugins: Path, tmp_path: Path
    ) -> None:
        user_file = None
        if clash == "user":
            user_file = write_sources(
                tmp_path / "conf" / "sources.toml", {"taken": lakes_table(mkdir(tmp_path / "u"))}
            )
        if clash == "earlier plugin":
            plugin(plugins, {"taken": lakes_table(mkdir(tmp_path / "a"))},
                   dist="aplug", package="aplug_data")  # fmt: skip
        taken = "glo30" if clash == "package" else "taken"
        theirs = {
            "plug-other": lakes_table(mkdir(tmp_path / "b")),
            taken: lakes_table(tmp_path / "b"),
        }
        plugin(plugins, theirs, dist="bplug", package="bplug_data")
        files = sources.plugin_source_files()
        without = sources.catalogue(user_file, tuple(p for p in files if p.dist != "bplug"))
        got = sources.catalogue(user_file, files)
        assert "plug-other" not in got.sources, sorted(got.sources)
        assert dict(got.sources) == dict(without.sources)
        assert dict(got.origins) == dict(without.origins)
        (line,) = got.refused
        assert "bplug" in line and taken in line, line

    def test_an_unreadable_plugin_file_is_one_refused_line(self, plugins: Path) -> None:
        write_source_plugin(plugins, sources='[source.bad]\nkind = "local"\ncolour = "red"\n')
        got = sources.catalogue(None, sources.plugin_source_files())
        (line,) = got.refused
        assert "fakeplug" in line, line
        assert "bad" not in got.sources

    @pytest.mark.parametrize(
        ("value", "init", "words"),
        [
            ("fakeplug_data:thing", "thing = 1\n", "fakeplug"),
            ("fakeplug_data", BROKEN, "broken on purpose"),
            ("fakeplug_data", "", "sources.toml"),
        ],
        ids=["attr", "import-error", "no-file"],
    )  # fmt: skip
    def test_a_faulty_entry_point_is_a_refused_line_never_an_exception(
        self, value: str, init: str, words: str, plugins: Path
    ) -> None:
        write_source_plugin(plugins, value=value, init=init)
        found = the_plugin_file("fakeplug")
        assert found.path is None and found.error, found
        got = sources.catalogue(None, sources.plugin_source_files())
        (line,) = got.refused
        assert "fakeplug" in line and words in line, line
        assert set(SHIPPED) <= set(got.sources)

    def test_plugin_files_come_in_distribution_name_order(self, plugins: Path) -> None:
        write_source_plugin(plugins, dist="bplug", package="bplug_data", sources="")
        write_source_plugin(plugins, dist="aplug", package="aplug_data", sources="")
        dists = [p.dist for p in sources.plugin_source_files()]
        mine = [d for d in dists if d in ("aplug", "bplug")]
        assert mine == ["aplug", "bplug"], dists

    def test_mesh_with_a_broken_plugin_says_so_on_stderr_and_runs(
        self, plugins: Path, tmp_path: Path
    ) -> None:
        write_source_plugin(plugins, init=BROKEN)
        dem = write_coarse(tmp_path / "dem" / "bowl.tif")
        out = tmp_path / "m.vtk"
        done = run_rasputin(
            "mesh", "--dem", str(dem), "--tolerance", "5", "--out", str(out), path=[plugins]
        )
        assert done.returncode == 0, done.stderr
        lines = [line for line in plain(done.stderr).splitlines() if "fakeplug" in line]
        assert lines, done.stderr
        assert all(line.startswith("rasputin: ") for line in lines), lines
        assert "broken on purpose" in " ".join(lines), lines
        assert out.is_file()


def test_the_land_cover_geopackage_entry_holds_features(tmp_path: Path) -> None:
    """A GeoPackage entry with its layer and class map, as section 4.2 writes it."""
    path = landcover_gpkg(tmp_path / "landcover.gpkg")
    got = entry(str(path), layer="landcover", class_map="corine")
    assert (got.values, got.layer, got.class_map) == ("features", "landcover", "corine")
    assert got.z_unit is None
