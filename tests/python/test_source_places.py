"""A local source in the catalogue, and where the user's sources file is found.

`docs/increments/55-dataset-discovery.md`, sections 4.2 and 4.3, tests T1
(catalogue rules) and T2 (places). Pure: `RemoteSource` and `catalogue` on
files under `tmp_path`, and the two place functions on environment mappings,
so no test here reads the process environment.

The new names (`user_source_candidates`, `user_source_file`) are read off the
module at test time, so their absence fails each test, not collection.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from dataset_fixtures import SHIPPED, cached_entry, local_entry, write_sources

ID = "my-local"


@pytest.fixture(scope="module")
def sources() -> ModuleType:
    return importlib.import_module("tin_engine.sources")


def entry(sources: ModuleType, table: dict[str, Any]) -> Any:
    return sources.RemoteSource.model_validate({**table, "id": ID})


def refused(sources: ModuleType, table: dict[str, Any], *says: str) -> str:
    """`table` refused with a `ValueError` whose message carries each of `says`."""
    with pytest.raises(ValueError) as caught:
        entry(sources, table)
    message = str(caught.value)
    for word in says:
        assert word in message, (word, message)
    return message


# ------------------------------------------------------------ T1


class TestT1LocalEntries:
    def test_a_local_entry_with_a_path_loads(self, sources: ModuleType, tmp_path: Path) -> None:
        got = entry(sources, local_entry(tmp_path / "tiles"))
        assert got.kind == "local"
        assert got.path == tmp_path / "tiles"
        assert got.crs is None
        assert got.read_unit == "metre"

    def test_a_relative_path_is_taken_from_the_files_directory(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        user = write_sources(
            tmp_path / "conf" / "sources.toml", {ID: local_entry(Path("data/tiles"))}
        )
        got = sources.catalogue(user).sources[ID]
        assert Path(got.path).resolve() == (tmp_path / "conf" / "data" / "tiles").resolve()

    def test_an_absolute_path_is_kept(self, sources: ModuleType, tmp_path: Path) -> None:
        user = write_sources(tmp_path / "conf" / "sources.toml", {ID: local_entry(tmp_path / "t")})
        assert sources.catalogue(user).sources[ID].path == tmp_path / "t"

    def test_a_missing_path_is_not_the_catalogues_problem(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        """The listing reports it, and a run refuses it as any unreadable --dem."""
        user = write_sources(tmp_path / "sources.toml", {ID: local_entry(tmp_path / "nowhere")})
        assert not (tmp_path / "nowhere").exists()
        assert sources.catalogue(user).sources[ID].kind == "local"

    def test_without_a_path_refused(self, sources: ModuleType) -> None:
        refused(sources, local_entry(None), "path")

    @pytest.mark.parametrize(
        ("field", "value"),
        [
            ("url", "https://example.invalid/a.tif"),
            ("url_template", "https://example.invalid/{name}.tif"),
            ("tile_list_url", "https://example.invalid/list.txt"),
            ("tile_name", "T_{ns}{lat}_{ew}{lon}"),
            ("tile_spacing", (0.1, 0.1)),
            ("years", (2020, 2021)),
            ("crs", "EPSG:25833"),
        ],
    )
    def test_each_remote_field_is_refused_naming_it(
        self, sources: ModuleType, tmp_path: Path, field: str, value: Any
    ) -> None:
        refused(sources, local_entry(tmp_path, **{field: value}), field)

    def test_a_local_class_raster_is_refused_pointing_at_land_cover(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        table = local_entry(tmp_path, values="classes", z_unit=None)
        refused(sources, table, "values", "--land-cover")

    def test_the_refusal_names_the_file_and_the_entry(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        user = write_sources(
            tmp_path / "my-sources.toml", {ID: local_entry(tmp_path, crs="EPSG:25833")}
        )
        with pytest.raises(ValueError) as caught:
            sources.catalogue(user)
        message = str(caught.value)
        for word in ("my-sources.toml", f"source.{ID}", "crs"):
            assert word in message, (word, message)


class TestT1RemoteEntries:
    def test_a_remote_entry_without_crs_is_refused(self, sources: ModuleType) -> None:
        refused(sources, cached_entry(crs=None), "crs")

    def test_a_remote_entry_with_a_path_is_refused(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        refused(sources, cached_entry(path=str(tmp_path)), "path")

    def test_control_the_remote_entry_loads(self, sources: ModuleType) -> None:
        got = entry(sources, cached_entry())
        assert (got.kind, got.crs, got.path) == ("one-cog", "EPSG:25833", None)


class TestT1Description:
    def test_description_defaults_to_empty(self, sources: ModuleType, tmp_path: Path) -> None:
        assert entry(sources, local_entry(tmp_path, description=None)).description == ""
        assert entry(sources, cached_entry(description=None)).description == ""

    def test_description_is_kept(self, sources: ModuleType, tmp_path: Path) -> None:
        got = entry(sources, local_entry(tmp_path, description="Bare ground, 10 m."))
        assert got.description == "Bare ground, 10 m."

    @pytest.mark.parametrize("key", SHIPPED)
    def test_each_shipped_source_says_what_it_holds(self, sources: ModuleType, key: str) -> None:
        """The agent chooses by the description; a shipped source without one
        gives it nothing to choose by. Its wording is the review's to check."""
        assert sources.SOURCES[key].description.strip(), key


# ------------------------------------------------------------ T2

H, X, R = "/home/ola", "/xdg/config", "/data/root"


class TestT2Candidates:
    @pytest.mark.parametrize(
        ("environ", "expected"),
        [
            pytest.param(
                {"RASPUTIN_DATA": R, "XDG_CONFIG_HOME": X, "HOME": H},
                (f"{R}/sources.toml", f"{X}/rasputin/sources.toml"),
                id="data root, then XDG; HOME not looked at",
            ),
            pytest.param(
                {"RASPUTIN_DATA": R, "HOME": H},
                (f"{R}/sources.toml", f"{H}/.config/rasputin/sources.toml"),
                id="data root, then HOME",
            ),
            pytest.param({"HOME": H}, (f"{H}/.config/rasputin/sources.toml",), id="HOME alone"),
            pytest.param({"XDG_CONFIG_HOME": X}, (f"{X}/rasputin/sources.toml",), id="XDG alone"),
            pytest.param(
                {"XDG_CONFIG_HOME": "rel/config", "HOME": H},
                (f"{H}/.config/rasputin/sources.toml",),
                id="relative XDG ignored",
            ),
            pytest.param(
                {"XDG_CONFIG_HOME": "", "HOME": H},
                (f"{H}/.config/rasputin/sources.toml",),
                id="empty XDG ignored",
            ),
            pytest.param({"HOME": "rel/home"}, (), id="relative HOME ignored"),
            pytest.param(
                {"RASPUTIN_DATA": "", "HOME": H},
                (f"{H}/.config/rasputin/sources.toml",),
                id="empty data root ignored",
            ),
            pytest.param({}, (), id="nothing set"),
        ],
    )
    def test_the_order(
        self, sources: ModuleType, environ: dict[str, str], expected: tuple[str, ...]
    ) -> None:
        got = sources.user_source_candidates(environ)
        assert isinstance(got, tuple)
        assert got == tuple(Path(p) for p in expected)


class TestT2TheFileRead:
    @pytest.fixture
    def places(self, tmp_path: Path) -> dict[str, Path]:
        return {"data": tmp_path / "data", "xdg": tmp_path / "xdg", "home": tmp_path / "home"}

    def environ(self, places: dict[str, Path]) -> dict[str, str]:
        return {
            "RASPUTIN_DATA": str(places["data"]),
            "XDG_CONFIG_HOME": str(places["xdg"]),
            "HOME": str(places["home"]),
        }

    def test_the_data_roots_file_when_it_exists(
        self, sources: ModuleType, places: dict[str, Path]
    ) -> None:
        first = write_sources(places["data"] / "sources.toml", {})
        write_sources(places["xdg"] / "rasputin" / "sources.toml", {})
        assert sources.user_source_file(self.environ(places)) == first

    def test_a_data_root_without_the_file_falls_through_to_xdg(
        self, sources: ModuleType, places: dict[str, Path]
    ) -> None:
        """'Set' is not 'exists': a shell with RASPUTIN_DATA for the cache and
        an agent's server without it see the same file."""
        places["data"].mkdir()
        second = write_sources(places["xdg"] / "rasputin" / "sources.toml", {})
        assert sources.user_source_file(self.environ(places)) == second

    def test_the_home_config_when_xdg_is_unset(
        self, sources: ModuleType, places: dict[str, Path]
    ) -> None:
        home = write_sources(places["home"] / ".config" / "rasputin" / "sources.toml", {})
        assert sources.user_source_file({"HOME": str(places["home"])}) == home

    def test_a_directory_named_sources_toml_is_not_a_file(
        self, sources: ModuleType, places: dict[str, Path]
    ) -> None:
        (places["data"] / "sources.toml").mkdir(parents=True)
        second = write_sources(places["xdg"] / "rasputin" / "sources.toml", {})
        assert sources.user_source_file(self.environ(places)) == second

    def test_none_when_no_candidate_exists(
        self, sources: ModuleType, places: dict[str, Path]
    ) -> None:
        assert sources.user_source_file(self.environ(places)) is None
        assert sources.user_source_file({}) is None
