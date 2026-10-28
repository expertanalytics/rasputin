"""The source catalogue as data: the shipped TOML file and the user's own.

`docs/increments/44-catalogue-as-data.md`, sections 4.2 to 4.4 and section 7's
T1 (the shipped catalogue unchanged), T2 (file refusals), T3 (user entries)
and T6 (tile names). The oracle for T1 is a snapshot of every shipped entry's
`model_dump(mode="json")` and `notice()` text, taken from the Python entries
at `b1353075` (the merge base) before they moved into TOML:
`tests/python/data/catalogue_b1353075.json`. The snapshot's keys are spelled
as the model spells its fields now (a renamed field is renamed there too);
its values are the entries' values at that commit.

User files are written by this suite under `tmp_path`; T2's is named
`my-sources.toml`, so a refusal that names "the file" is told apart from one
naming the shipped `sources.toml`. No network: the one fetch is served by
`fetch_fixtures.RangeServer` on 127.0.0.1, and the mesh after it runs with
sockets closed.

PINNED HERE, where the design leaves it open (listed for `@architect`):

- `catalogue(path)` for a path that does not exist gives the shipped
  entries alone, as `catalogue(None)` does.
- A refusal's message carries the user file's name and the entry's key path
  `source.<id>`, plus the name of the offending field (`crs`, `tile_name`,
  ...); a TOML syntax error carries "line N".
- A user id equal to a shipped one is refused with a message naming
  `sources.toml` twice (both files have that name) and the id.
- A user entry's relative `rules` resolves to the path of the file
  beside the user file (compared after `Path.resolve()`).
- `sources.data_root` treats an empty `RASPUTIN_DATA` as unset, as `cache_root`
  does today.
- Every shipped entry other than `glo30` has `tile_name` and `tile_spacing`
  None; T1 compares the rest of each dump with the snapshot, less `path`
  (None on every shipped entry) and `description`, which came later.
"""

from __future__ import annotations

import importlib
import json
import socket
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import tin_engine.cli as cli
from cli_helpers import USAGE, geojson, invoke, plain, runner, squashed
from fetch_fixtures import PROJECTED_CRS, RangeServer, projected
from test_cli_fetch import TRIANGLE, no_network
from test_cli_mesh_mosaic import same_mesh
from tin_engine.cli import app
from vtkread import read_vtk

SNAPSHOT = Path(__file__).parent / "data" / "catalogue_b1353075.json"
SHIPPED = ("anadem-v1", "glo30", "mapbiomas-c11")
STATIONS = ("nve-hrd",)
#: The entries added after the snapshot: vector layers fetched from a map service.
FEATURE_LAYERS = ("nve-lakes", "nve-glaciers", "clc2018")
GLO30_TILE_NAME = "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM"
#: 10" and 1" in degrees, each rounded up in the fifth significant digit.
GLO30_TILE_SPACING = (0.0027778, 0.00027778)
TILE_FIELDS = ("tile_name", "tile_spacing")
#: Fields every entry gained with local sources; `description` is checked
#: non-empty in `test_source_places.py`, `path` is None on every remote kind.
LOCAL_FIELDS = ("path", "description")
#: The field every entry gained with vector sources: None on every shipped one.
VECTOR_FIELDS = ("layer",)
#: The field the map-service kind added: the attributes asked of the layer.
FEATURE_FIELDS = ("fields",)
#: The field the named filter added: None on every shipped entry.
NAMED_FIELDS = ("named_by",)
USER_ID = "my-x"


@pytest.fixture(scope="module")
def sources() -> ModuleType:
    return importlib.import_module("tin_engine.sources")


@pytest.fixture(scope="module")
def snapshot() -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    return loaded


# ------------------------------------------------------------ TOML text


def value(v: Any) -> str:
    """One TOML value: JSON's string escapes are TOML's basic-string escapes."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(value(x) for x in v) + "]"
    return repr(v)


def toml(tables: dict[str, dict[str, Any]]) -> str:
    """`{"source.my-x": {...}}` as TOML tables, keys in the order given."""
    out: list[str] = []
    for header, fields in tables.items():
        out.append(f"[{header}]")
        out += [f"{k} = {value(v)}" for k, v in fields.items()]
        out.append("")
    return "\n".join(out)


def write(path: Path, tables: dict[str, dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(toml(tables), encoding="utf-8")
    return path


ONE_COG: dict[str, Any] = {
    "kind": "one-cog",
    "url": "https://example.invalid/x.tif",
    "crs": "EPSG:25833",
    "credit": "Test credit",
    "licence_note": "Test licence",
}
TILES: dict[str, Any] = {
    "kind": "cog-tiles",
    "url_template": "https://example.invalid/{name}/{name}.tif",
    "tile_list_url": "https://example.invalid/tileList.txt",
    "tile_name": GLO30_TILE_NAME,
    "tile_spacing": list(GLO30_TILE_SPACING),
    "crs": "EPSG:4326",
    "credit": "Test credit, fetched {date} from {url}",
    "licence_note": "Test licence",
}
PER_YEAR: dict[str, Any] = {
    "kind": "one-cog-per-year",
    "url_template": "https://example.invalid/classes_{year}.tif",
    "years": [2020, 2021],
    "values": "classes",
    "crs": "EPSG:4326",
    "credit": "Test classes, accessed on {date} through the link: {url}",
    "licence_note": "Test licence",
}
BASES = {"one-cog": ONE_COG, "cog-tiles": TILES, "one-cog-per-year": PER_YEAR}


def changed(base: dict[str, Any], *drop: str, **update: Any) -> dict[str, Any]:
    return {k: v for k, v in base.items() if k not in drop} | update


# ------------------------------------------------------------ T1


class TestT1TheShippedCatalogue:
    def test_the_shipped_keys(self, sources: ModuleType) -> None:
        assert sorted(sources.SOURCES) == sorted([*SHIPPED, *FEATURE_LAYERS])
        assert sorted(sources.STATION_SOURCES) == sorted(STATIONS)

    @pytest.mark.parametrize("key", SHIPPED)
    def test_each_source_dumps_as_on_master_but_the_added_fields(
        self, sources: ModuleType, snapshot: dict[str, Any], key: str
    ) -> None:
        dumped = sources.SOURCES[key].model_dump(mode="json")
        added = TILE_FIELDS + LOCAL_FIELDS + VECTOR_FIELDS + FEATURE_FIELDS + NAMED_FIELDS
        rest = {k: v for k, v in dumped.items() if k not in added}
        assert rest == snapshot["sources"][key]
        assert dumped["path"] is None
        assert dumped.get("layer") is None
        assert dumped["named_by"] is None

    def test_glo30_gains_its_tile_name_and_spacing(self, sources: ModuleType) -> None:
        glo30 = sources.SOURCES["glo30"]
        assert glo30.tile_name == GLO30_TILE_NAME
        assert tuple(glo30.tile_spacing) == GLO30_TILE_SPACING
        for key in set(SHIPPED) - {"glo30"}:
            entry = sources.SOURCES[key]
            assert (entry.tile_name, entry.tile_spacing) == (None, None), key

    @pytest.mark.parametrize("key", STATIONS)
    def test_each_station_list_dumps_as_on_master(
        self, sources: ModuleType, snapshot: dict[str, Any], key: str
    ) -> None:
        dumped = sources.STATION_SOURCES[key].model_dump(mode="json")
        assert dumped == snapshot["stations"][key]

    @pytest.mark.parametrize("key", (*SHIPPED, *STATIONS))
    def test_each_notice_is_masters_byte_for_byte(
        self, sources: ModuleType, snapshot: dict[str, Any], key: str
    ) -> None:
        entry = sources.SOURCES.get(key) or sources.STATION_SOURCES[key]
        assert sources.notice(entry) == snapshot["notices"][key]

    def test_the_entries_live_in_the_packaged_toml_file(self, sources: ModuleType) -> None:
        shipped = Path(sources.__file__).parent / "data" / "sources.toml"
        text = shipped.read_text(encoding="utf-8")
        for key in SHIPPED:
            assert f"[source.{key}]" in text, key
        assert "[station.nve-hrd]" in text

    def test_read_unit_is_metre_for_classes_and_z_unit_for_heights(
        self, sources: ModuleType
    ) -> None:
        assert sources.SOURCES["mapbiomas-c11"].z_unit is None
        assert sources.SOURCES["mapbiomas-c11"].read_unit == "metre"
        for key in ("anadem-v1", "glo30"):
            entry = sources.SOURCES[key]
            assert entry.read_unit == entry.z_unit == "metre", key


# ------------------------------------------------------------ T2


def refusal(sources: ModuleType, path: Path, *says: str) -> str:
    """`catalogue(path)` refused with a `ValueError` naming the file, the
    entry's key path and each of `says`; the message."""
    with pytest.raises(ValueError) as caught:
        sources.catalogue(path)
    message = str(caught.value)
    assert path.name in message, message
    assert f"source.{USER_ID}" in message, message
    for word in says:
        assert word in message, (word, message)
    return message


class TestT2FileRefusals:
    @pytest.fixture
    def user(self, tmp_path: Path) -> Path:
        return tmp_path / "my-sources.toml"

    @pytest.mark.parametrize("kind", BASES)
    def test_control_each_base_entry_loads(
        self, sources: ModuleType, user: Path, kind: str
    ) -> None:
        """The refusals below each change one thing in one of these."""
        got = sources.catalogue(write(user, {f"source.{USER_ID}": BASES[kind]}))
        assert got.sources[USER_ID].id == USER_ID
        assert got.sources[USER_ID].kind == kind

    @pytest.mark.parametrize(
        ("base", "drop", "update", "says"),
        [
            pytest.param(ONE_COG, (), {"colour": "red"}, ("colour",), id="an unknown key"),
            pytest.param(ONE_COG, (), {"id": USER_ID}, ("id",), id="an id key in the table"),
            pytest.param(ONE_COG, ("crs",), {}, ("crs",), id="crs missing"),
            pytest.param(ONE_COG, (), {"kind": "two-cogs"}, ("kind",), id="an unknown kind"),
            pytest.param(PER_YEAR, ("years",), {}, ("years",), id="per-year without years"),
            pytest.param(
                PER_YEAR, (), {"url_template": "https://example.invalid/c.tif"}, ("year",),
                id="per-year without {year}",
            ),
            pytest.param(TILES, ("tile_list_url",), {}, ("tile_list_url",), id="tiles, no list"),
            pytest.param(TILES, ("tile_name",), {}, ("tile_name",), id="tiles, no tile_name"),
            pytest.param(
                TILES, ("tile_spacing",), {}, ("tile_spacing",), id="tiles, no tile_spacing"
            ),
            pytest.param(
                TILES, (), {"url_template": "https://example.invalid/t.tif"}, ("url_template",),
                id="tiles without {name}",
            ),
            pytest.param(
                ONE_COG, (), {"tile_name": GLO30_TILE_NAME}, ("tile_name",),
                id="tile_name on one-cog",
            ),
            pytest.param(
                ONE_COG, (), {"tile_spacing": list(GLO30_TILE_SPACING)}, ("tile_spacing",),
                id="tile_spacing on one-cog",
            ),
            pytest.param(
                TILES, (), {"tile_spacing": [0.0, 0.00027778]}, ("tile_spacing",),
                id="tile_spacing 0",
            ),
            pytest.param(
                TILES, (), {"tile_spacing": [0.0027778, -0.00027778]}, ("tile_spacing",),
                id="tile_spacing negative",
            ),
            pytest.param(
                ONE_COG, (), {"credit": "fetched {when}"}, ("credit", "when"), id="credit {when}"
            ),
            pytest.param(
                ONE_COG, (), {"credit": "in {date.year}"}, ("credit",), id="credit {date.year}"
            ),
            pytest.param(
                ONE_COG, (), {"credit": "from {url[0]}"}, ("credit",), id="credit {url[0]}"
            ),
            pytest.param(
                TILES, (), {"tile_name": "tile_{x}"}, ("tile_name",), id="tile_name {x}"
            ),
            pytest.param(
                PER_YEAR, (), {"z_unit": "metre"}, ("z_unit",), id="z_unit on a classes entry"
            ),
        ],
    )  # fmt: skip
    def test_the_entry_is_refused_naming_file_key_path_and_field(
        self,
        sources: ModuleType,
        user: Path,
        base: dict[str, Any],
        drop: tuple[str, ...],
        update: dict[str, Any],
        says: tuple[str, ...],
    ) -> None:
        refusal(sources, write(user, {f"source.{USER_ID}": changed(base, *drop, **update)}), *says)

    def test_a_toml_syntax_error_carries_the_line(self, sources: ModuleType, user: Path) -> None:
        user.write_text('[source.my-x]\nkind = "one-cog"\ncrs = EPSG:25833\n', encoding="utf-8")
        with pytest.raises(ValueError) as caught:
            sources.catalogue(user)
        message = str(caught.value)
        assert user.name in message and "line 3" in message, message


# ------------------------------------------------------------ T3


class TestT3UserEntries:
    def test_no_user_file_is_the_shipped_catalogue(self, sources: ModuleType) -> None:
        got = sources.catalogue(None)
        assert dict(got.sources) == dict(sources.SOURCES)
        assert dict(got.stations) == dict(sources.STATION_SOURCES)

    def test_an_absent_user_file_is_the_shipped_catalogue(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        got = sources.catalogue(tmp_path / "sources.toml")
        assert dict(got.sources) == dict(sources.SOURCES)

    def test_a_user_entry_joins_the_shipped_ones(self, sources: ModuleType, tmp_path: Path) -> None:
        user = write(tmp_path / "sources.toml", {f"source.{USER_ID}": ONE_COG})
        got = sources.catalogue(user)
        assert sorted(got.sources) == sorted([*SHIPPED, *FEATURE_LAYERS, USER_ID])
        assert got.sources[USER_ID].crs == "EPSG:25833"
        assert got.sources[USER_ID].nodata is None
        assert sorted(got.stations) == sorted(STATIONS)

    def test_the_catalogue_reads_sources_at_call_time(
        self, sources: ModuleType, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """4.3: the suites that put an entry into `SOURCES` keep working."""
        extra = sources.RemoteSource(id="planted", **ONE_COG)
        monkeypatch.setattr(sources, "SOURCES", {**sources.SOURCES, "planted": extra})
        assert sources.catalogue(None).sources["planted"] == extra

    def test_a_shipped_id_is_refused_naming_both_files(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        user = write(tmp_path / "sources.toml", {"source.glo30": TILES})
        with pytest.raises(ValueError) as caught:
            sources.catalogue(user)
        message = str(caught.value)
        assert "glo30" in message and message.count("sources.toml") >= 2, message

    def test_a_station_table_is_refused(self, sources: ModuleType, tmp_path: Path) -> None:
        station = {
            "service_url": "https://example.invalid/services", "list_file": "x.csv",
            "crs": "EPSG:25833", "credit": "c", "licence_note": "l",
        }  # fmt: skip
        user = write(tmp_path / "sources.toml", {"station.x": station})
        with pytest.raises(ValueError) as caught:
            sources.catalogue(user)
        message = str(caught.value)
        assert user.name in message and "station" in message, message

    def test_a_relative_classes_is_beside_the_user_file(
        self, sources: ModuleType, tmp_path: Path
    ) -> None:
        user = write(
            tmp_path / "root" / "sources.toml",
            {f"source.{USER_ID}": changed(PER_YEAR, classes="maps/x.toml")},
        )
        got = sources.catalogue(user).sources[USER_ID]
        assert Path(got.classes).resolve() == (tmp_path / "root" / "maps" / "x.toml").resolve()


class TestT3DataRoot:
    def test_data_root_reads_rasputin_data(self) -> None:
        assert cli.data_root({"RASPUTIN_DATA": "/some/root"}) == Path("/some/root")
        assert cli.data_root({}) is None
        assert cli.data_root({"RASPUTIN_DATA": ""}) is None

    def test_the_cache_root_derives_from_it(self) -> None:
        assert cli.cache_root(None, {"RASPUTIN_DATA": "/r"}) == Path("/r") / "cache"
        assert cli.cache_root(Path("/c"), {"RASPUTIN_DATA": "/r"}) == Path("/c")

    @pytest.mark.parametrize("command", ["fetch", "mesh"])
    def test_the_help_names_the_user_file(self, command: str) -> None:
        """The user's file is found in more than one place, so the help points
        at the listing that names the one read, not at one variable."""
        text = plain(runner.invoke(app, [command, "--help"]).output)
        assert squashed("your sources file (see rasputin sources)") in squashed(text), text
        assert squashed("$RASPUTIN_DATA/sources.toml") not in squashed(text), text


class TestT3ThroughTheCli:
    @pytest.fixture
    def server(self) -> Any:
        with RangeServer() as running:
            yield running

    def test_fetch_and_mesh_a_user_entry(
        self,
        server: RangeServer,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """`fetch my-dem` copies the served file into `$RASPUTIN_DATA/cache`;
        `mesh --dem my-dem` reads it offline and equals the mesh of the file
        by path."""
        data = projected()
        root = tmp_path / "data"
        entry = {
            "kind": "one-cog", "url": server.put("p.tif", data), "crs": PROJECTED_CRS,
            "nodata": -9999.0, "credit": "My credit", "licence_note": "My licence",
        }  # fmt: skip
        write(root / "sources.toml", {"source.my-dem": entry})
        monkeypatch.setenv("RASPUTIN_DATA", str(root))
        domain = geojson(tmp_path / "d.geojson", TRIANGLE)
        code, output = invoke("fetch", "my-dem", "--domain", str(domain))
        assert code == 0, output
        assert (root / "cache" / "my-dem" / "manifest.json").is_file()
        assert "My credit" in (root / "cache" / "my-dem" / "NOTICE.txt").read_text()
        tif = tmp_path / "dem.tif"
        tif.write_bytes(data)
        by_path, by_key = tmp_path / "path.vtk", tmp_path / "key.vtk"
        mesh = ("--domain", str(domain), "--tolerance", "1")
        code, output = invoke("mesh", "--dem", str(tif), *mesh, "--out", str(by_path))
        assert code == 0, output
        monkeypatch.setattr(socket, "socket", no_network)
        monkeypatch.setattr(socket, "create_connection", no_network)
        code, output = invoke("mesh", "--dem", "my-dem", *mesh, "--out", str(by_key))
        assert code == 0, output
        same_mesh(read_vtk(by_key.read_bytes()), read_vtk(by_path.read_bytes()))

    def test_a_bad_user_file_is_a_usage_error_naming_it(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        root = tmp_path / "data"
        write(root / "sources.toml", {"source.my-dem": changed(ONE_COG, "crs")})
        monkeypatch.setenv("RASPUTIN_DATA", str(root))
        domain = geojson(tmp_path / "d.geojson", TRIANGLE)
        target = tmp_path / "x.vtk"
        code, output = invoke(
            "mesh", "--dem", "my-dem", "--domain", str(domain), "--tolerance", "1",
            "--out", str(target),
        )  # fmt: skip
        assert code == USAGE, output
        # Squashed: Rich breaks the long temporary path inside its panel at a
        # column set by the terminal width ("sources. toml" on CI's runner).
        assert "sources.toml" in squashed(output), output
        assert "source.my-dem" in squashed(output), output
        assert "Traceback" not in squashed(output), output
        assert not target.exists()

    def test_cache_alone_loads_no_user_entries(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """4.3: `--cache` says where tiles are, not where the user file is."""
        monkeypatch.delenv("RASPUTIN_DATA", raising=False)
        root = tmp_path / "data"
        write(root / "sources.toml", {"source.my-dem": ONE_COG})
        domain = geojson(tmp_path / "d.geojson", TRIANGLE)
        code, output = invoke(
            "mesh", "--dem", "my-dem", "--cache", str(root / "cache"), "--domain", str(domain),
            "--tolerance", "1", "--out", str(tmp_path / "x.vtk"),
        )  # fmt: skip
        assert code == USAGE, output
        assert "my-dem" in squashed(output), output


# ------------------------------------------------------------ T6


@pytest.mark.parametrize(
    ("lat", "lon", "name"),
    [
        (60, 10, "Copernicus_DSM_COG_10_N60_00_E010_00_DEM"),
        (-15, -45, "Copernicus_DSM_COG_10_S15_00_W045_00_DEM"),
        (0, 0, "Copernicus_DSM_COG_10_N00_00_E000_00_DEM"),
        (-1, -1, "Copernicus_DSM_COG_10_S01_00_W001_00_DEM"),
        (89, 179, "Copernicus_DSM_COG_10_N89_00_E179_00_DEM"),
        (-90, -180, "Copernicus_DSM_COG_10_S90_00_W180_00_DEM"),
    ],
)
def test_t6_glo30_tile_names(sources: ModuleType, lat: int, lon: int, name: str) -> None:
    assert sources.SOURCES["glo30"].tile_for(lat, lon) == name
