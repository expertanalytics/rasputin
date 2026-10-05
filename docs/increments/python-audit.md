# Audit: the Python layer — architecture, duplication, dead code

Status: **audit, not a design.** Written by `@architect` against master at
`12dace7` (its `src_python/` is byte-identical to `bc01cd8`). Read-only on
code. Every `file:line` citation is pinned to `12dace7`, so it stays true
while the code moves. It proposes a sequence of refactor PRs; each still runs the normal loop
(`docs/increments/README.md`): design note, red tests where behaviour changes,
code, review, and `@perf` where the diff touches what drives refine or mesh.

Ola's ask: look at the architecture, modularization and class hierarchies,
and above all at code doing almost the same thing in different places without
sharing logic.

**Prior art.** None applies: this audits shipped code and claims nothing new.
A PR below that designs a new module carries its own prior-art section.

## Verdict in five lines

1. The duplication is real but smaller in lines than it looks: about **450
   production lines** (6 % of 7,732 code lines) and **460 test lines** would
   go. The larger cost is **drift**: copies have already diverged, once into a
   bug (F3) and once into an inconsistency with the C++ core (F9).
2. `cli.py` (1,684 code lines, 22 % of the package) is a pipeline host, not a
   composition root: the whole `mesh` run, the catchment placement and five
   file encoders live in it. That is the main structural fault.
3. There is no shared "lattice" vocabulary: node-to-coordinate arithmetic on
   `RasterMeta` is spelt out by hand at about 30 sites in 12 modules.
4. Three GeoJSON readers apply three different rules to the `crs` member, and
   two functions named `read_lakes` read lakes two ways.
5. The tests repeat a CLI harness in about 20 files, and pin six private `cli`
   names, which blocks the cli split until they are re-pointed.

## 1. Size, re-measured

Non-blank lines (comments and docstrings included), on `12dace7`:
`src_python` 9,924; `tests/python` 44,988; `tools/` 3,180 and
`.claude/hooks/` 653 (3,833 together). The main session's figures reproduce.

Code lines only (no blank, comment or docstring line; the count of
`CLAUDE.md` §2), by a tokenizer script: `src_python` **7,732**, of which
`cli.py` 1,684, `mosaic.py` 404, `catchment.py` 378, `feature_input.py` 370,
`io/geotiff.py` 265, `viz/fixtures.py` 244, `io/repository.py` 233. Tests
36,523. Tools and hooks 3,231.

Where the tests are (non-blank lines): `test_cli_*` 9,533; `test_io_*` 4,306;
`test_core_*` 4,135; `test_viz_*` 2,392; `test_mosaic.py` 1,546; the tests of
the harness tools (`bench`, `brief`, `away`, `session_state`, the guards,
`ci_changes`, `check_citations`, settings, rule sizes) **6,550**, more than the
3,833 lines they test.

## 2. Findings, ranked by lines saved times risk

Each finding: evidence (file:lines on `12dace7`), proposed shape, lines that
would go (net production lines, the `CLAUDE.md` §2 count), risk, and the PR
in section 6 that takes it.

### F1. CLI options and refusal plumbing repeated — about 100 lines, low risk (PR G)

- 17 options are declared in two to four commands: `--out-parent` four times
  (`src_python/tin_engine/cli.py@12dace7:421-424, 620-623, 1829-1832, 2105-2108`), `--dem` three times,
  `--snap-spacing`, `--outline-tolerance`, `--map-radius`, `--reach-up`,
  `--rivers`, `--lakes`, `--domain`, `--domain-crs`, `--bbox`, `--out-crs`,
  `--cache`, `--refresh`, `--out-dir`, `--delaunay` twice each. 92 declaration
  lines beyond each first (AST count).
- "Finite and positive / at least 0" is checked by hand 7 times
  (`src_python/tin_engine/cli.py@12dace7:381, 835, 846-852, 1866-1870, 1875-1878, 2114-2116, 2117-2120`),
  and the outline tolerance a third time in `CatchmentRequest`
  (`src_python/tin_engine/catchment.py@12dace7:106-108`).
- 16 `try / except X / raise typer.BadParameter(str(exc), param_hint=...)`
  blocks (`src_python/tin_engine/cli.py@12dace7:865-868, 1278-1281, 1389-1412, 1459-1462, 1735-1736,
  1753-1756, 1882-1885, 1903-1907, 2122-2125, 2140-2145, ...`).
- The "no cache: set RASPUTIN_DATA ..." refusal twice (`src_python/tin_engine/cli.py@12dace7:1264-1269,
  1363-1368`); `--bbox` tuple to `Bounds` twice (`src_python/tin_engine/cli.py@12dace7:1272-1273,
  1383-1388`); `--map-radius`/`--reach-up` defaults given twice, as `or 500.0`
  in `src_python/tin_engine/cli.py@12dace7:1892` and as option defaults in `src_python/tin_engine/cli.py@12dace7:2083-2092`.

Shape: one `cli_options.py` (or a block at the top of `cli.py`) of
`Annotated` aliases (`OutParent`, `DemPaths`, `OutlineTolerance`,
`MapRadius`, ...) with Typer callbacks for the numeric rules; a
`refusing(hint, *errors)` context manager that turns a library refusal into
`BadParameter` naming the flag; `Bounds.of(tuple)`. Risk: low; the CLI suites
pin every wording and `param_hint`, so a slip shows.

### F2. Lattice arithmetic spelt out at ~30 sites — about 100 lines, medium risk (PR A)

`RasterMeta` (`src_python/tin_engine/io/models.py@12dace7:39-82`) holds the grid but no method, so every
module re-derives it:

- node rectangle `x_min + (cols-1)·dx`, `y_max - (rows-1)·dy`:
  `src_python/tin_engine/domain.py@12dace7:146-147`, `src_python/tin_engine/dem_input.py@12dace7:170-171, 257`, `src_python/tin_engine/catchment.py@12dace7:399-409`,
  `src_python/tin_engine/catchment_batch.py@12dace7:159-164`, `src_python/tin_engine/mosaic.py@12dace7:375-378`; and `TargetGrid`'s own
  twice (`src_python/tin_engine/target_grid.py@12dace7:117-118, 232-233`);
- node to coordinate: `src_python/tin_engine/grid_domain.py@12dace7:68-69`, `src_python/tin_engine/catchment.py@12dace7:323, 394, 461,
  467`, `src_python/tin_engine/burn.py@12dace7:169, 223`, `src_python/tin_engine/mosaic.py@12dace7:458, 466-467, 563`,
  `src_python/tin_engine/target_grid.py@12dace7:244, 251`, `src_python/tin_engine/cli.py@12dace7:1703-1705`;
- coordinate to index, with three rounding rules: `src_python/tin_engine/cli.py@12dace7:1701-1702`,
  `src_python/tin_engine/catchment.py@12dace7:376-377, 383-390`, `src_python/tin_engine/burn.py@12dace7:71-72, 151`,
  `src_python/tin_engine/target_grid.py@12dace7:179-181`, `src_python/tin_engine/reference.py@12dace7:76-77`, and `src_python/tin_engine/fetch/plan.py@12dace7:149-152`
  (the box snapped to the lattice a fourth way, though the fetch must fetch a
  superset of what `mosaic._snapped` reads, `src_python/tin_engine/mosaic.py@12dace7:340-352`);
- `Bounds` built from four keys by `dict(zip(...))` four times
  (`src_python/tin_engine/dem_input.py@12dace7:232-234`, `src_python/tin_engine/catchment.py@12dace7:408-409`, `src_python/tin_engine/cli.py@12dace7:1272-1273,
  1387`).

`src_python/tin_engine/grid_domain.py@12dace7:53-57` warns that another spelling of the node expression
can land an ulp off the node. That is the argument for one spelling.

Shape: on `RasterMeta`: `node_xy(rows, cols)` (vectorised, the exact
`x_min + c·dx`, `y_max - r·dy` spelling), `index_of(x, y)` (fractional),
`node_bounds() -> Bounds`, `cell_area`; on `IndexWindow`: `meet`, `within`,
`clamp`, `from_inclusive` (now `src_python/tin_engine/mosaic.py@12dace7:355-362, 491-503`, which also mixes
inclusive `(r0, c0, r1, c1)` tuples with `IndexWindow`). Move `Bounds`
(`src_python/tin_engine/mosaic.py@12dace7:64-81`) and `TileFootprint` (`src_python/tin_engine/io/repository.py@12dace7:53-61`) into
`io/models.py`: that removes `fetch/plan.py`'s import of `mosaic` for a value
type and the `io.repository` <-> `mosaic` import cycle (type-only today,
`src_python/tin_engine/io/repository.py@12dace7:44-47`, `src_python/tin_engine/mosaic.py@12dace7:33-34`). Risk: medium. Bit-exact node
coordinates feed `refine`, so the PR must leave every mesh byte-identical
(refine golden suite, `@perf`'s run on the 1 m set).

### F3. CRS comparisons repeated, and one has become a bug — about 25 lines, low risk (PR B)

- `parse_crs(a) != parse_crs(b)` at 8 sites (`src_python/tin_engine/catchment.py@12dace7:185, 202`,
  `src_python/tin_engine/feature_input.py@12dace7:179, 423`, `src_python/tin_engine/cli.py@12dace7:932, 954, 2045`, `src_python/tin_engine/fetch/run.py@12dace7:222`).
- `"none" if same else transform_description(...)` twice in one function
  (`src_python/tin_engine/cli.py@12dace7:931-934, 952-956`).
- "The tiles are in one CRS" three times, keyed two ways: by `crs` text in
  `src_python/tin_engine/dem_input.py@12dace7:188-190` and `src_python/tin_engine/catchment.py@12dace7:194-196`, by `epsg` in
  `src_python/tin_engine/dem_input.py@12dace7:227-231`. **The `epsg` copy is a bug:** a DEM whose CRS has
  no EPSG code, meshed with `--domain`, is refused with `cannot read the CRS
  'EPSG:None'`. Probe (run with the project venv, which loads this code):
  `_domain_plan` on a one-tile footprint with `epsg=None` and the WKT of UTM
  33N without its ID raises exactly that `ValueError`; the control with
  `epsg=25833` plans a 5 x 5 grid. The box path (`plan_mosaic` without a
  domain) does not have the bug.

Shape: `crs.same_crs(a, b)`, `crs.transform_label(src, dst)` ("none" or
PROJ's description), and one `single_crs(footprints) -> str` used by all
three. The bug fix needs a red test first (`@tester`).

### F4. Catchment placement written twice — about 50 lines, low risk (PR F)

`rasputin catchment --rivers` and `run_batch` do the same per-station work
in two places:

- place the gauge and build the request: `src_python/tin_engine/cli.py@12dace7:1709-1728, 1886-1902`
  against `src_python/tin_engine/catchment_batch.py@12dace7:177-204`;
- the gauge's report fields, built from `GaugeResult.__slots__` with the same
  `skip` tuple: `src_python/tin_engine/cli.py@12dace7:1782-1785` against `src_python/tin_engine/catchment_batch.py@12dace7:132-137`;
- `check_reach_crs` runs twice on `station-catchments` (`src_python/tin_engine/cli.py@12dace7:2146`, then
  `src_python/tin_engine/catchment_batch.py@12dace7:157`), and the unknown `--only` check twice
  (`src_python/tin_engine/cli.py@12dace7:2126-2130`, `src_python/tin_engine/catchment_batch.py@12dace7:154-156`);
- `Gauge` (`src_python/tin_engine/gauge.py@12dace7:45-53`) is a projection of `Station`
  (`src_python/tin_engine/io/station_set.py@12dace7:26-39`), copied field by field (`src_python/tin_engine/catchment_batch.py@12dace7:178`).

Shape: `GaugeResult.fields()` (one method), and one
`catchment_batch.request_for(x, y, placement, seed, ...)` that the single
command also calls. The CLI keeps the usage-error versions of the two checks;
`run_batch` keeps its own as a library precondition, or takes them as
already checked — one or the other, stated.

### F5. Three GeoJSON readers, three `crs` rules; two `read_lakes` — about 40 lines, medium risk (PR C)

- `domain._geojson` (`src_python/tin_engine/domain.py@12dace7:126-139`): no `crs` member means EPSG:4326;
  reads the file itself (`src_python/tin_engine/domain.py@12dace7:92`).
- `feature_input.read_source` (`src_python/tin_engine/feature_input.py@12dace7:403-414`): same default,
  its own parse of the member and of the features, reads the file itself.
- `io/station_set.features_of` (`src_python/tin_engine/io/station_set.py@12dace7:52-65`): the member is
  required; reads through `repository.read_json`.
- Two functions named `read_lakes`: `src_python/tin_engine/feature_input.py@12dace7:429-446` (any polygon
  source near a point, returns bare geometries) and `src_python/tin_engine/io/station_set.py@12dace7:126-142`
  (NVE's lakes file, returns `Lake` records, a MultiPolygon split). The CLI
  imports the second under an alias to dodge the clash (`src_python/tin_engine/cli.py@12dace7:2137`), and
  `--lakes` means a different file format in `catchment` and in
  `station-catchments` (`src_python/tin_engine/cli.py@12dace7:1808-1815` against `2075-2082`).
- Writing: `io/geojson.catchment_geojson` (`src_python/tin_engine/io/geojson.py@12dace7:18-35`) and
  `fetch/nve._collection` (`src_python/tin_engine/fetch/nve.py@12dace7:115-121`) build the same
  FeatureCollection-with-`crs`-member document.

Shape: `io/geojson.py` owns both directions: `read_collection(data,
*, crs_required) -> (features, crs_text)` and `feature_collection(crs,
features) -> dict`. The two lake readers keep their behaviours, renamed by
what they read (`read_lake_polygons` near a point; `read_nve_lakes`), and
share the polygon-part filter. Risk: medium; the refusal wordings differ today
and the suites pin them, so the PR either keeps each wording or lands a
`@tester` amendment.

### F6. Mesh topology joins duplicated between `cli.py` and `viz/scene.py` — about 35 lines, low risk (PR E)

`cli._undirected`, `_masked_pairs`, `_chain_masks` (`src_python/tin_engine/cli.py@12dace7:504-552`) do what
`viz/scene._key`, `_mesh_edges`, `_chain_edges` (`src_python/tin_engine/viz/scene.py@12dace7:111-176`) do:
the same mask-bit convention (bit `e` is edge `(v[e], v[(e+1)%3])`), the same
undirected key, the same union of properties. `src_python/tin_engine/cli.py@12dace7:560-563` says it may not
call `viz/` (ruling 6), so it copied instead.

Shape: a neutral pure module, `topology.py` (numpy, typed against
`viz/protocols.py`'s `MeshLike`/`PslgLike`, no `_core`, no `viz`), used by
both. The convention is then written once. This changes `viz/scene.py`'s
import rule (today it imports only `protocols`), so the isolation test in
`tests/python/test_viz_scene.py@12dace7:366-405` changes in the same PR; the rule's purpose (no
`_core` in `viz/`) still holds.

### F7. File encoding still in `cli.py`, and the two mesh writers share no gate — about 30 lines, low risk (PR D)

- `cli.py` encodes five outputs itself: the palette JSON (`src_python/tin_engine/cli.py@12dace7:1105`),
  `results.csv` rows and header (`src_python/tin_engine/cli.py@12dace7:1985-1995, 2029-2033`),
  `summary.json` (`src_python/tin_engine/cli.py@12dace7:2182-2184`), the fetched station files
  (`src_python/tin_engine/cli.py@12dace7:1337-1339`), and the `.ply` comment lines (`src_python/tin_engine/cli.py@12dace7:1065-1067`).
  `project_structure.md` lists the first three as known exceptions to "all
  encoding lives in `io/`".
- `src_python/tin_engine/io/ply.py@12dace7:114-131` and `src_python/tin_engine/io/vtk_legacy.py@12dace7:192-200` each have the
  control-character-and-ASCII gate; both check land-cover codes fit int32
  (`src_python/tin_engine/io/ply.py@12dace7:174-177`, `src_python/tin_engine/io/vtk_legacy.py@12dace7:111-115`); both build the sorted
  vocabulary table. `cli._ascii` (`src_python/tin_engine/cli.py@12dace7:1475-1477`) and `crs.crs_label`
  (`src_python/tin_engine/crs.py@12dace7:80`) both escape to ASCII by `backslashreplace`.

Shape: `io/_text.py` (`checked_ascii`, `escaped_ascii`, `int32_codes`), and
the five encoders moved into `io/` (`io/tables.py` or beside their readers).

### F8. `cli.py` hosts the mesh pipeline — about 60 lines saved, more moved; high risk (PR H)

`src_python/tin_engine/cli.py@12dace7:776-1062` (the `mesh` body), `1150-1217` (report and fixture mesh),
`1487-1706` (`_DemMesh`, `_Once`, `_dem_mesh`, `_refine_phases`,
`_off_node`) and `277-317` (`_engine`) are the whole DEM-to-mesh run: start
chains, engine, sample or refine, edge strip, final check, trim, and the
assembly of ~45 record values (`src_python/tin_engine/cli.py@12dace7:906-970, 1551-1672`). None of it is
about flags. It cannot run from a GUI or API worker without Typer, which the
architect's pillar 4 and the python skill both require, and its refusals are
`typer.BadParameter` raised deep inside (`src_python/tin_engine/cli.py@12dace7:1581, 1612, 1631, 1669`).

Also: `_Run` and `Attempt` (`src_python/tin_engine/cli.py@12dace7:235-275`) are two records of one engine
run; `_DemMesh`, `Sizes` (`src_python/tin_engine/stats.py@12dace7:108-126`) and the `values` dict carry the
same run facts three ways; `_write_report` rebuilds `run_record.stats_rows`
by hand (`src_python/tin_engine/cli.py@12dace7:1177`, against `src_python/tin_engine/run_record.py@12dace7:168-170`, which only tests use).

Shape: `mesh_run.py`: `MeshRequest` (frozen Pydantic: tile handle, stride,
tolerance, domain, features, flags as values) and `run(request, clock) ->
MeshResult` (the trimmed mesh, the record, the report sizes), raising
`MeshRefusal(flag, message)`; `start_mesh.py` for `_engine` and the
fixture path. `cli.mesh` becomes flags to request, run, write, print.
Moving ~550 lines is near net zero under the net rule. Risks: `tools/bench.py`
patches `refine` by name inside `cli`'s namespace (`tools/bench.py@12dace7:750-761`),
so the seam moves with the code and `@perf` re-checks the bench tool;
tests reach `cli._dem_mesh`, `_engine`, `_triangulated`, `_chain_masks`
(section 3, X2); and the in-flight 23c-2 branch (`worktree-23c`, not merged)
adds about 250 lines to `cli.py` and a `pieces.py`. **This PR waits for
23c-2 to merge.**

### F9. Three NoData rules — about 8 lines, low risk, one behaviour question (PR A)

`mosaic._valid` (`src_python/tin_engine/mosaic.py@12dace7:506-510`) and `src_python/tin_engine/burn.py@12dace7:140` treat NaN and the
sentinel as NoData, as the core does (`include/terrain/raster/raster.hpp@12dace7:63-66`).
`target_grid._valid` (`src_python/tin_engine/target_grid.py@12dace7:153-155`) uses `isfinite`, so it also
treats +-inf as NoData. Probe: on `[1.0, inf, nan]` with no sentinel,
`target_grid._valid` gives `[True, False, False]`, `mosaic._valid`
`[True, True, False]`. Shape: one `valid_mask(values, nodata)` beside
`RasterMeta`. Which rule wins is a question for Ola (default below).

### F10. Types that say `Any` where a type exists — 0 lines, low risk (in PRs A and F)

`DemRepository` (`src_python/tin_engine/io/repository.py@12dace7:64-74`) lacks `load_window` and `check`,
which every caller uses, so `dem_input._open_reprojected` takes
`repository: Any, footprints: Any` (`src_python/tin_engine/dem_input.py@12dace7:183`), and `catchment`
passes `footprints: Any` four times (`src_python/tin_engine/catchment.py@12dace7:226, 285, 362`);
`cli._placed(repository: Any)` (`src_python/tin_engine/cli.py@12dace7:1711, 1731`);
`CatchmentRequest.lakes: tuple[Any, ...]` (`src_python/tin_engine/catchment.py@12dace7:89`);
`reference.summarise(rows: Sequence[Any])` reads `StationResult` by attribute
(`src_python/tin_engine/reference.py@12dace7:198-265`). Widen the Protocol and name the types.

### F11. Dead code and checks that guard nothing reachable — about 15 lines, low risk (in PRs G, H)

- `import asyncio` inside `fetch` (`src_python/tin_engine/cli.py@12dace7:1252`) repeats the module import
  (`src_python/tin_engine/cli.py@12dace7:42`).
- The two-code-systems refusal (`src_python/tin_engine/cli.py@12dace7:1439-1448`) cannot fire: its own
  comment says every coded map is CORINE, and `CLASS_MAPS` is closed in code
  (`src_python/tin_engine/feature_input.py@12dace7:102-118`). A one-line test over `CLASS_MAPS` guards the
  same future change at no runtime cost.
- `PIECE_VOCABULARY` (`src_python/tin_engine/features.py@12dace7:189-191`) has no user on master; 23c-2
  uses it. Not dead, in flight.
- `RemoteSource` and `StationSource` (`src_python/tin_engine/sources.py@12dace7:22-40, 83-96`) repeat five
  fields; a base model `CreditedSource` holds them (about 6 lines).
- Results are frozen dataclasses in some modules (`src_python/tin_engine/stats.py@12dace7:7-8` states the
  rule) and frozen Pydantic models in others (`gauge.py`, `reference.py`,
  `src_python/tin_engine/fetch/run.py@12dace7:47`). Not wrong; the rule should be written once in the
  python skill rather than in a module docstring.

## 3. Tests

- **X1. A CLI harness in about 20 files — about 350 test lines (PR T1).**
  Bit-identical helper copies alone are 156 lines (AST comparison over
  `tests/python/`, docstrings ignored): `bumpy` 8 copies differing only in the
  seed (`tests/python/test_cli_mesh_domain.py@12dace7:134-136`, `tests/python/test_cli_mesh_features.py@12dace7:126-128`,
  ...), `invoke` 8 (`tests/python/test_cli_catchment.py@12dace7:93-95`, `tests/python/test_cli_mesh_dem.py@12dace7:53-55`,
  ...), `run` 7 variants, `refused` 5 variants
  (`tests/python/test_cli_catchment.py@12dace7:137-144`, `tests/python/test_cli_mesh_domain_crs.py@12dace7:135-143`,
  `tests/python/test_cli_mesh_geographic.py@12dace7:149-157`, `tests/python/test_cli_mesh_plain_output.py@12dace7:108-115`),
  `squashed` 3, `plain_square`/`square` 6, `write_geojson` 3, a module-level
  `CliRunner()` in 17 files. Shape: `tests/python/cli_harness.py` with
  `invoke(command, *args)`, `ran(...)`, `refused(..., says=...)`, and a
  `rough_dem(seed)` fixture factory.
- **X2. Tests that pin private `cli` names — 0 lines, but they block F6 and
  F8.** `cli._triangulated` (4 uses), `cli._chain_masks` (4), `cli._engine`
  (2), `cli._destination` (2), `cli._dem_mesh` (monkeypatched,
  `tests/python/test_cli_mesh_domain_crs.py@12dace7:299`), `cli._undirected` (1), in
  `feature_fixtures.py`, `test_cli_mesh.py`, `test_cli_mesh_vtk.py`,
  `test_cli_draw.py`, `test_cli_mesh_domain_crs.py`. Each PR that moves the
  function re-points the test to the public name in the same PR, as a
  `@tester` commit.
- **X3. Import-firewall tests scattered over 8 files — about 110 test lines
  (PR T2).** `tests/python/test_viz_svg.py@12dace7:536-594`, `tests/python/test_viz_scene.py@12dace7:366-405`,
  `tests/python/test_features.py@12dace7:580-608`, `tests/python/test_io_geotiff.py@12dace7:1303`,
  `tests/python/test_landcover.py@12dace7:450-453`, `tests/python/test_outline.py@12dace7:139`, `tests/python/test_io_ply.py@12dace7:439`,
  `tests/python/test_io_vtk_legacy.py@12dace7:492`; most of the 52 modules have none.
  `tests/python/test_features.py@12dace7:608` greps the source text for `_core`, which pins prose,
  not imports. Shape: one table-driven `test_layering.py` that declares, per
  module, the first-party modules it may import (section 5's layers) and
  checks it with `importscan.first_party_imports`. It becomes the executable
  dependency map, and every later PR in section 6 is checked by it.
- Not proposed: rewriting the 9,533 lines of `test_cli_*` against
  `mesh_run.py` after F8. They test behaviour through the command, which is
  what Ola runs; moving them is churn.

## 4. tools/ (light)

- Five private git runners: `_git` in `tools/brief.py@12dace7:115-118`,
  `tools/rule_sizes.py@12dace7:63-67`, `tools/session_state.py@12dace7:119-124`; `git` in
  `tools/check_citations.py@12dace7:111-112` and `.claude/hooks/gates_after_commit.py@12dace7:35-41`
  (about 20 lines). A shared `tools/_git.py`.
- Two lists of governed files that already differ:
  `tools/check_citations.py@12dace7:92-99` (has `testing.md`) and
  `.claude/hooks/guard_governance.py@12dace7:51-76` (has `tools/away.py` and more,
  not `testing.md`). They may answer different questions; if so each list
  says so, else one imports the other. Governed files: Ola approves each edit.
- Two VTK readers: `tools/bench.py@12dace7:336` (`read_vtk_ascii`) and
  `tests/python/vtkread.py`.
- No tool duplicates `src_python` logic except `bench.py`, which reuses
  `tin_engine.stats.quality` (good) and patches `cli.refine` (finding F8).

## 5. Target architecture (one page)

Six layers. A module imports only from its own layer or below;
`test_layering.py` (X3) enforces it.

```
L5  cli.py            flags -> request; refusal -> BadParameter(flag);
                      result -> files and stderr lines. No algorithm.
L4  pipelines         dem_input, feature_input, mesh_run (new), catchment,
                      catchment_batch, fetch/run, fetch/nve
                      request (frozen Pydantic) in, result out; no typer,
                      no print, no path but what the request names
L3  core adapters     raster.to_core (the one raster adapter), start_mesh
                      (build_pslg -> node -> triangulate), edge_strip,
                      final_check: the only importers of _core
L2  io/ codecs        bytes <-> values: geotiff, cog, geopackage, gml,
                      geojson (read AND write), ply, vtk_legacy, tables
                      (csv/json/palette), mesh_index, station_set, rivers;
                      repository.py the one module that opens files
L1  pure algorithms   crs (same_crs, transform_label, single_crs), mosaic,
                      target_grid, grid_domain, domain, chains, elevation,
                      outline, burn, gauge, sensitivity, reference,
                      landcover, decompose, topology (new), viz/*
L0  values            io/models (RasterMeta with node methods, IndexWindow
                      with window methods, Bounds, TileFootprint, DemTile,
                      valid_mask), features, sources, run_record, stats,
                      palettes
```

What changes against today: `Bounds`, `TileFootprint` and the lattice
methods move down to L0, so `fetch/`, `target_grid`, `catchment` and
`io/repository` stop importing `mosaic` for value types; the GeoJSON reading
in `domain.py` and `feature_input.py` moves into `io/geojson.py`; `cli.py`
loses the mesh run (to `mesh_run.py`), the engine wrapper (to
`start_mesh.py`), the topology joins (to `topology.py`), the catchment
placement (to `catchment_batch`) and five encoders (to `io/`). `cli.py`
should end near 900 code lines, all of it flags, refusal mapping and
printing. Stderr wording that is domain text (the catchment's window lines,
the land-cover line) moves to pure wording functions beside the result type,
as `run_record.summary` already does for the mesh.

## 6. PR order

Each PR is a refactor: behaviour-preserving unless it says otherwise, and
well under 700 net lines (most are net negative). The order puts the
safety net first and the `cli.py` surgery after 23c-2 merges.

| # | PR (branch name) | Takes | Net production lines | Waits for | Gates beyond review |
|---|---|---|---|---|---|
| T2 | `audit-layering-test` | X3 | 0 (tests only, -110) | nothing | none |
| T1 | `audit-cli-test-harness` | X1 | 0 (tests only, -350) | nothing | none |
| B | `audit-crs-helpers` | F3, with the `EPSG:None` fix | about -25 | nothing | red test for the fix |
| A | `audit-lattice` | F2, F9, F10 (repository Protocol) | about -100 | B | `@perf` run: meshes byte-identical |
| F | `audit-catchment-shared` | F4, F10 (catchment types) | about -50 | nothing | none |
| C | `audit-geojson-io` | F5 | about -40 | B | `@tester` amendment if wordings move |
| D | `audit-encoders` | F7 | about -30 | C (shares `io/geojson.py`) | none |
| E | `audit-topology` | F6, X2 for `_chain_masks`/`_undirected` | about -35 | 23c-2 merged | none |
| G | `audit-cli-options` | F1, F11 | about -100 | 23c-2 merged | none |
| H | `audit-mesh-run` | F8, X2 for the rest | about -60 (about 550 moved) | G, E | `@perf`: bench tool seam and byte-identical meshes |
| tools | `audit-tools-git` | section 4 | about -20 (tools are not production; governed files need Ola) | nothing | Ola's approval per governed file |

Total: about -450 production lines, -460 test lines, and the drift points
(lattice spelling, NoData rule, CRS checks, GeoJSON `crs` rules, mask
convention) each written once.

If 23c-2 is not near merging, E, G and H can instead be folded into its
follow-up: 23c-2 should put its own run in `pieces.py`, not in `cli.py`, so
the cli split does not grow by its 250 lines.
