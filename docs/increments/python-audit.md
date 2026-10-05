# Audit: the Python layer — architecture, duplication, dead code

Status: **audit, not a design.** Written by `@architect` against master at
`12dace7` (its `src_python/` is byte-identical to `bc01cd8`). Read-only on
code. Every `file:line` citation is pinned to `12dace7`, so it stays true
while the code moves. It proposes a sequence of refactor PRs; each still runs the normal loop
(`docs/increments/README.md`): design note, red tests where behaviour changes,
code, review, and `@perf` where the diff touches what drives refine or mesh.

Accepted by Ola on 2026-10-05, with the defaults to its three questions
(section 7). Its first PR, T2, is designed in section 8. Status of T2:
`@tester`'s tests are in `2f47ebb`: 184 non-blank test lines added and 185
removed, -1 net against the design's about -35 (`git diff -U0 44fa7f5 97eea35
-- tests`, non-blank lines). The gap is what section 8 did not cost:
`test_layering.py` adds 152 against about 120, because its 13-line module
docstring and the 12-line `UPWARD` dict with its comment were not costed;
`importscan.py` adds 19 and removes 5 against about 8; reworded references
to `test_layering.py` in six other test files add 13, not costed; and 185
go against about 165. `4686456` adds 4 and removes 2, so the branch is +1
net. Code review rounds 1 and 2 (below) are recorded and their blocking
items fixed (`4686456`, this paragraph, and the round-1 record's two
citations, now pinned to `97eea35`). Next `@reviewer` round 3: run
`python3 tools/check_citations.py` and re-read the round-1 record (tests
only: no `@developer` step, no `@perf` run); then push on Ola's yes.
PR F (`audit-catchment-shared`, branch `worktree-audit-catchment`, section
11): red `84e12fa`, green `8d9e2c5` and its trim `5c6a9f9`; code review
round 1 asked for changes, made (Review, below); next `@reviewer` round 2.

Re-checked against master `44fa7f5`: `git diff --stat 12dace7 44fa7f5 --
src_python` is empty, and of the files cited below only `tools/brief.py`
changed, outside the cited lines (`git show 44fa7f5:tools/brief.py | sed -n
115,118p` is still `_git`). The `@12dace7` pins therefore also read true at
`44fa7f5`, and the findings stand unchanged; 23c-2 (`worktree-23c`) has still
not merged, so section 6's waits stand too.

Ola's ask: look at the architecture, modularization and class hierarchies,
and above all at code doing almost the same thing in different places without
sharing logic.

**Prior art.** None applies: this audits shipped code and claims nothing new.
A PR below that designs a new module carries its own prior-art section.

## Verdict in five lines

1. The duplication is real but smaller in lines than it looks: about **440
   production lines** (6 % of 7,732 code lines) and **385 test lines** would
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
three. The bug fix needs a red test first (`@tester`). It is a product fix in
its own right: the Austrian data work that first met such a CRS stays out of
git (Ola's ruling), so the red test builds its own EPSG-less tile.

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
share the polygon-part filter. Ola's ruling (section 7): both `--lakes`
behaviours stay, and each command's `--help` for `--lakes` says which file it
takes. Risk: medium; the refusal wordings differ today
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
both. The convention is then written once. Ola's ruling (section 7): `viz/`
may import this one shared module. This changes `viz/scene.py`'s import rule
(today it imports only `protocols`); after T2 that rule is `viz.scene`'s row
in `test_layering.py`, which the PR's `@tester` commit widens by `topology`.
The rule's purpose (no `_core` in `viz/`) still holds.

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
`RasterMeta`. Ola's ruling (section 7): +-infinity is not NoData, as in the
C++ core, so `target_grid` changes behaviour and PR A needs a red test for it.

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

### F12. Imports that point up a layer — 0 to 20 lines, low risk (found while designing T2)

Tabling every module's imports against section 5's layers (section 8) finds
eight edges that point upward. Each is allowed in T2's table as a named
exception, with the PR that removes it:

| Edge (importer -> imported) | Why it exists | Removed by |
|---|---|---|
| `mosaic` -> `io.cog` | `window_meta`, which is lattice arithmetic on `RasterMeta` | A (becomes a `RasterMeta` method) |
| `mosaic` -> `io.repository` | `TileFootprint`, under `TYPE_CHECKING` | A (`TileFootprint` to `io/models`) |
| `chains` -> `feature_input` | the value type `TerrainFeature` lives in a pipeline module | C (the type moves down, to `features`) |
| `gauge` -> `io.rivers`, `io.station_set` | the value types `RiverSegment` and `Lake` live in codecs | F (the types move to L0) |
| `fetch.http` -> `tin_engine` (the package root) | `installed_version`, in an `__init__` that also imports `_core` | G (`installed_version` to its own L0 module) |
| `catchment` -> `_core` | `upstream`, `accumulate`, `reduce_ring` called from a pipeline | F (a small L3 adapter beside `raster`) |
| `cli` -> `_core` | the engine and `ChainRole` | H (engine to `start_mesh`); the `ChainRole` use may stay, and H says so |

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
- **X3. Import-firewall tests scattered over 8 files — about 165 test lines,
  about 35 net after the new table (PR T2, section 8).** `tests/python/test_viz_svg.py@12dace7:536-594`, `tests/python/test_viz_scene.py@12dace7:366-405`,
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
                      catchment_batch, fetch/run, fetch/nve, fetch/plan
                      request (frozen Pydantic) in, result out; no typer,
                      no print, no path but what the request names
L3  core adapters     raster.to_core (the one raster adapter), start_mesh
                      (build_pslg -> node -> triangulate), edge_strip,
                      catchment_core (PR F), final_check, and the package
                      root (it re-exports Point2/Point3): the only
                      importers of _core
L2  io/ codecs        bytes <-> values: geotiff, cog, geopackage, gml,
                      geojson (read AND write), ply, vtk_legacy, tables
                      (csv/json/palette), mesh_index, station_set, rivers;
                      repository.py the one module that opens files;
                      fetch/http, the one module that opens a connection
L1  pure algorithms   crs (same_crs, transform_label, single_crs), mosaic,
                      target_grid, grid_domain, domain, chains, elevation,
                      outline, burn, gauge, sensitivity, reference,
                      landcover, decompose, topology (new), viz/*
L0  values            io/models (RasterMeta with node methods, IndexWindow
                      with window methods, Bounds, TileFootprint, DemTile,
                      valid_mask), features, sources, run_record, stats,
                      palettes, hydrography (PR F: RiverSegment, Station,
                      Lake)
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
safety net first and the `cli.py` surgery after 23c-2 merges. After T2, a PR
that adds, drops or moves a first-party import edits `test_layering.py`'s
table in a `@tester` commit, and deletes the section 8 exception it removes
(F12 says which).

| # | PR (branch name) | Takes | Net production lines | Waits for | Gates beyond review |
|---|---|---|---|---|---|
| T2 | `audit-layering-test` | X3 | 0 (tests only, about -35) | nothing | none; `@tester` then `@reviewer`, no `@developer` |
| T1 | `audit-cli-test-harness` | X1 | 0 (tests only, -350) | nothing | none |
| B | `audit-crs-helpers` | F3, with the `EPSG:None` fix | about -25 | nothing | red test for the fix |
| A | `audit-lattice` | F2, F9, F10 (repository Protocol), F12 (`mosaic`'s two) | about -100 | B | red test for the +-inf ruling; `@perf` run: meshes byte-identical |
| F | `audit-catchment-shared` | F4, F10 (catchment types), F12 (`gauge`'s two, `catchment` -> `_core`) | about +15 (section 11; was about -40) | nothing | red test for the lakes type; the byte probe (section 11) |
| C | `audit-geojson-io` | F5, F12 (`chains` -> `feature_input`) | about -40 | B | `@tester` amendment if wordings move, and for the two `--help` texts |
| D | `audit-encoders` | F7 | about -30 | C (shares `io/geojson.py`) | none |
| E | `audit-topology` | F6, X2 for `_chain_masks`/`_undirected` | about -35 | 23c-2 merged | none |
| G | `audit-cli-options` | F1, F11, F12 (`installed_version`) | about -95 | 23c-2 merged | none |
| H | `audit-mesh-run` | F8, X2 for the rest | about -60 (about 550 moved) | G, E | `@perf`: bench tool seam and byte-identical meshes |
| tools | `audit-tools-git` | section 4 | about -20 (tools are not production; governed files need Ola) | nothing | Ola's approval per governed file |

Total: about -385 production lines (-440 before section 11 re-measured F), -385 test lines, and the drift points
(lattice spelling, NoData rule, CRS checks, GeoJSON `crs` rules, mask
convention) each written once.

If 23c-2 is not near merging, E, G and H can instead be folded into its
follow-up: 23c-2 should put its own run in `pieces.py`, not in `cli.py`, so
the cli split does not grow by its 250 lines.

## 7. Ola's rulings (2026-10-05)

Ola accepted the audit and the defaults to its three questions:

1. **+-infinity is not NoData**, as in the C++ core
   (`include/terrain/raster/raster.hpp@12dace7:63-66`). NoData is NaN or the
   sentinel. `target_grid._valid` follows the others (F9, PR A, red test first).
2. **Both `--lakes` behaviours stay.** The two readers get names that say what
   they read (`read_lake_polygons`, `read_nve_lakes`), and each command's
   `--help` for `--lakes` says which file it takes (F5, PR C).
3. **`viz/` may import one shared pure-Python module** for the mesh-edge
   joins, `topology.py` (F6, PR E). Still no `_core` in `viz/`.

The EPSG-less CRS fix (F3) stands as a product fix on its own: the Austrian
side-work that first met such a file stays out of git, by Ola's ruling.

## 8. T2 design: `tests/python/test_layering.py`

Tests only; no production file changes, so no red step and no `@developer`.
`@tester` writes it and `@reviewer` audits it. It must pass on master.

### Where it lives, and what it reads

- `tests/python/test_layering.py`: the table and the checks, in one file, so
  the dependency map is what a reviewer reads.
- `tests/python/importscan.py` stays the one scanner, with two fixes (tests
  only). Both are bugs today; the probe for each, run with the project venv
  from `tests/python/`:
  - **(a) A package `__init__` resolves its relative imports against its
    parent.** `first_party_imports(tin_engine.io)` returns
    `{'tin_engine.ply', 'tin_engine.vtk_legacy'}`; it should return
    `tin_engine.io.ply` and `tin_engine.io.vtk_legacy`. Fix: the base package
    is `module.__name__` when the module has `__path__`, else its
    `rpartition(".")[0]`.
  - **(b) `from <package> import <submodule>` reports only the package.**
    `src_python/tin_engine/cli.py@12dace7:62` is `from tin_engine import edge_strip, final_check,
    installed_version`, and `first_party_imports(tin_engine.cli)` contains
    neither `tin_engine.edge_strip` nor `tin_engine.final_check`. Fix: for each
    imported name, report `<package>.<name>` when the package has `__path__`
    and `importlib.util.find_spec("<package>.<name>")` finds it; report the
    package itself only when some name is not a submodule (here
    `installed_version`, so `cli` also imports the root `tin_engine`).
  - Unchanged: `TYPE_CHECKING` imports count (they are where the
    `io.repository` <-> `mosaic` cycle lives); a mention in prose is not an
    import.
- The scanner imports each module (`inspect.getsource`), and importing any
  `tin_engine` module runs the package root, which imports `_core`; so the
  suite needs the built extension, as the eight tests it replaces already do.

### The table

One dict, module (short name, `tin_engine.` dropped; the package root is
`tin_engine`, the extension `_core`) to (layer, the exact set of first-party
modules it imports). Rows as today's code has them, read with the fixed
scanner: equality, not a subset, so the table is the map and a new edge is a
visible table edit. Layers are section 5's:

- L0: `io.models`, `features`, `sources`, `run_record`, `stats`, `palettes`
- L1: `crs`, `mosaic`, `target_grid`, `grid_domain`, `domain`, `chains`,
  `elevation`, `outline`, `burn`, `gauge`, `sensitivity`, `reference`,
  `landcover`, `decompose`, `viz`, `viz.fixtures`, `viz.protocols`,
  `viz.scene`, `viz.style`, `viz.svg`
- L2: `io`, `io.cog`, `io.geojson`, `io.geopackage`, `io.geotiff`, `io.gml`,
  `io.mesh_index`, `io.ply`, `io.repository`, `io.rivers`, `io.station_set`,
  `io.vtk_legacy`, `fetch.http`
- L3: `raster`, `edge_strip`, `final_check`, `tin_engine`, `_core`
- L4: `dem_input`, `feature_input`, `catchment`, `catchment_batch`, `fetch`,
  `fetch.plan`, `fetch.run`, `fetch.nve`
- L5: `cli`

That is all 52 modules plus `_core`. A second dict, `UPWARD`, holds F12's
eight edges, each with the PR that removes it as its value.

The rows that carry the eight old tests' rules, as they are today (the
reviewer checks these against the old assertions): `viz.svg` {`viz.scene`,
`viz.style`}; `viz.scene` {`viz.protocols`}; `viz.style`, `viz.fixtures`,
`viz.protocols`, `features`, `landcover`, `outline`, `io.models` {}; and
`io.geotiff` {`io.models`}; `io.ply`, `io.vtk_legacy` {`features`}. Each is
equal to or stricter than the test it replaces.

### The checks

1. **The table names every module.** The `*.py` files under
   `Path(tin_engine.__file__).parent`, as dotted names, plus `_core`, equal the
   table's keys. A new module without a row fails, and so does a row for a
   module that is gone.
2. **Each module imports exactly its row** (parametrised by module, `_core`
   excluded): `first_party_imports(import_module(...))`, shortened, equals the
   row. The failure names the extra and the missing edges.
3. **Every edge goes down or sideways** (the table only): layer of the
   imported <= layer of the importer, and an edge into `_core` comes from L3,
   unless the edge is in `UPWARD`.
4. **No `UPWARD` entry is stale**: each is an edge of the table that breaks
   the rule of check 3, so the PR that removes the edge, or makes it legal,
   must delete its exception.
5. **No module imports by name**: no call to `importlib.import_module` or
   `__import__` in any module's AST (none today; the package's `importlib`
   uses are `metadata`, `util.find_spec` and `resources`). This replaces the
   source-text grep `tests/python/test_features.py@12dace7:606-608`, whose
   reason was the import `ast` cannot see, and it now covers every module
   rather than one, without pinning prose.

Before committing, `@tester` shows each of checks 1, 2, 3 and 5 can fail by one
plant each in a scratch copy (for example `import tin_engine._core` added to
`viz/style.py` fails check 2), restored afterwards, and names the plants in the
commit message. No mutation round beyond that.

### What it replaces

Delete, in the same commit:

| File (lines at `12dace7`) | What goes |
|---|---|
| `tests/python/test_viz_svg.py@12dace7:536-592` | `TestModuleIsolation`, including the `ChainRole` text check on `cli.py` (prose; the mapping's totality is `test_cli_draw.py`'s) |
| `tests/python/test_viz_scene.py@12dace7:366-403` | `TestModuleIsolation` |
| `tests/python/test_features.py@12dace7:575-608` | both firewall tests, the `_core` source-text grep among them, and `FEATURES_SOURCE` |
| `tests/python/test_io_geotiff.py@12dace7:1284-1310` | `_imported_modules` and `test_module_never_imports_core` |
| `tests/python/test_landcover.py@12dace7:449-453` | `TestPurity` |
| `tests/python/test_outline.py@12dace7:136-139` | `test_the_tracer_imports_no_core` |
| `tests/python/test_io_ply.py@12dace7:428-439` | `test_the_only_first_party_import_is_the_vocabulary` (the rest of `TestPurity` stays) |
| `tests/python/test_io_vtk_legacy.py@12dace7:487-492` | `test_the_only_first_party_import_is_the_vocabulary` |

and any import or path constant (`ast`, `importscan`, `SCENE_PATH`, `VIZ`,
`REPO_ROOT`) left unused; `ruff check` finds them.

Expected test-line delta: about 165 non-blank lines go, about 130 come
(table about 75, checks about 45, `importscan` about 8): about -35 net. The
gain is the map, not the lines.

### Merge-order hazard

23c-2 (`worktree-23c`) adds `pieces.py` and new `cli` imports. Whichever of
T2 and 23c-2 merges second must add those rows; the merge queue tests each on
top of the other, so the second one goes red there rather than on master.

### Left for a later PR

- `src_python/tin_engine/features.py@2f47ebb:11` still names
  `test_viz_svg.py::TestModuleIsolation`, which T2 deleted; the rule now
  lives in `tests/python/test_layering.py`. A docstring fix for PR C
  (`audit-geojson-io`) or PR G (`audit-cli-options`), whichever touches
  `features.py` first.
- `tools/scratch_copy.py@2f47ebb:19` cites line 1271 of
  `tests/python/test_io_geotiff.py` for the child process that replaces
  `PYTHONPATH`. T2 removed one line above it, so line 1271 is still inside
  the same `subprocess.run` call (`tests/python/test_io_geotiff.py@2f47ebb:1269-1276`)
  but one line lower in it. The quotation still holds; re-cite it in
  whichever PR next edits `tools/`.

### Lesson

A design that deletes test lines lists the line citations into those files
that the deletion breaks or moves (`python3 tools/check_citations.py` on a
scratch copy with the lines removed), and says which get pinned. T2's design
did not, and its tests commit left nine broken citations in dated review
records.

## 11. PR F design: `audit-catchment-shared` (F4; F10's catchment types; F12's three edges)

Branch `worktree-audit-catchment`, from T2's head `b63132e`. The section
number assumes T1's design (section 9 on `worktree-audit-t1`) and PR B's
(section 9 on `worktree-audit-crs`, 10 once T1 lands) come first; the PR
that lands last renumbers. `src_python/` at `b63132e` is byte-identical to
`44fa7f5` and to `12dace7` (`git diff --stat 12dace7 b63132e -- src_python`
is empty), so findings F4, F10 and F12 stand as written, and the citations
here are pinned to `44fa7f5`, which is on master.

What it does, in one line each:

1. `rasputin catchment --rivers` and `run_batch` build a station's catchment
   request with one function, `catchment_batch.seed_for`, and report the
   placement and the gauge with one method each (F4).
2. The `Any`s in the catchment code become the types that exist (F10).
3. `gauge` stops importing two codecs, and `catchment` stops importing
   `_core` (F12): the three value types move to a new layer-0 module
   (`hydrography.py`), and the three `_core` calls move to a new layer-3
   adapter (`catchment_core.py`).

**No `@perf` run.** The acceptance rule covers refine and mesh code and what
drives it (`docs/increments/README.md`, "Acceptance"). Catchment
delineation does not drive refine or mesh: it writes a polygon file, which
a later, separate `mesh --domain` run reads as its input. The PR touches no
C++, no `cli.mesh` code path, and nothing `tools/bench.py` runs
(`grep -n -i "catchment\|station" tools/bench.py` is empty). Its safety net
is the byte comparison below.

### Prior art: legacy and literature

*Literature.* None applies. This moves code that exists; it adds no method
and claims nothing new.

*Legacy.* Nothing. The legacy tree has no catchment, gauge or station code:

```
$ git grep -c -i -E "gauge|catchment|station" legacy-archive -- legacy
(no output, exit 1)
```

The tag holds 30 files under `legacy/` (`git ls-tree -r --name-only
legacy-archive -- legacy | wc -l`), and the same grep for `raster` lists
`legacy/bindings.cpp`, `legacy/rasputin/application.py` and
`legacy/rasputin/geo_tiff_reader.py`, so the empty result is not a grep that
could not match.

### The shared placement (F4)

Today, at `44fa7f5`:

- `run_batch` places the gauge, picks its lake seed, and builds the
  `CatchmentRequest` (`src_python/tin_engine/catchment_batch.py@44fa7f5:177-209`); the single command does the
  same for the river path in `_placed` and in `catchment`'s body
  (`src_python/tin_engine/cli.py@44fa7f5:1709-1728`, `src_python/tin_engine/cli.py@44fa7f5:1888-1901`). The refusal "no mapped river line within
  N m of the station" is written in both
  (`src_python/tin_engine/catchment_batch.py@44fa7f5:207-209`, `src_python/tin_engine/cli.py@44fa7f5:1724-1726`).
- The report fields of the placement (its dump less `reach` and `position`,
  plus `uncertainty_m`) and of the gauge (`GaugeResult`'s slots less four,
  and the sensitivity) are built twice:
  `src_python/tin_engine/catchment_batch.py@44fa7f5:132-137, 184-185` and `src_python/tin_engine/cli.py@44fa7f5:1782-1785`.

**The two gauge copies are not the same, and the difference is kept.** The
CLI's copy is `{**keep, "uncertainty_m": ..., **burn, **asdict(s),
"causes": causes}`: it keeps the sensitivity's `well_posed`, and its
`causes` key sits where `Sensitivity` declares it (before `well_posed`),
holding the gauge's joined causes. The batch's copy drops `causes` and
`well_posed` from the sensitivity and appends `causes` last; it has no
`well_posed` because `StationResult` has no such column. The catchment
file's bytes depend on the key order (the probe below changes hash when
only `well_posed` and `causes` swap places), so the shared method returns
the CLI's dict, in the CLI's order, and the batch drops `well_posed`. The
batch's order does not matter: it goes into `StationResult(**fields)`.

The shape:

- **`gauge.Placement.report_fields(self) -> dict[str, Any]`**:
  `{**self.model_dump(exclude={"reach", "position"}), "uncertainty_m":
  self.reach.uncertainty}`. Docstring: the placement's columns, as the
  catchment file and `StationResult` name them.
- **`catchment.GaugeResult.report_fields(self) -> dict[str, Any]`**:
  the slots less `node`, `chain`, `sensitivity` and `causes`, then
  `asdict(self.sensitivity)` with its `causes` replaced in place by
  `self.causes`. The docstring says the order is the catchment file's and
  that `causes` is the joined one. `dict[str, Any]`, not `object`: the batch
  spreads it into `StationResult(**fields)`, which mypy checks per keyword.
- **`catchment_batch.NoRiverLine(CatchmentError)`**: no mapped line within
  the map radius. Its one message is today's, unchanged:
  `f"no mapped river line within {map_radius:g} m of the station"`.
- **`catchment_batch.seed_for`**:

  ```python
  def seed_for(
      gauge: Gauge,
      segments: Sequence[RiverSegment],
      crs: str,
      request: BatchRequest,
      lakes: Sequence[Lake] | None = None,
  ) -> tuple[Placement | None, LakeSeed | None, CatchmentRequest]:
  ```

  `gauge` and `segments` are in `crs`, the river file's CRS. It places the
  gauge (`place`, with `request.map_radius` and `request.reach_up`), picks
  the lake seed when `lakes` is given (`lake_seed`), and returns the
  placement, the seed and the request: seeded by the lake when there is a
  seed (`seed=seed.point`, `lakes` the seed's polygons, `lakes_crs=crs`),
  else by the reach (`seed=(gauge.x, gauge.y)`, `reach=placement.reach`),
  with `seed_crs=crs` and `request.outline_tolerance` either way, exactly
  the two requests `src_python/tin_engine/catchment_batch.py@44fa7f5:190-204` builds. With neither, it raises
  `NoRiverLine`. It reads no file and checks no CRS: `check_reach_crs`
  stays with the callers, which have the repository.
  `BatchRequest` carries the three numbers for both callers: the single
  command is a batch of one.

The callers:

- `run_batch`: one `try` holds `seed_for`, the placement's and the seed's
  fields, and the `delineate` call; `except NoRiverLine` sets
  `refusal_cause` "no_river" and its message, before the existing
  `MixedGridRefusal` and `CatchmentError` clauses (`NoRiverLine` is a
  `CatchmentError`, so it must come first). The placement's fields are
  recorded before `delineate` runs, so a refusal by `delineate` keeps them,
  as today. The gauge's fields are `report_fields()` less `well_posed`.
  `delineate` is still called by its name in `catchment_batch`'s namespace:
  `tests/python/test_cli_station_catchments.py@44fa7f5:328, 367` patch it
  there. `_gauge_fields` goes.
- `cli.catchment` with `--rivers`: `_placed` keeps reading the river file,
  `_reach_crs`, and moving the seed into the river file's CRS; it builds a
  `BatchRequest` from `--map-radius`, `--reach-up` (the `or 500.0` and
  `or 1000.0` defaults stay; F1 and PR G own them) and
  `--outline-tolerance`, calls `seed_for` with no lakes, and turns
  `NoRiverLine` into `BadParameter(str(exc), param_hint="--rivers")`, the
  wording and hint of today. It returns the placement, the line's river
  name, and the request, which `catchment` passes to `delineate`. Without
  `--rivers`, `catchment` builds its request as today. `_placement_report`
  returns `{**p.report_fields(), **gauge.report_fields()}`.

**The two checks that run twice stay twice, on purpose.** `check_reach_crs`
(`src_python/tin_engine/cli.py@44fa7f5:2146` then `src_python/tin_engine/catchment_batch.py@44fa7f5:157`) and the unknown `--only` check
(`src_python/tin_engine/cli.py@44fa7f5:2126-2130`, `src_python/tin_engine/catchment_batch.py@44fa7f5:154-156`): `run_batch` keeps both as
preconditions of a library call, since a GUI or API caller has no CLI in
front of it; the CLI keeps both as usage errors, which name the flag and
refuse before `--out-dir` is created (the suite pins that order,
`tests/python/test_cli_station_catchments.py@44fa7f5:289`). The rule itself
is written once (`check_reach_crs`); the `--only` wording is written twice
and both suites pin it.

**`Gauge` stays a separate type.** It is what `place` and `lake_seed` need
(a point, a watercourse number, a river name), and the single command has a
point but no station. After this PR the field-by-field copy from `Station`
is in one place, `run_batch`.

### Types (F10, the catchment part)

| Where (at `44fa7f5`) | Today | After |
|---|---|---|
| `src_python/tin_engine/catchment.py@44fa7f5:89` | `lakes: tuple[Any, ...] \| None` | `tuple[BaseGeometry, ...] \| None` (shapely's base class: `_lake` calls `get_parts` on any geometry, then requires a `Polygon`) |
| `src_python/tin_engine/catchment.py@44fa7f5:178` | `Flood = Callable[[Any], tuple[UpstreamOutcome, Any]]` | `type Flood[T] = Callable[[DemTile], tuple[UpstreamOutcome, T]]`; `_grow[T]` returns `T` where it returned `Any` |
| `src_python/tin_engine/catchment.py@44fa7f5:210, 270` | `flood(tile: Any)` | `tile: DemTile` (what `assemble(...).tile` is) |
| `src_python/tin_engine/catchment.py@44fa7f5:226, 285, 362` | `footprints: Any` | `Sequence[TileFootprint]` (what `plan_mosaic` takes) |
| `src_python/tin_engine/catchment.py@44fa7f5:266, 299` | `pick: Callable[[Any], int]`, `at_u(path: Any)` | `GaugePath` (`burn.py`, already imported from) |
| `src_python/tin_engine/catchment.py@44fa7f5:371, 452, 521` | seed and mask `Any` | `npt.NDArray[np.uint8]`, through one alias `Mask` |
| `src_python/tin_engine/cli.py@44fa7f5:1711, 1731` | `repository: Any` | `DemRepository` |

`catchment.py` then imports nothing from `typing` but `Any` (for
`GaugeResult.report_fields`), `Literal` and `Self`.
Left as `Any`, on purpose: `run_batch`'s `fields: dict[str, Any]` (spread
into `StationResult(**fields)`) and the two `report_fields` returns, for the
same reason; and `sensitivity.assess(path: Any)`, which is not catchment
code and would add a `sensitivity -> burn` edge.

### The upward edges (F12)

- **`gauge` -> `io.rivers`, `io.station_set`.** New module
  `src_python/tin_engine/hydrography.py`, layer 0, no first-party import:
  `RiverSegment` (from `io/rivers.py`), `Station` and `Lake` (from
  `io/station_set.py`), moved unchanged (same fields, validators, frozen
  configs, docstrings). Every first-party importer imports them from
  `hydrography`: `gauge`, `catchment_batch`, `cli`, and the two codecs. The
  codecs drop them from their `__all__`. Three test imports move with them
  (`tests/python/test_gauge.py@44fa7f5:46, 539`, `tests/python/test_catchment_batch.py@44fa7f5:500`), in `@tester`'s
  commit. Why a new module rather than `io/models.py`: that file holds the
  raster's values, and PR A rewrites it (`Bounds`, `TileFootprint`, the
  lattice methods); and why not `gauge.py`: a codec would then import an
  algorithm module for a value type.
- **`catchment` -> `_core`.** New module
  `src_python/tin_engine/catchment_core.py`, layer 3, the catchment's
  `_core` calls, as `edge_strip.py` is the edge strip's:
  `upstream(tile: DemTile, seed: npt.NDArray[np.uint8]) -> UpstreamOutcome`
  and `accumulate(tile: DemTile) -> AccumulateOutcome`, each the core call
  on `raster.to_core(tile)`; and `reduce_ring`, `ReduceStatus`,
  `UpstreamOutcome` re-exported for `catchment`'s use and annotations
  (`__all__`, which mypy's strict mode needs for a re-export). `catchment`
  then imports neither `_core` nor `raster`, and never holds a core raster
  view. `raster.py` stays the one module that builds one
  (`project_structure.md`, "Boundary contract"). Not a pure re-export of
  `_core`'s `upstream`: that would pass the table and leave `catchment`
  calling the core's interface directly.

### `tests/python/test_layering.py` (a `@tester` commit)

Rows, after:

- layer 0: new `"hydrography": ""`.
- layer 3: new `"catchment_core": "_core io.models raster"`.
- `"gauge": "hydrography"` (was `io.rivers io.station_set`).
- `"io.rivers": "hydrography io.station_set"` (it still imports the
  station set's private helpers).
- `"io.station_set": "crs hydrography io.repository"`.
- `"catchment": "burn catchment_core crs gauge io.models io.repository
  mosaic outline sensitivity"` (loses `_core` and `raster`).
- `"catchment_batch": "catchment crs gauge hydrography io.repository
  reference"` (loses `io.rivers` and `io.station_set`).
- `"cli"`: gains `hydrography`; nothing else changes.

`UPWARD` loses three entries: `("gauge", "io.rivers")`, `("gauge",
"io.station_set")` and `("catchment", "_core")`. Check 4
(`test_no_upward_exception_is_stale`) fails until they are deleted, so the
table commit and the code must land together in the PR; the table commit is
red against `b63132e`'s code, which is its red step for the import changes.

### Red tests (`@tester`, one commit with the table, before any code)

Behaviour changes in one place only: `CatchmentRequest` refuses lakes that
are not shapely geometries, at construction. Probe (project venv, at
`44fa7f5`): `CatchmentRequest(seed=(0.0, 0.0), seed_crs="EPSG:25833",
lakes=("not a polygon",), lakes_crs="EPSG:25833")` is accepted today; a
model with `tuple[BaseGeometry, ...] | None` and
`arbitrary_types_allowed` raises `ValidationError` (a `ValueError`) on it,
keeps a `Polygon` and a `MultiPolygon` as they are (the same objects), and
turns a list into a tuple, as today. No caller passes anything else: the
CLI's lakes come from `feature_input.read_lakes`
(`src_python/tin_engine/feature_input.py@44fa7f5:429-446`, shapely geometries) and the batch's from `Lake.polygon`.

1. **`tests/python/test_catchment.py`**: a request whose `lakes` holds a
   string raises `pydantic.ValidationError` naming `lakes`. Red today.
2. **`tests/python/test_layering.py`**: the rows and `UPWARD` above. Red
   today (the new modules have no file; three exceptions are not stale yet).
3. The three test imports re-pointed to `tin_engine.hydrography`. Red today
   (no such module).

Nothing else is a new test: the shared function changes no output, and the
suites that pin today's outputs (`test_cli_catchment.py`,
`test_catchment_batch.py`, `test_cli_station_catchments.py`) stay as they
are. No mutation round.

### The safety net: byte-identical outputs

`docs/increments/python-audit-probes/catchment_bytes.py` builds the suites'
own fixtures (`batch_fixtures`, `gauge_fixtures`) and runs eight commands:
`station-catchments` with references, with the lake file (one station seeded
inside its lake), with a lake-line river file and its lake (one station
seeded by the `lake_line` rule), and on the mixed-grid DEM (one `mixed_grid`
refusal); and `catchment` with `--rivers`, with `--rivers` and a 15 m radius
(the no-river refusal), plain (the outlet node), and with `--lakes`. Each
run's rows cover the river seed and the `no_river` refusal too. It prints a
SHA-256 for each written file and for the command's output, with timings
masked (the `seconds` column, and every `<number> s` on stderr) and the
temporary path replaced. Its first line names the `tin_engine` it loaded.

Run at `b63132e`'s code (the main checkout's venv, whose `src_python/` is
`44fa7f5`'s, byte-identical), twice: both outputs are identical (33 lines:
the package's path, then 32 hashes).
It can fail: with two plants applied at run time (the CLI's `well_posed` and
`causes` keys swapped, nothing else; the batch's `lowered_nodes` plus one),
16 of the 32 hashes change, the single command's catchment file among them,
and the lake-seeded files, which carry no gauge, do not. Before the push,
`@reviewer` runs it on the branch with the branch's own venv and compares
it with a run at `b63132e`: every line after the first is equal.

```bash
PYTHONPATH=tests/python .venv/bin/python docs/increments/python-audit-probes/catchment_bytes.py
```

### Net production lines

| File | Added | Removed | Net | Measured (`b63132e..5c6a9f9`) |
|---|---|---|---|---|
| `hydrography.py` (new) | about 30 | 0 | about +30 | +30 |
| `io/rivers.py`, `io/station_set.py` | about 2 | about 29 | about -27 | -27 (3 added, 30 removed) |
| `gauge.py` | about 3 | about 1 | about +2 | +2 |
| `catchment.py` | about 6 | 0 | about +6 | +14 (39 added, 25 removed) |
| `catchment_core.py` (new) | about 13 | 0 | about +13 | +12 |
| `catchment_batch.py` | about 43 | about 48 | about -5 | -1 (49 added, 50 removed) |
| `cli.py` | about 18 | about 22 | about -4 | -1 (33 added, 34 removed) |
| **Total** | | | **about +15** | **+29** (171 added, 142 removed) |

**Measured, after the green commits: +29, not about +15**, by `python3
tools/count_loc.py b63132e 5c6a9f9`; under the +40 that would be a finding.
Of the 29, 16 are import lines (top-level `import` and `from` statements,
counted by line at both revisions): `catchment.py` +6, `catchment_core.py`
+7, `hydrography.py` +5, `cli.py` +2, `catchment_batch.py` -2, `gauge.py`
and `io/station_set.py` -1 each. The three files off their estimate:

- `catchment.py`, +14 against about +6: the `catchment_core` import takes
  seven lines where the `_core` one took one (on one line it is 102
  characters, over the formatter's 100), so imports are +6;
  `GaugeResult.report_fields` is 4 counted lines, which the table did not
  cost in this file; and the types are 4 more than costed (two more
  aliases, `Mask` and `Footprints`, and `_burnt_flood`'s signature over
  three lines).
- `catchment_batch.py`, -1 against about -5: imports -2; the batch's copy
  of the gauge's fields drops `well_posed` in its own statement, and the
  `try` round `seed_for` holds the placement's and the seed's fields.
- `cli.py`, -1 against about -4: imports +2 (`hydrography`, and
  `DemRepository` for the retyped `repository`), and the `BatchRequest`
  the single command builds takes five lines.

**Line count is the proxy, complexity the measure.** Ola, 2026-10-05,
answering the main session's question D4: "The LOC is basically a proxy
for complexity. Imports add very little. So it sounds like a would be
right here." By that measure the +29 is about +13 outside import
lines, for two placement copies made one, typed catchment signatures, and
three upward imports removed.

**Not the audit's about -50 (section 2) or -40 (section 6).** The two
placement copies share about 25 lines a side, and the shared function and
its exception cost about 21, so F4 saves about 10. F10 is about 0. F12 costs
about +16: two new modules, of which the adapter's 13 lines are the price of
the rule that only layer 3 imports `_core`, and the value types' move about
+3. Section 6's row F says so. The reviewer counts with `python3
tools/count_loc.py b63132e <head>`; a result above +40 is a finding.

`project_structure.md` (not counted) gains rows for `hydrography.py` and
`catchment_core.py`, and the rows for `catchment.py` (its `_core.upstream`
and `_core.reduce_ring` become `catchment_core`'s), `catchment_batch.py`
(`seed_for`), `io/rivers.py` and `io/station_set.py` change; section 5's
picture gains `hydrography` in L0 and `catchment_core` in L3. The green
commits left both out; code review round 1 asked for them, and the commit
that records that round makes them. `io/station_set.py`'s row needed no
change: it names what the readers return, which still holds.

### Citations this PR moves, pinned now

Every live unpinned line citation into a file this PR changes, found with
`git grep -nE "(catchment|catchment_batch|gauge|cli|raster|rivers|station_set|test_layering|test_catchment|test_catchment_batch|test_gauge)\.py:[0-9]+|project_structure\.md:[0-9]+" -- '*.md' '*.py'`,
less the pinned ones. Pinned to `44fa7f5` in this design commit, each quotation
re-read there:

- `docs/benchmarks/2026-10-05/nve-hrd/README.md:53` and
  `docs/increments/29-nve-reference-catchments.md`'s Sagafoss line
  (`catchment.py:259`, the `_plan` call in `_grow`);
- `docs/increments/15c-geographic-dem.md`'s `catchment.py:194-195`;
- `docs/increments/15f-edge-strip.md`'s `cli.py:1612`;
- `docs/increments/25-plain-output.md`'s `cli.py:1966-1969`, `:1040`,
  `:961`, `:1113-1125`, `:1145-1147` and `:835`;
- `docs/increments/27-node-sampling.md`'s `cli.py:845`.

All but `:835` and `:845` are the same edit, character for character, as PR B's
(`worktree-audit-crs`), so the two branches merge them without a conflict.
Left as written: `15f-edge-strip.md`'s `cli.py:1417` (it says it describes
`17c2d14` and is left as history), and the dated review records that cite
`cli.py` or `project_structure.md` lines; `check_citations.py` lists them as
at risk once the code moves, and they record what was true then.

### Overlap with PR B and T1

Whichever of PR F and PR B lands second merges master into its branch (a
merge, not a rebase: a rebase rewrites history and needs Ola's yes), then
reruns the whole Python suite, `check_citations.py` and the probe. The files
both change:

- `src_python/tin_engine/catchment.py`: B changes the `crs` import
  (`src_python/tin_engine/catchment.py@44fa7f5:39`), `check_reach_crs`'s test (`:185`) and `delineate`'s
  CRS checks (`:194-197`, `:202`); F changes the imports at `:29`, `:37`,
  `:53`, and lines from `:89` down, none of them B's. Line 38 is the only
  unchanged line between B's `:39` and F's `:37`.
- `src_python/tin_engine/cli.py`: B changes the `crs` import
  (`src_python/tin_engine/cli.py@44fa7f5:86`) and `:932-956`, `:2045`; F changes the `catchment_batch`
  import (`:84`), `:108`, `:114-115`, adds two imports, and `:1709-1735`,
  `:1766-1785`, `:1888-1902`. One unchanged line (`:85`) between `:84` and
  `:86`.
- `tests/python/test_catchment_batch.py`: B adds two tests after line 476; F
  re-points the import at line 500.
- `docs/increments/python-audit.md`: both add a section and edit the status
  paragraph, a certain conflict: keep both sections and both statuses.
- The nine citation pins B also makes: identical, so no conflict.

B's design names no site in `catchment_batch.py`, `gauge.py` or
`io/rivers.py` (its table, section 9 on `worktree-audit-crs`); the overlap
is the four files above. T1 (`worktree-audit-t1`) changes
`tests/python/test_cli_catchment.py` and `test_cli_station_catchments.py`,
which this PR does not touch, and `python-audit.md`, as above. PR A, later,
moves `TileFootprint` to `io/models.py` and rewrites the lattice arithmetic
in `_seed_mask`, `_outline`, `_extent` and `run_batch`'s boxes: it changes
the `TileFootprint` import this PR adds to `catchment.py`, and lines next
to the signatures this PR retypes.

## Review

**T2 (`audit-layering-test`), code review, round 1, 2026-10-05.** Range `44fa7f5..97eea35` (e86b86d audit, 8190438 rulings and T2 design, 2f47ebb tests, 97eea35 citation pins). Verdict: CHANGES REQUESTED. LOC: 0 production lines (`count_loc.py`); test lines +184 -185, -1 net against about -35. Not pushed; no CI. Whole Python suite 5132 passed, 17 skipped; ruff, format, mypy and check_citations clean. An independent AST resolver agrees with the table for all 52 modules. Seven planted breaks in a scratch copy each failed only their own check: a row for a missing module, a deferred import, a relative import in a package `__init__`, `_core` imported from layer 4, two stale exceptions, and `importlib.import_module`. Check 4's stricter reading is sound. Every deleted firewall assertion is carried by a row that is equal or stricter. All 16 new pins quote what their records say. Blocking: (1) `@tester`: the `# fmt: off` comment at `tests/python/test_layering.py@97eea35:27-28` describes the formatter wrongly; (2) `@architect`: the -1 net explanation at `docs/increments/python-audit.md@97eea35:12-15` names the wrong cause. Suggestions: check 4's wording in section 8; `testing.md@97eea35:183-185`'s "no compiled extension" claim, false before this branch.

**T2 (`audit-layering-test`), code review, round 2, 2026-10-05.** Range `97eea35..90cba64` (4686456 fmt-off comment, 90cba64 round 1 recorded and net explanation). Verdict: CHANGES REQUESTED. LOC: 0 production lines; branch test lines +186 −185, +1 net. Not pushed; no CI. Both round-1 blocking items fixed and true: `ruff format --diff` on a copy without the markers does what the new comment says (joined rows 98 and 99 characters, limit 100), and the per-file counts in the net explanation match `git diff -U0`. Check 4's wording matches the test. ruff, format, mypy and check_citations clean; test_layering 56 passed. Blocking: (1) `@architect`: the round-1 record's citations `tests/python/test_layering.py@97eea35:27-28` and `docs/increments/python-audit.md@97eea35:12-15` (pinned at recording) now resolve to the fixed text unless pinned. Suggestion: say "adds 4 and removes 2" at line 19.

**PR F (`audit-catchment-shared`), code review, round 1, 2026-10-05.** PR F code review r1 (`b63132e..5c6a9f9`): CHANGES REQUESTED — `project_structure.md` rows (`project_structure.md@5c6a9f9:150-153` still named `_core.upstream` and `_core.reduce_ring`, `:251` listed `RiverSegment` under `rivers.py`, no rows for `hydrography.py` and `catchment_core.py`) and section 5's picture, the status line (`docs/increments/python-audit.md@5c6a9f9:25-27`), and section 11's net lines against the measured +29; fixed by `@architect` in the commit that records this round.

Non-blocking, for a later `@tester` and `@developer` pass: tests still reach `RiverSegment`, `Station` and `Lake` through the codec modules rather than `tin_engine.hydrography` (`tests/python/test_station_set.py@5c6a9f9:51, 245`, `tests/python/test_rivers.py@5c6a9f9:89`); `catchment.py`'s module docstring still names `_core.upstream`, `_core.reduce_ring` and `_core.accumulate` where the calls now go through `catchment_core`.
