# Audit: the Python layer — architecture, duplication, dead code

Status: **audit, not a design.** Written by `@architect` against master at
`12dace7` (its `src_python/` is byte-identical to `bc01cd8`). Read-only on
code. Every `file:line` citation is pinned to `12dace7`, so it stays true
while the code moves. It proposes a sequence of refactor PRs; each still runs the normal loop
(`docs/increments/README.md`): design note, red tests where behaviour changes,
code, review, and `@perf` where the diff touches what drives refine or mesh.

Accepted by Ola on 2026-10-05, with the defaults to its three questions
(section 7). Its first PR, T2, is designed in section 8; the second, T1,
in section 10 (designed on T2's head `b63132e`; it was section 9 until PR
B's section 9 was merged in beside it). T2 merged as PR #188 and T1 as
PR #189. History of T2:
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
Since then T2 merged as #188 and T1 as #189.
PR F (`audit-catchment-shared`, branch `worktree-audit-catchment`, section
11): red `84e12fa`, green `8d9e2c5` and its trim `5c6a9f9`; code review
round 1 asked for changes, made in `e2baa5f`; round 2 approved `e2baa5f`
(Review, below). Master is merged into the branch twice: the first merge
(`e4f7a42`, test side done by `@tester`) was reviewed (Review, below); the
second (`8c2bc9a`) takes master with #188 to #192; `@reviewer` reviewed
and approved it (Review, below). F merged into master as #196.

PR B (`audit-crs-helpers`, F3) is designed in section 9, on branch
`worktree-audit-crs` after T2. `@tester`'s red commit `29aff00` has red
tests 1-7 (43 failing, each for its own reason); its step found three more
`==` comparisons, ruled in section 9 ("After the red step"). Red tests
8-10 are in `7dddda8`, which also showed rule 2 of `same_crs` wrong (a
west-pointing UTM 33 counted as EPSG:25833); rule 2 now also checks the
axes and the prime meridian (section 9, "After red tests 8-10"). `5074ef6`
adds the `+pm=paris` pair and found a thirteenth site (`domain.py:103`); a
final sweep over `src_python/` and `tools/` (section 9, "The final sweep")
finds no fourteenth. Red tests 11 and 12 are in `20bcf2c`; `@developer`'s
first green commit is `5197f9a`. Code review round 1 found a third false "the
same" (a `+lon_0` held only in PROJ's remark), so the rule is re-ruled from
what the transform does, not from CRS attributes (section 9, "After code
review round 1"). `@tester`'s re-spelt fixtures and red rows are in
`13b3e6d`, with Ola's ruling D15 b (a PROJ string naming no datum is not its
EPSG code, and the refusal says which code to write) folded into section 9
("After the round-1 red step"). `@developer`'s green commit for them is
`65cd528`. Code review round 2 asked for prose only (this paragraph, the
net lines, section 6's row B, and the hint's scope), fixed in section 9 and
the review record. Prose correction after approval: section 9's question 1
no longer claims B's rule accepts the Austrian Lambert built from its GeoKeys.
Pushed as PR #192 at `618328b`; B merged into master as #192. T1 (PR #189) then merged
into master and conflicted with it, so master is merged into the branch;
the merge renumbers T1's design to section 10.

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
eight edges that point upward; seven remain exceptions since #191 resolved
`fetch.http` -> `tin_engine` (its row below). Each is allowed in T2's table as a named
exception, with the PR that removes it:

| Edge (importer -> imported) | Why it exists | Removed by |
|---|---|---|
| `mosaic` -> `io.cog` | `window_meta`, which is lattice arithmetic on `RasterMeta` | A (becomes a `RasterMeta` method) |
| `mosaic` -> `io.repository` | `TileFootprint`, under `TYPE_CHECKING` | A (`TileFootprint` to `io/models`) |
| `chains` -> `feature_input` | the value type `TerrainFeature` lives in a pipeline module | C (the type moves down, to `features`) |
| `gauge` -> `io.rivers`, `io.station_set` | the value types `RiverSegment` and `Lake` live in codecs | F (the types move to L0) |
| `fetch.http` -> `tin_engine` (the package root) | `installed_version`, in an `__init__` that also imports `_core` | #191 (`worktree-cpp-dead`, the C++ dead-code removal): the root no longer imports `_core`, so it sits in L0 and the edge points down |
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
                      catchment_core (PR F), final_check: the only
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
                      Lake), and the package root (installed_version)
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
| T1 | `audit-cli-test-harness` | X1 (section 10) | 0 (tests only, about -100; X1's -350 is corrected there) | T2 | none; `@tester` then `@reviewer`, no `@developer` |
| B | `audit-crs-helpers` | F3, with the `EPSG:None` fix | +22, measured at `65cd528` (section 9; first estimated -25) | T2 | red tests for the rule and the fix |
| A | `audit-lattice` | F2, F9, F10 (repository Protocol), F12 (`mosaic`'s two) | about -100 | B | red test for the +-inf ruling; `@perf` run: meshes byte-identical |
| F | `audit-catchment-shared` | F4, F10 (catchment types), F12 (`gauge`'s two, `catchment` -> `_core`) | about +15 (section 11; was about -40) | nothing | red test for the lakes type; the byte probe (section 11) |
| C | `audit-geojson-io` | F5, F12 (`chains` -> `feature_input`) | about -40 | B | `@tester` amendment if wordings move, and for the two `--help` texts |
| D | `audit-encoders` | F7 | about -30 | C (shares `io/geojson.py`) | none |
| E | `audit-topology` | F6, X2 for `_chain_masks`/`_undirected` | about -35 | 23c-2 merged | none |
| G | `audit-cli-options` | F1, F11 | about -95 | 23c-2 merged | none |
| H | `audit-mesh-run` | F8, X2 for the rest | about -60 (about 550 moved) | G, E | `@perf`: bench tool seam and byte-identical meshes |
| tools | `audit-tools-git` | section 4 | about -20 (tools are not production; governed files need Ola) | nothing | Ola's approval per governed file |

Total: about -345 production lines, the sum of the rows above, tools
included (B at its measured +22, F at its estimated +15; F's measured +29,
section 11, makes it about -330; the first estimate was about -440), about
-100 test lines (T2 came out at +1 and T1 is re-estimated at about -100,
against the -385 first estimated), and the drift points
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
- The scanner imports each module (`inspect.getsource`), and every L3
  module imports `_core`; so the suite needs the built extension, as the eight tests it replaces already do.

### The table

One dict, module (short name, `tin_engine.` dropped; the package root is
`tin_engine`, the extension `_core`) to (layer, the exact set of first-party
modules it imports). Rows as today's code has them, read with the fixed
scanner: equality, not a subset, so the table is the map and a new edge is a
visible table edit. Layers are section 5's:

- L0: `io.models`, `features`, `sources`, `run_record`, `stats`, `palettes`,
  `tin_engine`
- L1: `crs`, `mosaic`, `target_grid`, `grid_domain`, `domain`, `chains`,
  `elevation`, `outline`, `burn`, `gauge`, `sensitivity`, `reference`,
  `landcover`, `decompose`, `viz`, `viz.fixtures`, `viz.protocols`,
  `viz.scene`, `viz.style`, `viz.svg`
- L2: `io`, `io.cog`, `io.geojson`, `io.geopackage`, `io.geotiff`, `io.gml`,
  `io.mesh_index`, `io.ply`, `io.repository`, `io.rivers`, `io.station_set`,
  `io.vtk_legacy`, `fetch.http`
- L3: `raster`, `edge_strip`, `final_check`, `_core`
- L4: `dem_input`, `feature_input`, `catchment`, `catchment_batch`, `fetch`,
  `fetch.plan`, `fetch.run`, `fetch.nve`
- L5: `cli`

That is all 52 modules plus `_core`. A second dict, `UPWARD`, holds F12's
edges, each with the PR that removes it as its value: eight found, seven since
#191 resolved `fetch.http` -> `tin_engine`.

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

## 9. PR B design: `audit-crs-helpers` (F3)

Branch `worktree-audit-crs`, from T2's head `b63132e`; it lands after T2
and does not depend on T1. `src_python/` at `b63132e` is byte-identical to
`44fa7f5` and to `12dace7` (`git diff --stat 12dace7 b63132e -- src_python`
is empty), so the citations below are pinned to `44fa7f5`, which is on
master. Not refine or mesh code, and the benchmark runs pass no `--out-crs`
(`grep -n out-crs tools/bench.py` is empty), so no `@perf` run.

### What F3 got wrong, re-measured

- **The `EPSG:None` refusal is latent, not live.** Every `TileFootprint`
  comes from the GeoTIFF reader (`src_python/tin_engine/io/repository.py@44fa7f5:112-114, 241`),
  and the reader resolves a CRS only from an EPSG code: a projected CRS given
  by parameters (`ProjectedCSTypeGeoKey` 3072 = 32767, "user-defined") is
  refused at `src_python/tin_engine/io/geotiff.py@44fa7f5:283-288`. So no real
  footprint has `epsg` None, and the `EPSG:None` path is reached only by a
  footprint built by hand, as F3's probe did. It is still fixed here: the
  three "one CRS" checks become one, and the day the reader admits such a CRS
  the domain path must not break.
- **Ola's Austrian openDEM file is refused by the reader, not by `EPSG:None`.**
  Probe, with the project venv: `io.geotiff.read_meta` on
  `../rasputin_data/austria_dgm10/dhm_at_lamb_10m_2018.tif` raises
  `GeoTiffError: ProjectedCSTypeGeoKey (3072) = 32767 is not a resolvable
  EPSG code`. Its GeoKeys give Lambert conic conformal (two standard
  parallels) by parameters: parallels 46 and 49, origin 47.5 N
  13.33333333300013 E, false easting and northing 400 000 m, on MGI (EPSG
  4312), metres. That is EPSG:31287 to within 3.3e-10 degrees of longitude
  (0.03 mm at the origin's latitude). `kaprun_dgm10_31287.tif`, the converted copy, reads
  as `EPSG:31287`. **PR B does not make the original file readable**;
  question 1 below asks whether a later PR should.
- **Nine comparison sites, not eight.** F3 missed
  `src_python/tin_engine/dem_input.py@44fa7f5:164` (`target != parse_crs(first.crs)`, with `target`
  parsed three lines above), which decides whether the DEM is resampled.
  The red step found three more, all in the `==` form, which this design
  had missed too (`grep -rnE "(==|!=)" src_python/tin_engine` over CRS
  values): `src_python/tin_engine/domain.py@44fa7f5:62`,
  `src_python/tin_engine/feature_input.py@44fa7f5:300` and
  `src_python/tin_engine/fetch/plan.py@44fa7f5:118`. `@tester`'s
  `5074ef6` step found a thirteenth, missed by that grep because it calls
  `_parsed`, not `parse_crs`: `src_python/tin_engine/domain.py@44fa7f5:103`.
  Thirteen in all; the final sweep below is typed, not a grep, and finds no
  more. `src_python/tin_engine/mosaic.py@44fa7f5:585` (`ma.crs != mb.crs`) is
  not one: it compares tile texts, as `single_crs` does.
- **The live bug is the comparison itself.** `parse_crs(a) != parse_crs(b)`
  is pyproj's `CRS.__eq__`, PROJ's equivalence with axis order and parameter
  layout counted. It calls a CRS different from the EPSG code it is by
  definition, when it is spelt as a PROJ string or as GDAL's WKT1. Probe
  (pyproj 3.8.0, PROJ 9.8.1): `+proj=lcc +lat_1=46 +lat_2=49 +lat_0=47.5
  +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=bessel +units=m
  +no_defs` `==` `EPSG:31287` is False, and so is EPSG:3035's own
  `to_wkt("WKT1_GDAL")` against `EPSG:3035` (WKT1 has no axis order, so
  PROJ reads it east-north; EPSG:3035 is north-east). Today that means:
  `--out-crs` given as a PROJ string of the DEM's own CRS resamples the DEM
  for nothing (`src_python/tin_engine/dem_input.py@44fa7f5:164`); a river
  file, feature file or seed in such a spelling is refused as "not the DEM's"
  (`src_python/tin_engine/catchment.py@44fa7f5:185, 202`,
  `src_python/tin_engine/cli.py@44fa7f5:2045`,
  `src_python/tin_engine/feature_input.py@44fa7f5:423`); a GeoJSON domain
  in EPSG:25833 given `--domain-crs` as the PROJ string of EPSG:25833 is
  refused as disagreeing with itself
  (`src_python/tin_engine/domain.py@44fa7f5:103`); and the record says a transform
  ran where none did (`src_python/tin_engine/cli.py@44fa7f5:932-933, 952-956`).

### Prior art: legacy and literature

*Literature.* No published method; the rule leans on PROJ's own notions.
Since code review round 1 it uses two: equivalence, and the operation PROJ
builds between the two CRSs being its `noop` ("Pass a coordinate through
unchanged", `proj.org/en/stable/operations/conversions/noop.html`).
Identification, below, was the rule's second leg until then and is kept
here as the record of why it was dropped.
Equivalence is PROJ's `isEquivalentTo`, which pyproj's `CRS.equals` calls.
Identification is PROJ's `identify`, which pyproj's `list_authority` and
`to_epsg` call. Its confidence levels for a projected CRS, from PROJ's C++
reference (`proj.org/en/stable/development/reference/cpp/crs.html`,
`ProjectedCRS::identify`): 100, name and definition match; 90, equivalent,
names not exactly the same; 70, "CRS are equivalent (equivalent base CRS,
conversion and coordinate system), but the names are not equivalent"; 50,
"equivalent base ellipsoid and conversion, but the coordinate system do not
match (e.g. different axis ordering or axis unit)"; 25, "not equivalent, but
there is some similarity in the names". Nothing new is claimed.

*Legacy.* `git grep -nE "to_epsg|CRS\.equals|\.equals\(|is_exact_same|same_crs|IsSame|crs ==|crs !=|epsg ==|epsg !=" legacy-archive -- legacy`
returns one line, `legacy/rasputin/globcov_repository.py:141`
(`if pts_crs != self.data_crs:`), pyproj's `!=` as today. Nothing to carry
across.

### The helpers, in `crs.py` (layer L1)

They go in the existing `src_python/tin_engine/crs.py`, the layer-1 module
that already holds `parse_crs` and the one `Transformer.from_crs`. No new module, no new
first-party import: `crs`'s row in `tests/python/test_layering.py` stays
empty, and every caller already imports `crs`, so **the layering table does
not change** and there is no `UPWARD` entry to add or remove. Check 2 of
that test confirms it.

```python
def same_crs(a: str | CRS, b: str | CRS) -> bool:
    """Whether coordinates in `a` and in `b` name the same points: PROJ calls
    them equivalent once both are in x-then-y order, and the operation it
    builds between them is its `noop`."""

def transform_label(src: str | CRS, dst: str | CRS) -> str:
    """'none' when same_crs(src, dst), else transform_description(src, dst)."""

def single_crs(texts: Iterable[str], refusal: type[ValueError] = ValueError) -> str:
    """The one CRS text among `texts` (a DEM's tiles' `meta.crs`), or
    `refusal` naming them all."""
```

**`same_crs`, the rule** (ruled after code review round 1; it replaces the
rule `5197f9a` implements, whose second leg was EPSG identification plus a
frame check). "The same" means: coordinates in `a`, read as coordinates in
`b`, name the same points, so the transform a site would run can be skipped.
Build the package's one transformer, `_transformer(a, b)` (`always_xy=True`;
unreadable text is `parse_crs`'s `ValueError`, wording unchanged). The two
are the same when both hold:

1. **Equivalent once both are in x-then-y order:**
   `source_crs.equals(target_crs)` on the transformer's own copies, which
   PROJ gives in x-then-y order (probe: EPSG:3035's two axes come back east,
   north). Axis order never matters in this package: every transform is
   `always_xy` (the module docstring of `crs.py`). This leg settles the
   datum. The operation cannot: PROJ's null and ballpark transformations
   move no point either (EPSG:4258 to 4326 is a `noop` of 1 m accuracy; a
   PROJ string with no datum reaches an EPSG datum through a "ballpark"
   `noop`).
2. **The operation is PROJ's `noop`:** the first word of
   `transformer.definition` is `proj=noop`. This leg settles the conversion,
   from what the transform does rather than from attributes: PROJ builds the
   pipeline from everything both definitions hold, including a parameter
   kept only in a PROJ string's remark, so it sees what the CRS object does
   not show. It reads the operation's text and moves no point, so red tests
   8-10's guard still holds. Only `noop` passes: under `always_xy` both
   sides are already x-then-y, so a remaining `axisswap` flips or swaps an
   axis and a `unitconvert` rescales, and both move points (UTM 33 with
   `+axis=wnu` against EPSG:25833 is one `axisswap order=-1,2`, 1 550 km).

A `ProjError` while building the transformer (PROJ finds no operation, for
example between Earth and Mars) is "not the same", never an exception. So
is a `source_crs` or `target_crs` of None: pyproj types both `CRS | None`
(None for a transformer built from a pipeline), so mypy needs the branch,
though no `from_crs` transformer has one (probe: 0 of 141 pairs over 13
CRSs, geocentric, compound, bound, rotated and Mars among them). "Not the
same" is the safe side: a transform runs and the record names it, or a
refusing site refuses, as on master.

Neither leg alone is enough (probe table below). Equivalence alone calls
`+proj=longlat +datum=WGS84 +lon_0=10` the same as `+lon_0=20`, 10 degrees
apart: PROJ keeps `+lon_0` of a `longlat` only in the remark, and pyproj's
`==` on master says the same, so that leak is older than PR B. The `noop`
alone calls EPSG:4258 the same as 4326, a `+towgs84=0` UTM string the same
as EPSG:25833, and a 3D CRS the same as its 2D one.

**What it drops from `5197f9a`:** `SAME_CONFIDENCE`, `FRAME_TOLERANCE`,
`_same_frame` and the EPSG code sets. Three rounds in a row (red tests
8-10, the derived probe set, code review round 1) each found a CRS
attribute that identification ignores (axis direction, prime meridian, and
now a conversion parameter held only in the remark); judging the
transform ends that line of leaks, rather than adding a fourth attribute.

**What "the same" tolerates: PROJ's own equivalence, no constant of
rasputin's.** Perturbing one parameter of EPSG:25833 or EPSG:31287, given as
WKT2 without its ID (datum named): the scale factor off by 2e-10 relative
is still the same, 1e-9 is not; the false easting, the longitude of origin
and the first standard parallel flip between 1.1e-10 and 2e-10.
Within that tolerance PROJ's own transforms from the geographic CRS into
either side give the same coordinates, 0.000 mm apart over a 15 by 15 grid
on the code's area of use, so skipping the transform changes nothing PROJ
would have produced. At the first step outside it the moves are 0 to
9.3 mm, and the rule says not the same. Checked on pyproj 3.8.0 / PROJ
9.8.1; the largest input is UTM 33N's area of use (12-18 E, 34.8-84 N). The
earlier "under 2 mm" bound (identification admitting a scale factor off by
2e-10, 1.9 mm at 84 N) is withdrawn: it measured identification, which the
rule no longer uses. The Austrian file's longitude of origin
(13.33333333300013 against EPSG's 13.3333333333333) is inside PROJ's
tolerance: EPSG:31287 written as WKT2 with that value is the same as
EPSG:31287 (by name; renamed "unknown" it is not, see question 1).

*Rejected, each on the probe set below ("breaks" counts triples where a is
the same as b and b as c, but a is not the same as c):*

- **Identification at confidence 70, with the frame check** (`5197f9a`): calls
  the `+lon_0` strings the same as EPSG:4326 and 4269; 164 breaks over
  the set's 85 CRSs. A PROJ
  string without a datum identifies as up to 11 codes on different datums
  (`+proj=utm +zone=33 +ellps=GRS80` as EPSG:25833, 3006, 3045 and eight
  more), so the old rule guessed the datum and the guess was not transitive.
- **The `noop` alone, or equivalence or the `noop`** (the "pipeline alone"
  reading): calls EPSG:4258 the same as 4326 and as 4269, and a
  `+towgs84=0` UTM string, EPSG:4937 (3D) and EPSG:25833+5941 (with
  height) the same as their plain codes.
- **The `noop` and either equivalence or a code in common**: keeps the
  datum-less strings the same, and is not transitive (pyproj's PROJ string of
  EPSG:3035 is the same as EPSG:3035, which is the same as its GDAL WKT1,
  but the string is not the same as the WKT1).
- **Check the points instead**: transforming a grid over the code's area of
  use and calling the pair the same when none moves more than 1 mm. It
  calls `Transformer.transform` inside `same_crs`, the very method red
  tests 8-10 refuse to prove that no point moved, and needs a tolerance and
  an area of use that a CRS without a code lacks.
- **"Only `unitconvert` and `axisswap` steps"**: wrong under `always_xy`
  (rule, leg 2 above).

Cost: one transformer, as before, without the identification calls: 2 to
12 ms a call (probe: EPSG:25833 against 3045 9 ms; the `+towgs84` Lambert
12 ms, 16 ms with identification). Calls happen once per site per run,
never per point.

**Transitive on the probe set:** 96 ordered same pairs over 85 CRSs, 0
asymmetric, 0 breaks. Not proven in general, and not needed for soundness:
each "the same" is PROJ's `noop` between the two, so a chain of them (a
waters file the same as the river file, the river file the same as the
DEM's CRS, `cli.py` and `catchment.py`) composes `noop`s and moves no point
either. **Not reflexive for one exotic spelling:** UTM 33 with `+axis=esu`
against itself is a pipeline of two `axisswap order=1,-2`, which PROJ does
not fold to `noop`, so it is "not the same" as itself (the safe side; no
reader writes such a CRS).

**The probe set.** Derived from what a CRS is made of (ISO 19111: a datum
with its ellipsoid and prime meridian, a coordinate system with its axes'
order, direction and unit, a conversion with its method and parameters, a
dimension, an epoch), not from the bug as found, plus the reviewer's pairs.
Each row perturbs one component of EPSG:25833 (UTM 33N), EPSG:31287 (the
Lambert), EPSG:4258 or EPSG:4326, as a PROJ string unless named. "At
`5197f9a`" is the rule that commit implements; "Now" is the rule above;
"Moves" is the largest distance, over a 15 by 15 grid on the EPSG code's
area of use, between a point and its image under the always_xy transform
from one to the other. Every row gives the same answer in both argument
orders, under both rules. The probe script is
`docs/increments/python-audit-probes/same_crs_probe.py` (run from the
repository root with the project venv; about 5 minutes).

| Component | Probe | At `5197f9a` | Now | Moves |
|---|---|---|---|---|
| axis order | EPSG:4326's GDAL WKT1; EPSG:25833, 31287 and 3035 as WKT2 without ID, axes swapped | same | same | 0 |
| axis order, no datum | UTM 33 `+axis=neu` | same | not | 0 |
| axis direction | UTM 33 `+axis=wnu`, `esu`, `wsu`; the Lambert `+axis=wnu`; geographic GRS80 `+axis=wnu` against EPSG:4258 | not | not | 1 390 to 18 700 km |
| linear unit | UTM 33 in US feet, in km, `+to_meter=1.0000001`, `1.000000001`; the Lambert in US feet | not | not | 9 mm and more |
| linear unit | EPSG:2263 (New York Long Island, US feet) against 32118 (the same in metres) | not | not | 1 140 km |
| angular unit | EPSG:4326 as PROJJSON without ID, axes in grads | not | not | 2 490 km |
| prime meridian | UTM 33 `+pm=paris`; geographic GRS80 `+pm=ferro` against EPSG:4258; EPSG:31251 (MGI Ferro) against 31254 (MGI) | not | not | 0 to 1 960 km |
| conversion parameter held in the remark | `+proj=longlat +datum=WGS84 +lon_0=10` and `+lon_0=-3` against EPSG:4326; `+proj=longlat +datum=NAD83 +lon_0=10` against EPSG:4269; `+lon_0=10` against `+lon_0=20` | **same** | not | 334 to 1 110 km |
| conversion parameter held in the remark | `+proj=longlat +datum=WGS84 +lon_0=10` against OGC:CRS84 | not | not | 1 110 km |
| datum | UTM 33 `+datum=WGS84`, `+ellps=WGS84`, `+towgs84=0,0,0,0,0,0,0`, `+towgs84=0,0,0`; the Lambert with `+nadgrids=@null` and with its `+towgs84` | not | not | 0 to 119 m |
| datum | geographic GRS80 against EPSG:4258, 4269, 4283; EPSG:4269 against 4258; 4258 against 4326; 26917 against 6346 | not | not | 0 |
| no datum named | the Lambert by parameters (either parallel order) and pyproj's PROJ strings of EPSG:31287, 25833 and 3035, each against its code | same | **not** | 0 |
| method, no datum | UTM 33 as `+proj=tmerc`, `+proj=etmerc`, `+proj=tmerc +approx` | same | not | 0 to 0.008 mm |
| parameters, no datum | the Lambert's longitude of origin off by 2e-9 degrees; the scale factor off by 2e-10; the false easting off by 1e-5 m; `+k=0.5` and `+x_0=1000` on `+proj=utm` (PROJ ignores both) | same | not | under 2 mm |
| parameters | the Lambert at 13.5 E; EPSG:25832 against 25833 | not | not | 12.8 km and more |
| hemisphere | UTM 33 `+south` | not | not | 10 000 km |
| dimension | EPSG:25833+5941 (with height) against 25833; EPSG:4937 (3D) and 4936 (geocentric) against 4258 | not | not | 0 or more |
| vertical | UTM 33 `+vunits=ft`; `+geoidgrids=@egm96_15.gtx` | not | not | 0 |
| longitude range | geographic GRS80 `+lon_wrap=180`, `+over`, against EPSG:4258 | not | not | 0 or 360 degrees |
| derived | rotated pole (`+proj=ob_tran`) against EPSG:4258 | not | not | 6 730 km |
| epoch | WGS 84 (G1762) with frame epoch 2020 against EPSG:9057 (frame epoch 2005) | not | not | 0 |
| same by definition | WKT2 without ID of EPSG:31287 and 25833; GDAL and ESRI WKT1 of 3035 and 25833; OGC:CRS84 and 4326; 3045 and 25833; `epsg:31287`; ESRI:102100 and EPSG:3857; `+proj=utm +zone=33 +datum=WGS84` and 32633; `+proj=longlat +datum=WGS84` and 4326; `+proj=longlat +datum=NAD83` and 4269; the Long Island Lambert with `+datum=NAD83 +units=us-ft` and 2263 | same | same | 0 |
| other body | Mars `+proj=longlat +a=3396190 +b=3376200` against EPSG:4326 | not | not | no operation |

The earlier table's claim that its two bold rows (axis direction, prime
meridian) were "the only ones where the code check alone calls the same
what is not" is **withdrawn**: the remark-held `+lon_0` row is a third, and
the "no datum named" row shows identification guessing a datum. The
earlier row "EPSG:25833 with coordinate epoch 2020.0 ... same (rule 1)" is
withdrawn too: pyproj 3.8.0 builds no CRS with a coordinate epoch
(`EPSG:25833@2020.0` is refused, "Coordinate epoch should not be provided
for a static CRS"), so no CRS text reaches `same_crs` with one; the epoch
row above perturbs a frame epoch instead.

**Pinned limits of the rule** (each a test below):

- **A PROJ string that names no datum is not the same as any EPSG code**,
  pyproj's own `to_proj4()` of EPSG:25833 and the Lambert by parameters
  included. Such a string says "on the GRS80 (or Bessel) ellipsoid", not
  "on ETRS89" (or MGI); PROJ reaches the code's datum through a ballpark.
  Where the site reprojects, that ballpark runs, moves nothing, and the
  record names it ("Ballpark geographic offset from unknown to ETRS89");
  where the site refuses, it refuses as on master. A PROJ string with
  `+datum=WGS84` or `+datum=NAD83` names its datum and is the same as its
  code.
- **A PROJ string with `+towgs84`** (EPSG:31287 written with
  `+towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232`, or UTM 33
  with `+towgs84=0,0,0` or seven zeros) is **not** the same as its code:
  pyproj reads it as a bound CRS on an unknown datum.
- **OGC:CRS84 is the same as EPSG:4326**: coordinates are always x then y
  here, so they name the same points. `crs_label` still writes each by its
  own name.
- **EPSG:3045 is the same as EPSG:25833**: one definition under two codes.
- **An axis pointing the other way, another prime meridian, or a
  `+lon_0` on a `longlat`, is not the same.**

**Labels do not change.** `crs_label` and `target_grid` keep
`to_epsg(min_confidence=100)` (`src_python/tin_engine/crs.py@44fa7f5:77`,
`src_python/tin_engine/target_grid.py@44fa7f5:200`): a record names the CRS as it
was given, not as identified. When `--out-crs` is the same as the DEM's CRS
the DEM is not resampled, and the record's `crs` is the DEM's own,
`opened.tile.meta.crs` (`src_python/tin_engine/cli.py@44fa7f5:883`), as it is today for
`--out-crs EPSG:25833` on an EPSG:25833 DEM.

**`single_crs`** is keyed on the `meta.crs` text, as two of today's three copies
are. Every tile the reader makes has `EPSG:n` text, so keying by
`same_crs` would buy nothing and cost a transformer per pair. It replaces:
`src_python/tin_engine/dem_input.py@44fa7f5:188-190` (`MosaicError`),
`src_python/tin_engine/dem_input.py@44fa7f5:227-229` (`MosaicError`; the
`epsg` copy, and with it the `EPSG:None` bug), and
`src_python/tin_engine/catchment.py@44fa7f5:194-197` (`CatchmentError`).
The domain path then moves the domain into that text
(`given.to_crs(crs)`, now `f"EPSG:{epsgs[0]}"`), and its re-raise
(`src_python/tin_engine/dem_input.py@44fa7f5:248`) names it
(`in the DEM's {crs}`). That reads `EPSG:25833` for every tile the reader
makes, so the wording is unchanged there.

### The sites, as they change

| Site (at `44fa7f5`) | Today | After |
|---|---|---|
| `src_python/tin_engine/dem_input.py@44fa7f5:164` | `target != parse_crs(first.crs)` | `not same_crs(target, first.crs)` |
| `src_python/tin_engine/dem_input.py@44fa7f5:188-190` | own "one CRS" check | `crs = single_crs((f.meta.crs for f in footprints), MosaicError)` |
| `src_python/tin_engine/dem_input.py@44fa7f5:227-231, 248` | keyed on `epsg` | `single_crs(..., MosaicError)`; the domain moved into and the re-raise naming that text |
| `src_python/tin_engine/catchment.py@44fa7f5:185` | `!=` | `not same_crs(crs, dem_crs)`; wording unchanged |
| `src_python/tin_engine/catchment.py@44fa7f5:194-197` | own "one CRS" check | `dem_crs = single_crs(..., CatchmentError)` |
| `src_python/tin_engine/catchment.py@44fa7f5:202` | `!=` | `not same_crs(...)`; wording unchanged |
| `src_python/tin_engine/feature_input.py@44fa7f5:179` | `!=` | `not same_crs(...)` |
| `src_python/tin_engine/feature_input.py@44fa7f5:423` | `!=` | `not same_crs(...)`; wording unchanged |
| `src_python/tin_engine/cli.py@44fa7f5:932-933` | `same = ...; how = ...` | `"domain_transform": transform_label(given.crs, dem_crs)` |
| `src_python/tin_engine/cli.py@44fa7f5:952-956` | five-line conditional | `transforms.append(transform_label(own, dem_crs))` |
| `src_python/tin_engine/cli.py@44fa7f5:2045` | `!=` | `not same_crs(...)`; wording unchanged |
| `src_python/tin_engine/fetch/run.py@44fa7f5:222` | `!=` | `not same_crs(...)`; wording unchanged |
| `src_python/tin_engine/domain.py@44fa7f5:62-63` | `if source == target: return self` | `if same_crs(source, target):` return the same polygon labelled `target.to_string()` (`self.model_copy(update=...)`) |
| `src_python/tin_engine/feature_input.py@44fa7f5:300` | `src == self.dem` | `same_crs(src, self.dem)` |
| `src_python/tin_engine/fetch/plan.py@44fa7f5:118` | `frame == source` | `same_crs(frame, source)` |
| `src_python/tin_engine/domain.py@44fa7f5:103` | `_parsed(crs) != _parsed(own)` | `not same_crs(_parsed(crs), _parsed(own))`; wording unchanged, and unreadable text is still `_parsed`'s `DomainError` |

`src_python/tin_engine/cli.py@44fa7f5:919` (`transform_description` on the resampled path) stays: a resampled DEM's
CRS is never the target's. `parse_crs` stays public; `cli` and `catchment`
stop importing it.

### The final sweep (after `5074ef6`)

Three rounds each missed a site, because each searched by spelling. This one
searches by type and by data flow. Rerun it from the repository root, at
any commit before the green one (`src_python/` and `tools/` on this branch
are `44fa7f5`'s until then; `git diff --stat 44fa7f5 HEAD -- src_python
tools` is empty):

```bash
.venv/bin/python docs/increments/python-audit-probes/crs_sweep.py
```

What it lists is in the script's docstring: every comparison operator, the
set, dict, subscript and method forms of keying, and `match`, wherever one
side may hold a CRS. "May hold a CRS" is: mypy infers a pyproj CRS type,
or a first-party model with a `crs`, `epsg` or `srs` field (whole-model
`==`); or the text names a CRS; or the value is derived from a pyproj CRS
(`.to_string()`, `.to_epsg()`), assigned from, passed into a parameter
from, or returned as such a value, to a fixed point across files. A text
with no CRS-named source is listed in the comparison, membership and set
forms anyway, as "no CRS name". It over-reports by design: 524 lines over 65
files, 133 of them "no CRS name". Two checks that it can fail:

- **It walks every comparison.** stderr prints `mypy walk: 1182
  comparisons` and `ast: 1182 comparisons`; the two must agree.
- **It finds planted shapes.** In a scratch copy, eight planted comparisons
  were each listed: two parsed CRSs under names `a`, `b`;
  `p.to_string() == q.to_string()`; two `RasterMeta`s by `!=`; two `str`
  parameters `u == v`; `u in seen`; `len({u for u in texts}) > 1`;
  `_k(a) == _k(b)` with `_k` returning `(m.crs, m.delta_x)`; and two CRS
  texts passed into parameters `x`, `y` used as dict keys. The earlier
  version of the script, without the function-return rule, missed the
  `_k` shape, and so missed the real `src_python/tin_engine/mosaic.py@44fa7f5:307`
  (`_key(a) == _key(b)`, whose tuple starts with `m.crs`).

Not listed: a CRS text with no CRS-named source, used only as a subscript or
a dict display's key, or typed `Any`. Every CRS source in the package has a
CRS name (a model field, a `--*-crs` option, a `"crs"` or `"srsName"`
member, `parse_crs`), so such a text would have to come from data under
another key.

**The full list of CRS comparisons**, sorted by hand from the 524 lines (each pinned to `44fa7f5`):

*Two CRSs by pyproj's `==` / `!=`: the thirteen sites, all to `same_crs`*
(the table above): `src_python/tin_engine/dem_input.py@44fa7f5:164`, `src_python/tin_engine/catchment.py@44fa7f5:185, 202`,
`src_python/tin_engine/feature_input.py@44fa7f5:179, 300, 423`, `src_python/tin_engine/cli.py@44fa7f5:932, 954, 2045`,
`src_python/tin_engine/fetch/run.py@44fa7f5:222`, `src_python/tin_engine/fetch/plan.py@44fa7f5:118`, `src_python/tin_engine/domain.py@44fa7f5:62, 103`. With the
script's `NAMES` and `FIELDS` set to match nothing (types and data flow
alone), it lists 21 comparisons outside "no CRS name": these thirteen, seven
that test one CRS or its parts against None or a constant, and one float
comparison (`src_python/tin_engine/io/geotiff.py@44fa7f5:384`) reached through a shared parameter name.

*Two CRS texts compared as text: they stay text, ruled here.*

| Site | What it compares | Ruling |
|---|---|---|
| `src_python/tin_engine/dem_input.py@44fa7f5:188-189`, `227-228`; `src_python/tin_engine/catchment.py@44fa7f5:194-195` | the tiles' texts, as a set | to `single_crs`, keyed on text (ruled above) |
| `src_python/tin_engine/mosaic.py@44fa7f5:307` (`_aligned`) | `_key(a) == _key(b)`, the tile text first | stays: tiles on one lattice; `single_crs` has already made them one text, so it never splits a CRS |
| `src_python/tin_engine/mosaic.py@44fa7f5:585` (`_mixed`) | `ma.crs != mb.crs` | stays (ruled above): words why `_aligned` said no |
| `src_python/tin_engine/mosaic.py@44fa7f5:267` | `plan.tiles[0].meta == plan.meta` (the text inside) | stays: a single-tile shortcut; `plan.meta` takes its text from the tiles (`src_python/tin_engine/mosaic.py@44fa7f5:416`), and a miss only takes the general path |
| `src_python/tin_engine/mosaic.py@44fa7f5:484` | `tile.meta != listed` | stays: the same file read twice by the same reader; any change, a respelt CRS included, means the file changed |
| `src_python/tin_engine/io/gml.py@44fa7f5:70-73` | the `srsName`s of one GML file, as a set | stays: one file, one writer; a file spelling one CRS two ways is refused naming both, never read wrongly |
| `src_python/tin_engine/io/geopackage.py@44fa7f5:171` | a geometry's `srs_id` against its layer's | stays: row ids in the file's own CRS table, and the GeoPackage standard requires them equal |

*One CRS against a constant or its own parts (not two CRSs):*
`src_python/tin_engine/io/geotiff.py@44fa7f5:137, 287, 293, 303, 321, 325, 331, 337`,
`src_python/tin_engine/io/geopackage.py@44fa7f5:121`, `src_python/tin_engine/io/models.py@44fa7f5:71-73`, `src_python/tin_engine/crs.py@44fa7f5:78`,
`src_python/tin_engine/dem_input.py@44fa7f5:154-156`, `src_python/tin_engine/feature_input.py@44fa7f5:302, 304, 487`,
`src_python/tin_engine/fetch/run.py@44fa7f5:245`.

Every other line is not a CRS: a name that matches (`source` is a catalogue
source, `target` an output path, `code` a land-cover code) or an unrelated
`str`. Nothing in `tools/` compares a CRS.

### Refusal wordings

One changes, by design: the three "one CRS" refusals become one, in plain words.

| Where | Today | After |
|---|---|---|
| `single_crs` (all three sites) | `the tiles are in 2 CRSs, ['EPSG:25832', 'EPSG:25833']; need one`, and on the domain path `... EPSG:[25832, 25833]; a domain needs one` | `the DEM files are in 2 different CRSs (EPSG:25832, EPSG:25833); all must be in one CRS` |

Two more gain a hint, by Ola's ruling D15 b ("After the round-1 red step"
below): the two refusals that compare a river file's or a request's CRS
against the DEM's.

| Where | Today | After, when the hint applies |
|---|---|---|
| `check_reach_crs` (`src_python/tin_engine/catchment.py@44fa7f5:185-186`) | `the river file's CRS, X, is not the DEM's, EPSG:n` | the same, then `; if you mean EPSG:n, write EPSG:n` |
| `delineate`'s reach check (`src_python/tin_engine/catchment.py@44fa7f5:202-203`) | `the river reach must be in the DEM's CRS, EPSG:n` | the same, then `; if you mean EPSG:n, write EPSG:n` |

Every other wording is unchanged and stays pinned by its suite:
`cannot read the CRS ...` (`parse_crs`), `the river file's CRS, X, is not the
DEM's, Y`, `the river reach must be in the DEM's CRS, Y`, `F is in X but the
given CRS is Y`, `the W file's CRS, X, is not the river file's, Y; polygons
are not reprojected`, `O: the header's CRS is X; the catalogue's S is Y`.
Question 2 asks Ola to approve the new one.

### Red tests (`@tester`, one commit, before any code)

`crs.same_crs`, `transform_label` and `single_crs` do not exist yet, so each
new `test_crs.py` test fails on its own with `AttributeError`, inside the
existing `crs` fixture. Lean: no throwaway implementation, no mutation round.

1. **`tests/python/test_crs.py`, `TestSameCrs`**, each pair in both orders:
   - the same: the Lambert PROJ string above and `EPSG:31287`; the WKT of
     EPSG:31287 with its ID removed; `EPSG:3035`'s `to_wkt("WKT1_GDAL")` and
     `EPSG:3035`; `CRS.from_epsg(25833).to_proj4()` and `EPSG:25833`;
     `OGC:CRS84` and `EPSG:4326`; `EPSG:3045` and `EPSG:25833`; and, as the
     control, `EPSG:31287` and `epsg:31287`.
   - not the same: the Lambert with `+units=us-ft`; the Lambert at
     `+lon_0=13.5`; `EPSG:25832` and `EPSG:25833`; `EPSG:4258` and
     `EPSG:4326`; that `+towgs84` PROJ string of EPSG:31287 (the pinned
     limit); a Mars `+proj=longlat +a=3396190 +b=3376200` and `EPSG:4326`
     (False, not an exception); UTM 33 with `+axis=wnu` and `EPSG:25833`
     (added in `7dddda8`); UTM 33 with `+pm=paris`
     (`+proj=utm +zone=33 +ellps=GRS80 +units=m +pm=paris +no_defs`) and
     `EPSG:25833` (added after `7dddda8`, the frame check's second half).
   - unreadable text is a `ValueError` matching `cannot read the CRS`.
2. **`TestTransformLabel`**: `"none"` for the Lambert string against
   `EPSG:31287`; `transform_description("EPSG:4326", "EPSG:25833")` for
   that pair.
3. **`TestSingleCrs`**: `["EPSG:25833"] * 3` gives `"EPSG:25833"`;
   `["EPSG:25833", "EPSG:25832", "EPSG:25833"]` raises the given type (a
   local `ValueError` subclass), with exactly the wording above; with no
   type given, `ValueError`.
4. **Amend `tests/python/test_dem_input_domain.py`, `TestOneCrs`**
   (`tests/python/test_dem_input_domain.py@44fa7f5:630-632`): the three
   `in message` assertions become the new wording,
   `"(EPSG:25832, EPSG:25833)"` and `"all must be in one CRS"` among them.
5. **The `EPSG:None` fix, `test_dem_input_domain.py`**: `dem_input.open_dem`
   with `dem_input.repository_for` monkeypatched to return an in-memory
   repository (`catchment_fixtures.MemoryRepository` with a no-op `check` and
   `load_window = None`). Its tiles have `epsg=None` and `crs` the WKT of
   EPSG:25833 with its ID removed. With a domain, the tile's array and its
   lattice (`x_min`, `y_max`, spacings, rows, columns) equal those of the
   same tiles with `epsg=25833`, and its `crs` is that WKT text. A domain
   outside the tiles is refused naming that CRS text, and no message
   contains `EPSG:None`.
6. **The parameter Lambert is accepted where `EPSG:31287` is,
   `tests/python/test_catchment_batch.py`** beside
   `test_check_reach_crs_is_the_one_rule`: a `MemoryRepository` whose first
   tile's meta has `epsg=31287`; `check_reach_crs(<Lambert string>, repo)` is
   None, as `check_reach_crs("EPSG:31287", repo)` is; the Lambert at 13.5 E
   is refused with today's wording (`river file's CRS ... DEM's, EPSG:31287`).
7. **`--out-crs` the same by definition, `tests/python/test_cli_mesh_geographic.py`,
   G7**: parametrise `test_out_crs_equal_to_the_dems_own_writes_the_same_bytes`
   over `"EPSG:25833"` (today's) and `CRS.from_epsg(25833).to_proj4()`. The
   second writes the same bytes as no `--out-crs`. Red today: the DEM is
   resampled.

Red today: 1 (all), 2, 3, 4, 5, 6 (the Lambert half) and 7 (the PROJ string).
The rest of the suite must stay green; no other test is expected to move
(the `domain_transform` and `features_transform` tests use CRS pairs that
differ).

**After the red step** (`29aff00`), one line each:

- `domain.py:62` uses `same_crs`: one rule everywhere; otherwise a domain spelt as a PROJ string of the DEM's CRS goes through a transform while the record says `domain_transform` "none".
- When the same, `to_crs` returns the same polygon labelled `target.to_string()`, not `self`: today's output exactly (probe: the PROJ-string domain comes back bit-identical, labelled `EPSG:25833`), and the result's `crs` is always `dst`'s.
- `feature_input.py:300` and `src_python/tin_engine/fetch/plan.py@618328b:118` use `same_crs` too, by the same rule; neither changes a wording.
- Red tests 8-10 below are needed: each site's output is the same today, so only a refused point-moving `Transformer` method can tell the fix from the bug.
- `@tester`'s departure, accepted: `TestTheSameCrs`'s guard refuses the point-moving methods (`transform`, `itransform`, `transform_bounds`), not `Transformer.from_crs`, since `same_crs` builds one to compare; the invariant (no point moved) is unchanged and the guard was shown still to catch a real transform.
- `@tester`'s departure, accepted: wording pins at the other two `single_crs` sites (the `--out-crs` path and `catchment.delineate`), beyond test 4's one.
- `29aff00` moved `tests/python/test_cli_mesh_geographic.py@44fa7f5:886`, cited by `docs/increments/h16-harness-fixes.md` line 579; that citation was pinned here to `44fa7f5` and by T1 to `97eea35`, and both read the same line; the merge with master kept T1's.

8. **`tests/python/test_domain.py`**: with `Transformer`'s `transform`,
   `itransform` and `transform_bounds` refused (as in `TestTheSameCrs`),
   `DomainPolygon(polygon=<a box in UTM 33>, crs=proj4_of(25833)).to_crs("EPSG:25833")`
   has `crs == "EPSG:25833"` and a polygon `equals_exact` to the given one
   at tolerance 0. Red today: the transform runs.
9. **`tests/python/test_feature_input.py`**: a feature source whose CRS is
   `proj4_of(<the DEM's EPSG>)`, read with the same three methods refused,
   gives the same geometries as the source spelt `EPSG:<n>`. Red today, and
   still red with only the table's earlier rows fixed (line 300 moves the
   points).
10. **`tests/python/test_fetch_plan.py`**: `source_box` with a `box` and
    `out_crs=proj4_of(<meta's EPSG>)`, with `transform_bounds` refused,
    equals `source_box` with `out_crs=None`. Red today.

**After red tests 8-10** (`7dddda8`), one line each:

- Ruled: rule 2 gains the frame check (axes and prime meridian, above); `@tester`'s candidate (axes only) is taken and widened, since the derived probe set found `+pm=paris` UTM 33 still called EPSG:25833 by it.
- Ruled: the "fraction of a millimetre" bound is withdrawn; identify admits a scale factor off by 2e-10, under 2 mm over UTM 33N's area of use (above).
- `@tester`'s departure, accepted: `test_domain.py`'s `no_transformer` guard refuses only the point-moving methods, as `TestTheSameCrs`'s does (one rule: building a transformer moves no point).
- `@tester`'s departure, accepted: `test_cli_mesh_domain_crs.py`'s `TestTheSameCrs.assert_as_16` compares the handed-on domain's CRS by pyproj equality, since `to_crs` now labels it with the target's text (ruled after the red step); the bits are still compared exactly.
- `@tester`'s departure, accepted: `proj4_of` and the guard move to a new shared helper, `tests/python/crs_fixtures.py`, used by five suites.
- The frame check needs one more not-the-same pair (`+pm=paris`, red test 1); `@tester` adds it before `@developer` starts, so it is red first.

11. **`tests/python/test_domain.py`, `TestReading`**, site 13: parametrise
    over `(member, flag)` in `("EPSG:25833", proj4_of(25833))`,
    `(UTM33, proj4_of(25833))` and `(proj4_of(25833), "EPSG:25833")`.
    `read_domain(path, flag)` on a GeoJSON square whose `crs` member is
    `member` is not refused; its `crs` is `member`, the file's own text (the
    flag never overrides it), and its polygon is `equals_exact`, tolerance 0,
    to that of `read_domain(path)` without the flag. Red today: `DomainError`
    `d.geojson is in EPSG:25833 but --domain-crs says +proj=utm ...`
    (probe, the project venv). Beside it, green today and to stay green: the
    flag `+proj=utm +zone=33 +ellps=GRS80 +units=m +pm=paris +no_defs`
    against the member `EPSG:25833` is refused with the wording
    `is in EPSG:25833 but --domain-crs says`, so the fix cannot widen past
    `same_crs`.
12. **`tests/python/test_crs.py`, `TestSameCrs`**, the None side: with
    `crs._transformer` monkeypatched to return an object whose `source_crs`
    is None and `target_crs` is `CRS.from_epsg(25833)`, and then the
    mirror, `same_crs("EPSG:25833", "EPSG:25833")` is False. Needed:
    `testing.md`'s edge-case rule (every condition the code handles
    specially has a named test), and no real input reaches the branch (the
    probe above). Red today: `same_crs` does not exist.

**After red test 1's `+pm=paris` pair** (`5074ef6`), one line each:

- Ruled: site 13, `domain.py:103`, uses `same_crs` with its wording unchanged (table above); red test 11.
- Ruled: the final sweep above replaces grep; its script is in `docs/increments/python-audit-probes/`, its list in this section, and it finds no fourteenth site; seven text comparisons stay text, each with its reason.
- Ruled: `float(pm.unit_conversion_factor)` is the accepted form for mypy (rule 2 above).
- Ruled: a None `source_crs` or `target_crs` is "not the same", with red test 12.
- Ruled: rule 2 reads its code sets from the CRSs as given, not from the x-then-y copy (rule 2 above; the other reading gets 6 of 44 same pairs wrong).

**After code review round 1** (`5197f9a`), one line each:

- Ruled: `same_crs` is PROJ equivalence once x-then-y **and** a `noop` operation (the rule above); identification, `SAME_CONFIDENCE`, `FRAME_TOLERANCE` and `_same_frame` go. The three previous "Ruled" lines about rule 2's code sets and `float()` lapse with it.
- Ruled: the reviewer's option ("only unit-convert and axis-swap steps") is narrowed to `noop` alone, since under `always_xy` any remaining step moves points; and it is a second leg, not the whole rule, since the `noop` alone calls EPSG:4258 the same as 4326.
- Ruled: a PROJ string that names no datum is no longer the same as its EPSG code (pinned limits above); question 3 is amended to say so.
- The probe set is a committed script, `docs/increments/python-audit-probes/same_crs_probe.py`; its 85 CRSs give 0 asymmetric pairs and 0 transitivity breaks under the rule.

`@tester`, one commit, before `@developer`. Rows marked *red* fail at `5197f9a`; the others are re-spellings that pass at `5197f9a` and on master would fail as the earlier red tests did.

- **`tests/python/crs_fixtures.py`**: add `axes_swapped(epsg)`, the WKT2 of `epsg` without its ID and with its two axes in the other order (pyproj's `==` on master calls it different; the rule calls it the same). Its docstring's first sentence ("a CRS spelt as a PROJ string of an EPSG code is that code") is no longer true; say "spelt by definition". Keep `proj4_of` for the not-the-same rows below.
- **Red test 1, `SAME`**: drop "Lambert by parameters" and "PROJ string of 25833"; add `CRS.from_epsg(31287).to_wkt("WKT1_GDAL")` and EPSG:31287, and `axes_swapped(25833)` and EPSG:25833.
- **Red test 1, `NOT_SAME`**, both orders: *red* the Lambert by parameters and EPSG:31287; *red* `proj4_of(25833)` and EPSG:25833; *red* `+proj=longlat +datum=WGS84 +lon_0=10 +no_defs` and EPSG:4326; *red* the same with `+lon_0=-3`; *red* `+proj=longlat +datum=NAD83 +lon_0=10 +no_defs` and EPSG:4269; *red* the `+lon_0=10` string and the `+lon_0=20` string.
- ***Red*, transitivity, `TestSameCrs`**: with S the `+lon_0=10` string, `same_crs("EPSG:4326", "OGC:CRS84")` is True, and `same_crs(S, "EPSG:4326")` and `same_crs(S, "OGC:CRS84")` are both False (at `5197f9a` S is the same as EPSG:4326 but not as CRS84).
- **`test_crs_objects_are_accepted`**: the True pair becomes `CRS.from_epsg(31287)` and `CRS.from_user_input(<31287's GDAL WKT1>)`.
- **Red test 2, `TestTransformLabel`**: the "none" pair becomes 31287's GDAL WKT1 and EPSG:31287, both orders; *red* `transform_label(<Lambert by parameters>, "EPSG:31287")` equals `transform_description` of that pair and is not "none".
- **Red test 6, `test_catchment_batch.py`**: `check_reach_crs(<31287's GDAL WKT1>, repo)` is None; *red* the Lambert by parameters is refused with today's wording (`the river file's CRS, <the string>, is not the DEM's, EPSG:31287`); the Lambert at 13.5 E stays refused.
- **Red tests 7 to 11**: every `proj4_of(n)` becomes `axes_swapped(n)` (in test 9 as the GeoJSON `crs` member; in test 11 as the member and as the flag). Each still passes with the point-moving methods refused. If a reader refuses WKT in one of those places, that is a departure to report, not to work around.
- ***Red*, `tests/python/test_domain.py`**, the leak at a site: `DomainPolygon(polygon=<a 1 by 1 degree box at 0-1 E, 50-51 N>, crs=<the +lon_0=10 string>).to_crs("EPSG:4326")` has its polygon's bounds 10 degrees further east (`minx` within 1e-9 of 10.0); at `5197f9a` it comes back unchanged.
- **Red test 12** stays: the None branch returns before the `noop` leg reads `definition`.

**After the round-1 red step** (`13b3e6d`), one line each:

- Ola's ruling, verbatim: "D15 b". The main session's option b: a PROJ string that names no datum is not its EPSG code (as on master), and the refusal adds a hint naming the code to write.
- The hint, appended to the refusal text: `; if you mean EPSG:n, write EPSG:n`, where `n` is the DEM's code.
- Two sites only, both in `catchment.py`: `check_reach_crs` and `delineate`'s reach check (refusal wordings above). Other refusing sites and every record stay as they are.
- When: `parse_crs(dem_crs).to_epsg(min_confidence=100)` is a code `n`, and the refused CRS's `to_epsg(min_confidence=100)` is None. A refused CRS with a code of its own (`EPSG:32633` against an EPSG:25833 DEM) gets no hint: its writer already chose a code.
- Ruled on the default to Ola's open question: the hint shows even when the refused CRS matches no code at all, so the Lambert at 13.5 E gets it too. The hint offers a code; it does not claim the two are the same.
- Ruled at code review round 2, intended under that reason, no code change: the hint's test is "the refused CRS has no exact code", so it also shows for `OGC:CRS84` against an EPSG:25833 DEM, and for `+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs` (which `same_crs` treats as EPSG:32633) against EPSG:25833; both get `; if you mean EPSG:25833, write EPSG:25833`, and `EPSG:32633` itself gets none (probe at `65cd528`: `_code_hint(g, "EPSG:25833")` for the three).
- Where: a private `_code_hint(given, dem_crs) -> str` in `catchment.py` (about 5 lines), returning the hint or `""`. Not in `crs.py`: both callers are in one module, and `crs.py`'s public surface stays `same_crs`, `transform_label`, `single_crs`. It runs only on the refusing path, so its identification calls (milliseconds) never touch an accepted run; `same_crs` has already parsed `given`, so it raises nothing new.
- `@tester`'s departure, accepted: test 6's refusal of the Lambert at 13.5 E, and `test_a_reach_in_another_crs_is_refused`, assert `startswith(today's words)`, not equality, since the hint may follow; the hint is pinned by `in` at the two no-datum tests.
- `@tester`'s departure, accepted: `crs_fixtures.axes_swapped` reverses the axis list of the code's PROJJSON with its `id` removed, then writes WKT2; it asserts the ID is gone and that pyproj's `==` calls the result different, so a no-op swap fails loudly.
- `@tester`'s departure, accepted: the `+lon_0=10` domain test bounds all four of the box's coordinates within 1e-9 degrees (the design named `minx` only); 1e-9 degrees is about 0.1 mm at longitudes up to 11, far inside the 10-degree move it detects.

### Net production lines

**Measured at `5197f9a`: +36** (`python3 tools/count_loc.py 44fa7f5
5197f9a`: 77 added, 41 removed), against the design's about +16. The +20
is all line-splitting the design did not cost, none of it new logic:

- `crs.py` +40 against about +29. `_same_frame` and `FRAME_TOLERANCE` came
  to 16 against about 7: the formatter splits the axis loop's `if` over
  three lines and the prime-meridian list over six. `same_crs` came to 13
  against 12, `single_crs` to 8 against 7; the two import lines changed in
  place (net 0).
- The sites -4 against about -13. `dem_input.py`'s `crs` import gains two
  names and passes 100 characters, so the formatter writes it as nine
  lines: `dem_input.py` is +4, not -4. `cli.py` is -5 against -6;
  `catchment.py` -3 as designed; the other four files 0.

**After code review round 1: about +20**, the D15 b hint's `_code_hint`
and its two call sites (about +5) included. The rule above drops
`_same_frame` (15 lines), the two constants (2), and 4 lines of
`same_crs` (the code sets, the separate `equals` return, and the parse
`_transformer` already does). The `dem_input.py` import may be packed onto
fewer lines under `# fmt: skip` only if the review says why (`CLAUDE.md`
section 2). Each of the thirteen comparisons is one line before and after,
so they save nothing; the audit's about -25 assumed `same_crs` was a
one-line alias of `!=`. Tests: red tests 1-7 came to 254 non-blank lines
added and 9 removed (`git diff -U0 b63132e 29aff00 -- tests`), 8-12 add
more, and the round-1 rows about 40. Section 6's row B and its total move by
about +40 accordingly; the drift point (one CRS rule) is still written once.

**Measured at `65cd528`: +22** (`python3 tools/count_loc.py 44fa7f5
65cd528`: 65 added, 43 removed), against about +20. `crs.py` +19 as
designed; `catchment.py` +4: -3 for `single_crs` as designed, and +7 for
the hint against about +5, the two `hint = _code_hint(...)` lines at the
call sites not costed; `cli.py` -5; `dem_input.py` +4 (its nine-line
import, not packed); the other four files 0.

### Citations this PR moves, pinned now

The `cli.py` edit moves every later line up about 6, `catchment.py` about
3, `dem_input.py` about 2. `feature_input.py` and `fetch/run.py` keep their
line count. Unpinned citations into those files, at or after the
edited lines, whose quotations hold at `44fa7f5`, were pinned to `44fa7f5`
in this design's commit, so the green commit breaks none:
`docs/increments/15c-geographic-dem.md` line 931 (`catchment.py:194-195`),
`docs/increments/15f-edge-strip.md` line 591 (`cli.py:1612`),
`docs/increments/25-plain-output.md` lines 53, 273, 274, 397 and 665
(`cli.py:1966-1969`, `:1040`, `:961`, `:1113-1125`, `:1145-1147`),
`docs/increments/29-nve-reference-catchments.md` lines 538 and 3130
(`dem_input.py:248`), 3134 (`cli.py:1956`) and 3380 (`catchment.py:259`), and
`docs/benchmarks/2026-10-05/nve-hrd/README.md` line 53 (`catchment.py:259`).
Left alone: unpinned citations whose text had already moved before this PR
(`15f-edge-strip.md` lines 1539 and 1904, `24-release-hardening.md` lines 515
and 517, `29-nve-reference-catchments.md` lines 3098 and 3114). They are
dated records of an older revision, and pinning them to `44fa7f5` would pin
a wrong line. No citation points into the test files this PR edits.

### Questions for Ola (defaults hold until he answers)

1. **Read GeoTIFFs whose CRS is given by parameters, like the Austrian
   openDEM file?** Default: yes, as a later small PR after B. The reader
   would build the CRS from the GeoKeys on the datum they name (MGI,
   EPSG:4312, for the Austrian file), not as a PROJ string, which names no
   datum and so is never the same as an EPSG code under B's rule. B's rule
   does not then accept it (the earlier probe kept EPSG's name in the WKT2,
   so PROJ matched by name; named "unknown", the pipeline is an
   inverse-then-forward Lambert, not a `noop`). The later PR matches by its
   own rule, same datum and PROJ equivalence on the east-north pair:
   `docs/increments/geotiff-crs-by-parameters.md` on branch
   `worktree-geotiff-param-crs`. One that matches no code stays refused as today. It overturns increment 11's rulings 6 and 7
   ("the CRS is resolved only through `pyproj.CRS.from_epsg`"), hence the
   question.
2. **The new wording** for DEM files in more than one CRS (table above).
   Default: as written.
3. **A CRS that is an EPSG code's definition under another spelling counts
   as that CRS:** no resampling, no refusal, transform "none" in the record.
   Default: yes. The limits it pins: a PROJ string that names no datum
   (`+ellps=GRS80` without `+datum`, as pyproj writes EPSG:25833) is not
   that CRS, so with it `--out-crs` resamples and a river file is refused,
   as on master; the `+towgs84` spelling is not the same; nor is one with
   an axis pointing the other way, another prime meridian, or a `+lon_0`
   on a `longlat`; CRS84 is the same as EPSG:4326. **Answered by Ola,
   "D15 b"**: the no-datum limit stands, and the river-file and reach
   refusals add `; if you mean EPSG:n, write EPSG:n` ("After the round-1
   red step").

## 10. T1 design: one way to drive the CLI in the tests

Tests only, so no red step and no `@developer`: `@tester` writes it,
`@reviewer` audits it, and it must pass on its base as well as after. The
branch starts at T2's head `b63132e` and lands after T2. Every citation in
this section into a test file is pinned, to `b63132e` unless it says otherwise.

### The rule

**Only helpers move; no test's assertion changes.** A test body changes only
where a helper it calls is renamed or takes the subcommand as its first
argument. The checks at the end make this mechanical.

### Re-measured at `b63132e`

X1's counts were taken at `12dace7`. Again, by an AST comparison of
top-level definitions over `tests/python/` (docstrings ignored):

| Helper | Copies | Distinct bodies | How they differ |
|---|---|---|---|
| `bumpy` | 9 (8 module-level, 1 inside a class) | 4 | only the seed: 16 (5 copies), 20 (2), 14 (1), 17 (1); all 17 x 21, uniform 0 to 50 m, `bumpy.tif` |
| `invoke` | 8 | 4 | the whole command line (3), `mesh` prepended (3), `draw` prepended and the raw `Result` returned (1), the whole command line with the output not passed through `plain` (1) |
| `run` (CLI) | 7 | 6 | two are identical; the others differ in return type and flags |
| `refused` (CLI, in-process) | 4 | 4 | see the mapping below |
| `squashed` | 3 | 2 | one passes the text through `plain` first |
| `square` / `plain_square` | 2 / 3 | 1 / 1 | the four other `square`s (`feature_fixtures.py`, `test_core_cdt.py`, `test_core_noding.py`, `test_domain.py`) are different things and stay |
| `write_geojson` | 3 | 3 | two write exactly what `test_cli_mesh_domain.geojson` writes, with a different `crs` default; `feature_fixtures.py`'s writes a FeatureCollection and stays |
| `plain`, `ANSI`, `BOX` | 2, 3, 2 | 1 each | none |
| module-level `CliRunner` | 17 | 2 | 16 pass `env={"NO_COLOR": "1", "TERM": "dumb"}`; `test_cli.py`'s passes nothing and stays |
| `USAGE = 2`, `ROWS, COLS = 17, 21`, `Ring` | 4, 5, 4 | 1 each | `test_dem_input_domain.py`'s `Ring` is not a CLI suite and stays |

The copies are also reached sideways: 39 import lines in 21 test files
import these names from other test modules (`test_cli_mesh`,
`test_cli_mesh_dem`, `test_cli_mesh_domain`, `test_cli_mesh_mosaic`,
`test_cli_catchment`, `test_cli_mesh_geographic`), so a test module is
also a helper library for six others.

The copies to delete hold 174 non-blank lines (measured per definition).
X1's "about 350 test lines" assumed the 38 inline
`code, output = invoke(...)` / `assert code == 0, output` pairs would be
folded into `ran`; this design does not do that (it would move assertions
out of test bodies), and X1 costed neither the new module nor the import
lines.

### The module: `tests/python/cli_driver.py`

Named `cli_driver`, not X1's `cli_harness`: since h17 (on branch
`worktree-ci-speed`), "harness test" means a test of `tools/` or the hooks,
selected by a `harness` marker, and this module serves the product suites.

Contents. A definition that moves unchanged keeps its name, so most call
sites do not change:

- `runner` (the 16 identical `CliRunner`s), `USAGE = 2`, `ANSI`, `BOX`,
  `plain(text)` (one docstring for the two copies), `squashed(text)` (the
  whitespace-only one).
- From `test_cli_mesh_domain.py`, moved: `UTM33`, `Ring`, `ROWS, COLS`,
  `SQUARE`, `HOLE` with their comments, and `geojson(path, outer, holes=(),
  crs=UTM33)`. From `test_cli_mesh_dem.py`: `write_tiff`.
- `COMMANDS`: the subcommand names, read once from the app
  (`typer.main.get_command(app).commands`; today `catchment`, `draw`,
  `fetch`, `fetch-stations`, `mesh`, `palette`, `station-catchments`,
  `version`).
- `invoke(command, *args) -> tuple[int, str]`: asserts `command in
  COMMANDS` (the message lists them), then runs `rasputin <command>
  <args>` and returns the exit code and `plain(output)`. The assertion is
  there because the mesh-only copies took no subcommand: a call written
  for them, such as `invoke("--dem", ...)`, would otherwise run `rasputin
  --dem ...`, exit 2 with "No such option", and pass a test that checks
  only for exit code 2. Probe, run with the main checkout's venv on a scratch copy
  of the module: `invoke("version")` returns `(0, ...)` and
  `invoke("--dem", "x")` raises the assertion naming the eight commands.
- `ran(command, *args) -> str`: `invoke`, then `assert code == 0, output`;
  returns the output. Used only inside helpers here.
- `refused(tmp_path, command, *args, says, squash=False) -> str`: the body
  of `tests/python/test_cli_mesh_domain_crs.py@b63132e:135-144`, with `command` passed
  through and, for each word, `assert squashed(word) in squashed(output)`
  when `squash` is set (the body of
  `tests/python/test_cli_mesh_geographic.py@b63132e:149-158`) and `assert word in
  output` otherwise. Both assert lines are kept verbatim, as two branches,
  so the assertion texts do not change.
- `mesh_to_vtk(tmp_path, *args, out="x.vtk") -> VtkFile`: the two
  identical `run`s (`tests/python/test_cli_mesh_domain_crs.py@b63132e:122-127`,
  `tests/python/test_cli_mesh_geographic.py@b63132e:132-137`), built on `ran("mesh", ...)`.
- Two fixture factories. Each returns a pytest fixture, and a suite binds it
  to the name its tests already request:
  - `rough_dem(seed)`: writes `tmp_path / "bumpy.tif"`, the array
    `default_rng(seed).uniform(0.0, 50.0, (ROWS, COLS)).astype(np.float32)`
    through `micro_tiff`. A suite writes `bumpy = rough_dem(16)`.
  - `polygon_file(outer, holes=())`: writes `tmp_path / "square.geojson"`
    with `geojson`. `square = polygon_file(SQUARE, (HOLE,))`,
    `plain_square = polygon_file(SQUARE)`.

  Probe of the pattern under this repo's pytest (9.1.1), in a scratch
  directory: two module-level bindings of one factory register under the
  names they are bound to (`pytest --fixtures` lists both), a test in a
  class sees the module-level one, a module that binds neither does not
  see it (its test errors, fixture not found), and `ruff check` and
  `ruff format --check` pass with the repo's settings. Importing a fixture
  function by name instead would leave an import that ruff reports unused.

A draft of the module is about 105 non-blank lines, about 30 of them moved
from `test_cli_mesh_domain.py` and `test_cli_mesh_dem.py`.

### How each copy maps

| Copy (at `b63132e`) | Becomes |
|---|---|
| `runner` in `tests/python/test_cli_catchment.py@b63132e:86`, `tests/python/test_cli_draw.py@b63132e:107`, `tests/python/test_cli_fetch.py@b63132e:61`, `tests/python/test_cli_mesh.py@b63132e:49`, `tests/python/test_cli_mesh_dem.py@b63132e:47`, `tests/python/test_cli_mesh_edge_strip.py@b63132e:103`, `tests/python/test_cli_mesh_landcover.py@b63132e:63`, `tests/python/test_cli_mesh_mosaic.py@b63132e:54`, `tests/python/test_cli_mesh_plain_output.py@b63132e:60`, `tests/python/test_cli_mesh_stats.py@b63132e:42`, `tests/python/test_cli_mesh_vtk.py@b63132e:40`, `tests/python/test_cli_station_catchments.py@b63132e:64`, `tests/python/test_fetch_nve.py@b63132e:78`, `tests/python/test_hardening.py@b63132e:50`, `tests/python/test_io_vtk_readback.py@b63132e:82`, `tests/python/test_palettes.py@b63132e:40` | `from cli_driver import runner`; every `runner.invoke(...)` unchanged |
| `plain`, `ANSI`, `BOX` in `tests/python/test_cli_mesh.py@b63132e:51-63`, `tests/python/test_cli_draw.py@b63132e:109-122`; `ANSI` in `tests/python/test_cli_mesh_geographic.py@b63132e:836` | imported from `cli_driver`; the 10 suites that import `plain` from `test_cli_mesh` import it from `cli_driver` |
| `invoke`, whole command line: `tests/python/test_cli_catchment.py@b63132e:93-95`, `tests/python/test_cli_fetch.py@b63132e:64-66`, `tests/python/test_fetch_nve.py@b63132e:81-83` | `cli_driver.invoke`; call sites unchanged (each already passes the subcommand first). `test_cli_mesh_geographic.py`'s `invoke_any` alias becomes `invoke` |
| `invoke`, mesh: `tests/python/test_cli_mesh_dem.py@b63132e:53-55`, `tests/python/test_cli_mesh_mosaic.py@b63132e:59-61`, `tests/python/test_cli_mesh_vtk.py@b63132e:48-50`, and its importers | `cli_driver.invoke`; every call gains `"mesh"` as its first argument (about 75 calls in 15 files, three of them inside asserts in `test_cli_mesh_vtk.py`) |
| `invoke` in `tests/python/test_cli_draw.py@b63132e:125-126` (returns the `Result`, about 45 calls) and `tests/python/test_cli_station_catchments.py@b63132e:76-78` (output not passed through `plain`; the suite splits it into lines) | stay, on `cli_driver.runner`: changing their return would change their tests |
| `bumpy`: seed 16 in `tests/python/test_cli_mesh_domain.py@b63132e:133-136`, `tests/python/test_cli_mesh_features.py@b63132e:125-128`, `tests/python/test_cli_mesh_landcover.py@b63132e:115-118`, `tests/python/test_cli_mesh_multi_features.py@b63132e:89-92`, and in class `TestTheSameCrs` at `tests/python/test_cli_mesh_domain_crs.py@b63132e:283-286`; seed 20 in `tests/python/test_cli_constraint_feet.py@b63132e:60-63`, `tests/python/test_cli_start_quality.py@b63132e:51-54`; seed 14 in `tests/python/test_cli_mesh_refine.py@b63132e:97-101`; seed 17 in `tests/python/test_cli_mesh_stats.py@b63132e:51-54` | `bumpy = rough_dem(<seed>)` at module level. In `test_cli_mesh_domain_crs.py` the binding goes at module level (no other test there requests `bumpy`): bound inside a class, a factory-made fixture receives the instance as `tmp_path` (probed: the test errors with `TypeError` on `TestInClass / str`) |
| `square` in `tests/python/test_cli_mesh_domain.py@b63132e:139-141`, `tests/python/test_cli_mesh_stats.py@b63132e:57-59`; `plain_square` in `tests/python/test_cli_mesh_features.py@b63132e:131-133`, `tests/python/test_cli_mesh_landcover.py@b63132e:121-123`, `tests/python/test_cli_mesh_multi_features.py@b63132e:95-97` | `square = polygon_file(SQUARE, (HOLE,))`, `plain_square = polygon_file(SQUARE)` |
| `squashed` in `tests/python/test_cli_mesh_cache.py@b63132e:48-50`, `tests/python/test_cli_mesh_geographic.py@b63132e:113-115` | imported from `cli_driver`. `tests/python/test_cli_station_catchments.py@b63132e:379-383` (through `plain` first) stays |
| `write_geojson` in `tests/python/test_cli_mesh_domain_crs.py@b63132e:104-111` (default `crs=None`) and `tests/python/test_cli_mesh_geographic.py@b63132e:118-125` (`crs` required, one ring); the alias `utm33_geojson` in `test_cli_mesh_domain_crs.py`; the alias `write_domain` in `test_cli_mesh_edge_strip.py` | `geojson`. A call that relied on `crs=None` passes `crs=None`; a positional `crs` becomes `crs=...`, because `geojson`'s third parameter is `holes`. Same dict, same key order, so the same bytes |
| `refused` in `tests/python/test_cli_mesh_domain_crs.py@b63132e:135-144` and `tests/python/test_cli_mesh_geographic.py@b63132e:149-158` | `cli_driver.refused(tmp_path, "mesh", *args, says=...)`, with `squash=True` in `test_cli_mesh_geographic.py` |
| `refused` in `tests/python/test_cli_catchment.py@b63132e:137-144` (any non-zero exit, not 2) and `tests/python/test_cli_mesh_plain_output.py@b63132e:108-115` (sets `sys.argv`, no `--out`) | stay: folding them in would change what they assert |
| `run` in `tests/python/test_cli_mesh_domain_crs.py@b63132e:122-127`, `tests/python/test_cli_mesh_geographic.py@b63132e:132-137` | `mesh_to_vtk` (31 calls renamed) |
| `run` in `tests/python/test_cli_catchment.py@b63132e:130-134`, `tests/python/test_cli_mesh_refine.py@b63132e:80-84` | stay, with their first two lines replaced by `ran(...)` |
| `run` in `tests/python/test_cli_mesh_landcover.py@b63132e:126-133`, `tests/python/test_cli_constraint_feet.py@b63132e:85-89`, `tests/python/test_cli_start_quality.py@b63132e:78-82` | stay (different returns; the last two spy on `cli.refine`, see "Not taken") |
| `USAGE` in `tests/python/test_cli_mesh_dem.py@b63132e:50`, `tests/python/test_cli_mesh_mosaic.py@b63132e:55`, `tests/python/test_cli_mesh_edge_strip.py@b63132e:104`, `tests/python/test_cli_mesh_plain_output.py@b63132e:61`; `ROWS, COLS` in `tests/python/test_cli_constraint_feet.py@b63132e:49`, `tests/python/test_cli_mesh_plain_output.py@b63132e:62`, `tests/python/test_cli_mesh_stats.py@b63132e:48`, `tests/python/test_cli_start_quality.py@b63132e:46`; `Ring` in `tests/python/test_cli_mesh_domain_crs.py@b63132e:88`, `tests/python/test_cli_mesh_geographic.py@b63132e:104` | imported from `cli_driver` |

After the change, no test module imports any of these names from another
test module. Names that are not helpers of the command (`quarter_circle`,
`SENTINEL`, `file_field`, `write_tiles`, `same_mesh`, `terrain`, ...) keep
their homes; moving them is not this PR.

Imports left unused (`CliRunner`, `re`, `json`, `app` where only the
deleted helper used it) go; `ruff check` finds them. Prose in test
docstrings and comments that names a moved helper's old home (for example
`tests/python/test_cli_mesh.py@b63132e:56`, "per `test_cli_draw.py`'s `plain`") is
reworded.

`tests/python/test_cli_mesh_geographic.py@b63132e:106` keeps its own `ROWS = COLS = 60` (another grid) and
does not import `cli_driver`'s.

### Taken from the test audit (`docs/increments/test-audit.md` on `worktree-ci-speed`)

- **R10, the GIL ticker:** `_Ticker` and `ticks_during` are identical in
  `tests/python/test_core_cdt.py@b63132e:751-790` and `tests/python/test_core_noding.py@b63132e:721-757`
  (64 non-blank lines). They move to `tests/python/gil_probe.py`; both
  suites import `ticks_during`. Each suite's thresholds and its control test
  stay as they are: the two suites' constants differ, and merging the two
  control tests into one would remove a test.
- **Not taken from R10:** the `valley` pair in `test_core_accumulate.py`
  and `test_core_upstream.py` (about -7 net once it has a module of its
  own); the seven identical two-line `repo` fixtures (each binding would
  save one line, and those files are the harness tests that h17 is
  re-marking on another branch); the viz fakes shared by
  `test_viz_scene.py` and `test_viz_svg.py` (they are half of a hand-drawn
  fixture whose vertex arrays differ between the two suites; moving half
  separates the drawing from its data); the lazy `importlib` imports in the
  five CLI suites this PR touches (each returns a module from a fixture
  that has the module's own name, so a plain import would clash with the
  fixture and the bodies would change); the shared quarter-circle run (it
  changes test structure).
- **Not taken, R1 (`tests/python/seams.py`):** it re-points what the tests
  patch (`cli.refine`, `cli._engine`), which is a change to what they test,
  and it needs the R9 audit of `raising=False` beside it. The two identical
  `calls`/`box` refine spies in `tests/python/test_cli_constraint_feet.py@b63132e:66-82`
  and `tests/python/test_cli_start_quality.py@b63132e:57-75` are R1's and wait for it
  (Q3, `audit-test-seams`).

### Expected test-line delta

Non-blank lines (the same count as section 8's): about -100. The CLI part
is about -80: 174 lines of copies and about 27 unused import lines go;
the new module adds about 75 beyond what it moves, the per-file
`cli_driver` imports and factory bindings about 25, and the call-site
edits about 25 (adding `"mesh"` to about 75 calls and renaming 31 pushes
some lines past 100 characters, and ruff format decides how they wrap).
The ticker is about -20. The spread is in the call-site wraps. Production
lines: 0 (`python3 tools/count_loc.py b63132e HEAD`).

Measured at `2816d41`: +457 −475, −18 net against about −100, because
adding `"mesh"` pushes calls past 100 characters and ruff format then puts
one argument per line, and `cli_driver.py` is 123 lines against about 105.

### Checks, run by `@tester` before committing and reported in the commit message

1. **Same tests.** `pytest --collect-only -q tests/python` at `b63132e`
   and at the head: the sorted node-ID lists are identical (`diff` empty).
   Baseline at `b63132e`, with the main checkout's venv: 5264 collected in
   the whole suite, 1123 in the 28 files this PR edits. The full suite's
   passed and skipped counts are equal on both, on the same machine.
2. **Same assertions.** A scratch script (not committed) takes every
   `assert` statement in `tests/python/*.py` at both revisions,
   `ast.unparse`d, after replacing `invoke('mesh', ` by `invoke(` in the
   head's text (the only rename inside an assert). For each `test_*.py`
   file, the multiset of asserts inside test functions is unchanged; and
   across the whole directory, the set of distinct assert texts is
   unchanged (the asserts of deleted copies are found in `cli_driver.py`,
   `gil_probe.py` or a copy that stays).
3. **Same inputs.** For seeds 14, 16, 17 and 20 the old `bumpy` and
   `rough_dem` write byte-identical `bumpy.tif` files; one call of each
   replaced `write_geojson` and its `geojson` replacement write identical
   bytes. Hashes in the commit message.
4. **The guard fails.** In a scratch copy, one mesh call left without
   `"mesh"` fails with the assertion naming the subcommands; restored
   afterwards. Name the plant in the commit message.
5. **No sideways imports.** `grep -nE '^from test_' tests/python/*.py`
   shows none of the moved names.
6. `ruff check .`, `ruff format --check .`, `python3
   tools/check_citations.py` (clean, and its at-risk list re-read), and
   `count_loc.py` at 0.

No mutation round beyond check 4.

### Citations into the edited files

Found with a scratch script that runs `tools/check_citations.py`'s own
scanner over every scanned file and keeps the unpinned line citations
whose target is one of the 28 files or the two `valley` files first
considered. Fourteen, all in `docs/increments/`;
none in code comments. Thirteen are pinned in this design's commit, each
to the commit that wrote the citing line (`git blame` on that line), where
the cited lines read as quoted. The exception: the code-review record at
`docs/increments/25-plain-output.md:944` names the design citation that its
branch had moved, so it is pinned with the design's, at `55c043e`.

| Citing line | Pinned to |
|---|---|
| `docs/increments/15e-memory-fixes.md:328` | `tests/python/test_cli_mesh_geographic.py@3e01580:327-328` |
| `docs/increments/15e-memory-fixes.md:386` | `tests/python/test_cli_mesh_multi_features.py@97eea35:17` |
| `docs/increments/25-plain-output.md:288`, `:433` | `tests/python/test_cli_mesh_refine.py@a2d3319:151-155` |
| `docs/increments/25-plain-output.md:585` | `tests/python/test_cli_mesh.py@ef52e8e:302` |
| `docs/increments/25-plain-output.md:587` | `tests/python/test_cli_mesh_landcover.py@ef52e8e:13` (and its `:62-66`) |
| `docs/increments/25-plain-output.md:590` | `tests/python/test_cli_mesh_stats.py@ef52e8e:262` (and its `:285-304`) |
| `docs/increments/25-plain-output.md:830` | `tests/python/test_cli_mesh_mosaic.py@5e235eb:278` (and its `:302`, `:331-336`) |
| `docs/increments/25-plain-output.md:942`, `:944` | `tests/python/test_cli_mesh_refine.py@55c043e:133` |
| `docs/increments/27-node-sampling.md:315` | `tests/python/test_cli_mesh_plain_output.py@80990e0:279` |
| `docs/increments/h16-harness-fixes.md@abd7c68:579` | `tests/python/test_cli_mesh_geographic.py@97eea35:886` |
| `docs/increments/h16-harness-fixes.md@abd7c68:615` | `tests/python/test_hardening.py@a61e848:247` |

Several of these were already stale at `b63132e` (for example
`tests/python/test_cli_mesh_landcover.py@b63132e:13` no longer quoted the line the record
discusses); pinning to the writing commit fixes those too. The fourteenth,
`docs/increments/29-nve-reference-catchments.md:3092` into
`tests/python/test_core_accumulate.py@b63132e:22`, points at a file this PR no longer edits.
`@tester` reruns the scan after the change: any new unpinned citation into
an edited file is `@tester`'s to pin in the same commit; one in a doc is
reported to the main session for `@architect`.

### Merge-order hazard

23c-2 (`worktree-23c`, not merged) adds suites that import `invoke` and
`USAGE` from `test_cli_mesh_dem` (`pieces_fixtures.py` and four
`test_cli_mesh_pieces_*.py`) and `write_geojson` from
`test_cli_mesh_geographic` (`test_cli_mesh_pieces_reprojected.py`).
Whichever of T1 and 23c-2 merges second re-points them to `cli_driver`,
adds `"mesh"` to its calls, and passes `crs=` by keyword. If it misses an
import, the import fails; if it misses a `"mesh"`, `invoke`'s assertion
fails, except before an argument that is itself a subcommand name (such as
the `catchment` gallery fixture), where the test's own message check fails
instead. Neither can pass silently, and the merge queue runs the second on
top of the first.

## 11. PR F design: `audit-catchment-shared` (F4; F10's catchment types; F12's three edges)

Branch `worktree-audit-catchment`, from T2's head `b63132e`. The section
number follows PR B's design (section 9, merged in #192) and T1's (section
10, merged in #189); F lands after both, and its merge of master checks the
numbers. The other open branch, PR C (`worktree-audit-geojson`), also takes
section 11 in its own master merge; F lands first (it is ready), so C takes
section 12 in its next master merge, and F's text stays 11. `src_python/` at `b63132e` is byte-identical to
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
`io/rivers.py` (its table, section 9, merged in #192); the overlap
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

T2 code review r3 (90cba64..b63132e): APPROVED.

T1 code review r1 (b63132e..2816d41): CHANGES REQUESTED — two unpinned citations in §9 (now section 10); fixed in 0d63d00.

**PR B (`audit-crs-helpers`), code review, round 1, 2026-10-05.** Head `5197f9a`. Verdict: CHANGES REQUESTED. LOC: +36 net production (`count_loc.py 44fa7f5 5197f9a`). Blocking, all `@architect`: (1) `same_crs` rule 2 calls `+proj=longlat +datum=WGS84 +lon_0=10` the same as EPSG:4326 (also `+lon_0=-3`, and `+datum=NAD83 +lon_0=10` against 4269) while the always_xy transform moves every point 10 degrees, and the rule is not transitive (EPSG:4326 = CRS84, the string = 4326, the string is not CRS84); (2) the status paragraph is stale; (3) the +36 against about +16 is not reconciled. Ruled in section 9, "After code review round 1".

**PR B (`audit-crs-helpers`), code review, round 2, 2026-10-05.** Head `65cd528`. Verdict: CHANGES REQUESTED, prose only. LOC: +22 net production (`count_loc.py 44fa7f5 65cd528`). Blocking, all `@architect`: the stale status paragraph, the unrecorded +22, section 6's row B; design note: the hint also shows for `OGC:CRS84` and a datum-WGS84 UTM string against an EPSG DEM. Fixed in this file (status, section 6, section 9 "Net production lines" and "After the round-1 red step"); no code change.

**PR F (`audit-catchment-shared`), code review, round 1, 2026-10-05.** PR F code review r1 (`b63132e..5c6a9f9`): CHANGES REQUESTED — `project_structure.md` rows (`project_structure.md@5c6a9f9:150-153` still named `_core.upstream` and `_core.reduce_ring`, `:251` listed `RiverSegment` under `rivers.py`, no rows for `hydrography.py` and `catchment_core.py`) and section 5's picture, the status line (`docs/increments/python-audit.md@5c6a9f9:25-27`), and section 11's net lines against the measured +29; fixed by `@architect` in the commit that records this round.

Non-blocking, for a later `@tester` and `@developer` pass: tests still reach `RiverSegment`, `Station` and `Lake` through the codec modules rather than `tin_engine.hydrography` (`tests/python/test_station_set.py@5c6a9f9:51, 245`, `tests/python/test_rivers.py@5c6a9f9:89`); `catchment.py`'s module docstring still names `_core.upstream`, `_core.reduce_ring` and `_core.accumulate` where the calls now go through `catchment_core`.

PR F code review r2 (5c6a9f9..e2baa5f): APPROVED.

PR F master merge review (e2baa5f..e4f7a42): CHANGES REQUESTED — resolution correct (+29 unchanged; layering rows and four UPWARD true; r2 APPROVED confirmed); `docs/increments/python-audit.md@e4f7a42:894-896` (B 'becomes section 10') and `:35-37` (status names @tester next; #190 missing) wrong; fixed in the next master merge.

PR F second master merge review (e4f7a42..8c2bc9a): APPROVED — resolution is only `catchment.py`'s import block (catchment_core and GaugePath, plus master's same_crs and single_crs, no _core) and `python-audit.md` (status, section 6 total, section 11 numbering, review records); +29 unchanged; both blockers from the first merge review closed; headings 8, 9 B, 10 T1, 11 F unique, and every section-N citation in src, tests and probes resolves; a planted `_core` import fails test_layering.
