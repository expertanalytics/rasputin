# Increment 54: the open core reads as one tool

**Status: designed** (`@architect`, 2026-10-16), one public PR, branch
`worktree-open-core-cleanup`. Not reviewed yet.

## 1. What Ola asked

Ola, 2026-10-15: "Tomorrow, we should look at the open core, and ensure it's a
logically piece of software." On 2026-10-16, on starting the five ranked fixes
of the read-only review as one public PR: "3: yes". His rule for the round:
"The best lines of codes are the ones you never built because you understood
they were not needed. The second best lines of codes are the ones you deleted
for the same reason :)".

The review is private, outside git
(`rasputin_scratch/open-core-review/2026-10-16.md`, master `9a7c312d`). Its
finding in one line: the split of increment 51 is clean, but a newcomer meets
developer machinery first (a gallery of test shapes inside `mesh`, design-record
labels in `--help`), and the plugin hook is documented only in an internal
design record. This increment takes its five ranked fixes, in one PR.

**Ola's rulings of 2026-10-16, recorded here:**

- The groundwork for cutting a huge domain into pieces (increment 23's cut
  runs) **stays**: basin scale is the standing fallback for unattended nights.
  Its docstrings say it is not wired to a command yet (fix 5).
- Renaming the import package `tin_engine` to `rasputin` is **decided later**,
  before the plugin hook is announced. Not in this increment; `docs/plugins.md`
  is written against `tin_engine` and says nothing of a rename.

## 2. Prior art: legacy and literature

Written before the design.

*Literature.* No algorithm here; the subject is command-line help and plugin
documentation.

- The Command Line Interface Guidelines (clig.dev, *Help*): "Display the most
  common flags and commands at the start of the help text", and the help
  should include "A description of what your program does". Their *Output*
  section: "By default, don't output information that's only understandable by
  the creators of the software." Fixes 1, 3 and 4 apply these as written. The
  guidelines have no rule on grouping commands; git's help groups them under
  headings, and that is the shape used here (Typer's `rich_help_panel`).
- Plugin discovery by entry points (Python Packaging User Guide, *Creating and
  discovering plugins*) and the CMake side (`cmake.root`, scikit-build-core's
  *Search paths*) are already cited in `docs/increments/51-open-core.md`,
  *2.1 Practice: how Python and C++ projects ship plugins*. `docs/plugins.md`
  documents that design; it changes none of it.

No novelty is claimed.

*Legacy.* The legacy tree had argparse scripts, one per task, with short help
strings and no grouping, no plugins, and no gallery. Nothing to carry across.

```
$ git grep -l -E 'argparse|click|typer|add_argument' legacy-archive -- legacy
legacy-archive:legacy/rasputin/application.py
legacy-archive:legacy/rasputin/geo_tiff_reader.py
legacy-archive:legacy/rasputin/web_visualize.py
```

## 3. Ideas file

Taken up (they are help strings this increment rewrites anyway):

- "`--start-quality-gain` help: say a negative gain turns the line split off":
  add ", with no line split even with constraint feet" to that option's help.
- "`--bbox` help: 'its nodes' is ambiguous with `--out-crs`": say "snapped
  outward to the DEM's nodes".

Overlapping, left out (handed back): the history checker's gaps (a capitalised
"Ruling", K and W IDs, bare "§N"), which `@reviewer` covers by eye in
`cli.py` for this PR; `app(prog_name="rasputin")` for `python -m tin_engine`'s
usage line.

## 4. The design, fix by fix

All production changes are in Python and Markdown. **No C++ file changes**, so
neither the C++ build nor `@perf`'s acceptance run is needed; the seam header
`include/terrain/refinement/seam.hpp` is deliberately left untouched for that
reason (fix 5).

### Fix 1. No history and no design notes in the help text

**Scope.** `src_python/tin_engine/cli.py`: every command docstring (they are
the help Typer prints and the MCP tool descriptions), every `help=` string,
every comment and the module docstring. At `9a7c312d`,
`python3 tools/check_code_history.py src_python/tin_engine/cli.py` prints 103
hits (measured), from line 7 to line 2912. After: 0, and `cli.py` is added to
`[tool.check_code_history] paths` in `pyproject.toml`, so the gate keeps it so.

**What a command's help is, after this fix:**

- First paragraph: one sentence saying what the command does, at most 120
  characters. It is also the MCP tool's description (`mcp/schema.py` takes the
  first paragraph), so it stands alone.
- Then at most a short paragraph each on what goes in and what comes out. No
  reStructuredText markup (Rich prints double backticks raw): write `--dem`,
  not ``` ``--dem`` ```. No module paths (`tin_engine.mosaic`), no design file
  names (`13-bundled-mesh.md`), no increment numbers, no "why it is not
  otherwise" (why there is no `--format`, MDAL's caveat, why `draw` exits 0 on
  a failure).
- `mesh`'s help (today about 60 lines, and its first line says "VTK or PLY"
  only) becomes at most 15: its first line names all three outputs, for
  example "Mesh a DEM into a TIN and write it as .vtk, .ply or .gpkg, chosen by
  the suffix of --out."

**Deleted, not rewritten.** The design rationale in the docstrings of `mesh`
and `draw`, the `cli.py` module docstring's account of the draw pipeline
(steps 1 to 5 and the `Attempt` argument; `Attempt`'s own docstring keeps what
a reader of that class needs), and the history half of each comment (the
"Increment 20, C2 (a):", "Grown once, at increment 8" kind). The design
records under `docs/increments/` already hold all of it. A comment whose only
content was the label goes entirely. Code is unchanged by this fix:
`python3 tools/check_code_history.py --same-code <base>` passes on the commit
that does it, which is therefore its own commit, before fixes 3 and 4.

### Fix 2. `docs/plugins.md`, for a third-party plugin author

**Scope.** A new `docs/plugins.md` (about 120 lines), generic: no private
package, repository or command name. A paragraph in `README.md`
(*Add a command of your own*, after *Connect an agent*) pointing to it. Three
docstring lines in `plugin_api.py`.

Its sections, each from the code it describes:

1. **Add a command.** Declare `[project.entry-points."rasputin.commands"]`
   with a `typer.Typer` as the value; each of its commands becomes a top-level
   `rasputin` command, and an MCP tool. The example is the hello plugin of
   `tests/python/plugin_fixtures.py`, as a package.
2. **Refusals.** A plugin that fails to import, or whose value is not a
   `typer.Typer`, is skipped with one line on stderr; a command with a core
   command's name is skipped (the core keeps it); of two plugins with one name,
   the first in distribution-name order wins. `rasputin version` lists the
   plugins that loaded. (`plugins.py`.)
3. **Where it appears in `--help`.** Under its own `rich_help_panel` if the
   command sets one, else under *Plugins* (fix 4).
4. **What you may import.** `tin_engine.plugin_api` and nothing else; every
   name in its `__all__`, one line each. "Promised" means: a name is removed or
   changed only in a release that says so; any other `tin_engine` module may
   change in any commit. The list is what plugins use today and grows on
   request. `plugin_api.py`'s module docstring gains that last sentence.
5. **C++.** The headers installed in the wheel; the promised ones
   (`terrain/core/point.hpp`, `terrain/raster/{geometry,raster,view,sample}.hpp`,
   `terrain/parallel_util/chunks.hpp`, `rasputin_bindings/numpy_raster.hpp`;
   the rest are installed but not promised); `find_package(rasputin CONFIG
   REQUIRED)` and the target `rasputin::terrain_headers`, found through the
   `cmake.root` entry point. The example is `tests/consumer/CMakeLists.txt`.
6. **The raster tag rule.** Declare your own `struct MyTag;` in your own
   source, bind `TaggedRasterView<MyTag>` with `py::module_local()`, build it
   from `core_args(tile)`, and never accept `_core.RasterView`: only NumPy
   arrays and Python scalars cross between two extension modules. One
   sentence of why, then a pointer to `docs/increments/51-open-core.md`,
   *4.3 The raster across two extension modules*.

Every claim in it is checked by `@reviewer` by running it
(`.claude/REQUIRED-READING.md`, *Before you publish*).

### Fix 3. The gallery leaves `mesh`; `draw` is hidden

**`mesh`: `--dem` becomes required.** `dem: list[str]` with no default, so a
missing `--dem` is Click's own "Missing option '--dem'" (exit 2). Today
`rasputin mesh --out x.vtk` answers "give a gallery fixture name or --dem PATH".

**Deleted from `cli.py`**, because with `--dem` required nothing reaches them:

- the positional `name` argument, `--flat` and `--crs`;
- the `else:` branch that writes a gallery fixture (its five "applies only
  with --dem" refusals, the `--flat` refusal, the `--crs` comment check), and
  `_fixture_mesh`;
- every other `if not dem` / "applies only with --dem" refusal in `mesh`
  (`--domain`, `--domain-crs`, `--bbox`), the `--flat`/`--crs` refusals under
  `if dem:`, the `if dem:` itself (its body unindented), the "exactly one of
  the two" check and the unknown-fixture check;
- the `has_dem`/`dem: bool` parameter of `_class_source`, `_tolerance_lines`
  and `_tolerance_slope` and the refusal each makes on it (`mesh` is their only
  caller: `grep -n` on the three names);
- the `None` cases of `dem_run` that become unreachable, where mypy then
  shows them dead;
- the now unused imports (`checked_ascii`, `flat_record`).

**Deleted from `run_record.py`:** `flat_record`, `FLAT_HEIGHTS` and the
"heights are not real" sentence in `summary`.

**`draw`** stays as a developer's bench (it is the only user of `viz/`), but
`hidden=True`: callable as `rasputin draw`, absent from `rasputin --help`.
`mcp/schema.py`'s `derive` skips hidden commands, so an agent is no longer
offered it (one condition beside `name != SERVER_COMMAND`).
Deleting `draw` and `viz/` is a question for Ola (section 7), not this PR.

**Docs.** `project_structure.md`: the `cli.py` entry loses "the gallery's
--crs goes through io/mesh_checks' checked_ascii", and the `viz/` entries say
`draw` is hidden. `testing.md` is unchanged (it describes the suites).

### Fix 4. Commands grouped and ordered; one description

**Registration in one place.** The `@app.command(...)` decorators in `cli.py`
are removed, and one block at the end of `cli.py` registers the commands in
the order `--help` shows them, each with its panel (constants, one per panel):

| Panel | Commands, in order |
|---|---|
| Meshing | `mesh`, `palette` |
| Catchments | `catchment`, `station-catchments` |
| Data | `fetch`, `fetch-stations` |
| Other | `version`, `mcp` |
| (hidden) | `draw` |

Rich prints panels in the order of their first command, and the unnamed
default panel (*Commands*) first of all (probed with Typer 0.27.2 on
2026-10-16). With every core command in a panel, only plugin commands could
land in *Commands*, at the top. So `plugins.application` gives a plugin
command that names no panel the panel *Plugins*, on a **copy** of its
`CommandInfo` (the plugin's own app is never changed, as today), and the
*Plugins* panel comes after the core's. The MCP tool list follows the new
order (it is the registration order).

**One description**, a constant `DESCRIPTION` in `cli.py`, used by
`typer.Typer(help=...)` and equal to `[project] description` in
`pyproject.toml`, taken from the README's first sentence: "Rasputin turns a
digital elevation model (DEM) into a triangle mesh of the terrain that stays
within a height tolerance you choose." (Today: "Parallel TIN engine for
terrain meshing." and "... CGAL-free, GDAL-free, pure C++20 core with async
Python bindings", which describes history and is wrong about async.) The
callback's docstring keeps only why the callback exists; Typer prints the
app's `help`, not the callback's (probed).

### Fix 5. Stale facts

- `INSTALL.md`, *Prerequisites*: drop "the compiler must provide the C++20
  `<chrono>` calendar types". `grep -rn -E 'year_month_day|sys_days' include
  src bindings` is empty at `9a7c312d`; `<chrono>` remains only for timing
  (`refinement/refine.hpp`, `refinement/refine_points.hpp`), which any C++20
  compiler has. `CLAUDE.md`'s "`std::chrono`" stays true for that reason and
  is not changed.
- `ROADMAP.md` row 51 and `docs/increments/51-open-core.md`'s status line: P2
  merged as #271 (`9a7c312d`); the private package's PRs are tracked in the
  private repository.
- The cut-run groundwork (Ola's ruling, section 1): one plain sentence that it
  is not yet called by any command, in the comment above
  `features.PIECE_VOCABULARY` (rewritten without its "Increment 23c" label),
  the module docstring of `io/mesh_index.py`, and `refine_seam`'s docstring in
  `_core.pyi`. Not in `seam.hpp` or `bindings/core.cpp`: a C++ change would
  bring a build and the acceptance rule for a comment.

## 5. Tests (`@tester`, red before green)

Fixtures: the CLI itself; `rough_dem` and `micro_tiff` of
`tests/python/cli_helpers.py` (17 rows by 21 columns, 357 nodes, 10 m by 5 m
spacing, EPSG:25833, z uniform in 0 to 50 m from a seed); the hello plugin of
`tests/python/plugin_fixtures.py`.

- **T1, gate scope.** `src_python/tin_engine/cli.py` is listed in
  `[tool.check_code_history] paths`, and the tool run on it exits 0 with no
  output. Red today: 103 hits.
- **T2, rendered help carries no history.** For the root and every command of
  `cli.app`, `draw` included, the rendered `--help` (wide columns) and each
  MCP tool description from `derive(cli.app)` match none of the patterns of
  `tools/check_code_history.py`'s `RULES` (imported from the tool by path, so
  there is one list), and contain no double backtick, no `tin_engine.` and no
  `.md`. This checks what a user sees, including text a `help=` builds from a
  constant, which the gate does not read.
- **T3, `mesh` without the gallery.** `rasputin mesh --out x.vtk` exits 2
  naming `--dem`, and its output has neither "gallery" nor "fixture";
  `rasputin mesh catchment --dem <rough> --out x.vtk` exits 2 (an extra
  argument); `--flat` and `--crs` are unknown options (exit 2); `mesh --help`
  contains none of "gallery", "fixture", "--flat", "--crs"; the first
  paragraph of `mesh --help` names `.vtk`, `.ply` and `.gpkg`.
- **T4, grouping and order.** Root help shows the panels Meshing, Catchments,
  Data, Other in that order, each with the commands of the table in that
  order; `draw` is absent from root help, `rasputin draw --list` still exits 0;
  `derive(cli.app)` names no `draw`. With the hello plugin passed to
  `application`, its command appears under *Plugins*, after *Other*, and the
  plugin app's own `CommandInfo` still names no panel afterwards; a plugin
  command that names a panel keeps it.
- **T5, one description.** Root help shows `cli.DESCRIPTION`, and it equals
  `pyproject.toml`'s `[project] description`.
- **T6, plugins page.** `docs/plugins.md` exists, names the group
  `rasputin.commands`, and names every entry of `plugin_api.__all__`, so a
  name added to the interface without a line in the page fails.
- **Moved, not new.** The tests that run `mesh <fixture> --flat` to exercise a
  writer (68 `--flat` occurrences in 19 files under `tests/python/`, counted at
  `9a7c312d`) move to `--dem` with `rough_dem`, or with `--domain` and
  `--features` where they need constraint edges with feature bits. Tests that
  pin only the deleted path go with it: `--flat`, `--crs`, the flat record and
  its sentence, the unknown fixture, the fixture with no mesh. Tests that call
  `cli._triangulated` or `cli.Attempt` directly (the mask tests in
  `test_cli_mesh.py`) stay: `draw` keeps that code. The MCP snapshots in
  `test_mcp_schema.py` change for the tool list and the rewritten help.

No suite here is invariant-critical: no mutation round.

## 6. Size against the 700-line rule (`CLAUDE.md` §2)

Docstrings and comments do not count; rewording a `help=` string adds and
removes the same lines. Estimated net production lines:

| Fix | Net |
|---|---|
| 1, help and comments | about 0 |
| 2, `plugins.md`, README, `plugin_api` docstring | 0 |
| 3, gallery out of `mesh`, `draw` hidden, `derive` skips hidden | about −80 |
| 4, registration block, panels, plugin panel, description | about +15 |
| 5, stale facts | 0 |
| **Total** | **about −65** |

One PR. Commit order on the branch: the red suite; fix 1 alone (the
`--same-code` check passes on it); then fixes 3, 4, 2 and 5.

## 7. Questions for Ola

1. Delete `rasputin draw` and the drawing package `viz/` (about 1,200 lines of
   synthetic test shapes and an SVG renderer, plus their tests) instead of
   hiding the command? **Default: hidden in this PR; deletion is a separate PR
   only on your yes.**
