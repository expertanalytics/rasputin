# Increment 51: open core, the snow and radiation code moved to a private package

**Status: designed** (`@architect`, 2026-10-15). Design review round 1
(changes requested) answered in this file. The five questions for Ola were
ruled at their defaults on 2026-10-15 (section 11).

## 1. What Ola asked

Ola owns Expert Analytics and has ruled that rasputin becomes *open core*: the
public repository `github.com/expertanalytics/rasputin` (MIT) keeps the
meshing core, and the snow and radiation modelling, with what is built on it,
moves into a package
that lives in the private repository `github.com/expertanalytics/rasputin-private`
and is not published.

Ola, verbatim, answering "keep the public repo public, move snow, radiation
and wind into a private package (open core)?": "yes".

What this increment delivers, in one line: a public core with no snow
or radiation code and a small, tested plugin interface; and a private
package, working name `rasputin-snow`, that builds against that interface
and makes `rasputin snow` appear when it is installed.

What it does not do: rewrite history. The snow pull requests merged on public
master (#257 to #267) stay in the public history and stay MIT for anyone who
has a copy; the removal is a pull request on public master like any other.

## 2. Prior art: legacy and literature

Written before the design.

### 2.1 Practice: how Python and C++ projects ship plugins

- **Entry points.** The Python Packaging User Guide, *Creating and discovering
  plugins* (packaging.python.org/en/latest/guides/creating-and-discovering-plugins/),
  names package metadata as the way a package "announces that it contains a
  specific kind of plugin", discovered with
  `importlib.metadata.entry_points(group=...)`. Used here unchanged: the group
  is `rasputin.commands`.
- **CLI commands as entry points.** Flask loads commands from the
  `flask.commands` entry-point group, declared in the plugin's
  `pyproject.toml` as `[project.entry-points."flask.commands"]`
  (flask.palletsprojects.com/en/stable/cli/, *Plugins*). The same shape here,
  with a Typer app in place of a Click command. **What differs:** Flask lets a
  later registration replace an earlier one; here a plugin can never replace a
  core command, and a broken plugin cannot stop the core CLI (section 4.1).
- **Typer.** `Typer.add_typer` without a name merges the sub-app's commands
  into the parent, later entries replacing earlier ones without a word, and
  refuses a sub-app that has a callback (read in Typer 0.27.2's
  `typer/main.py`, `get_group_from_info`). So the loader checks name clashes
  itself before merging, and never merges `cli.app` (which has a callback) as
  a nameless sub-app.
- **Headers shipped in the wheel.** NumPy's `numpy.get_include()` "Return[s]
  the directory that contains the NumPy *.h header files" for extension
  modules that compile against it (numpy.org, `numpy.get_include`); pybind11
  ships a CMake config in its wheel the same way. scikit-build-core's
  *Search paths* page calls `<PackageName>_ROOT`, filled from a dependency's
  `cmake.root` entry point, "the recommended interface to be used for
  importing dependent packages using `find_package`"
  (scikit-build-core.readthedocs.io/en/stable/configuration/search_paths.html).
  A downstream packager's changelog dates `cmake.root` to scikit-build-core
  0.11.0; not checked upstream, so pinned in section 10.
- **Two pybind11 modules.** pybind11 makes a class binding global across
  modules by default, shareable only with modules "built with sufficiently
  similar versions of pybind11 and of the C++ compiler and C++ standard
  library"; `py::module_local()` keeps a binding to one module
  (pybind11.readthedocs.io/en/stable/advanced/classes.html, *Module-local
  class bindings*). `module_local` alone does not keep two modules apart:
  pybind11 identifies a C++ type by its `typeid` name, and a module-local
  caster that finds no local match falls back to the global registration of
  the same name (section 4.3). The private module is built later, by another
  build, so this design shares **no** bound type between the two modules: only
  NumPy arrays and Python scalars cross, and the raster view each module binds
  is a distinct C++ type.

No novelty is claimed. This is packaging practice, applied.

### 2.2 Legacy

Nothing to carry across. The grep and the file list it returned:

```
$ git grep -l -iE 'entry_points|entry-points|plugin' legacy-archive -- legacy
legacy-archive:legacy/rasputin/reader.py
legacy-archive:legacy/tests/test_read_raster_file.py
```

Both are PIL's `TiffImagePlugin`, not a plugin mechanism. A second grep for
`snow|radiation` returned only land-cover class names
(`legacy/rasputin/globcov_repository.py`, `legacy/rasputin/gml_repository.py`).

## 3. The boundary

**Rule: private depends on public, never the reverse.** Nothing under the
public tree imports, includes, links or names the private package, apart from
the absence test of section 9.1 (T10), which names what must not be there.

### 3.1 What moves to the private package

Measured on master `a058bb6b`. The private increments stacked on it
(section 7) add files only under these paths, or to the snow command.

**This record names no private increment's content**: it is committed to the
public repository. What the private branches add is checked against this
table in the private repository, at Q1.

| Kind | Path on public master | In the private package |
|---|---|---|
| C++ headers | `include/terrain/snow/`, `include/terrain/radiation/` | same paths |
| Bindings | `bindings/snow.cpp`, `bindings/radiation.cpp` | same paths, plus `bindings/module.cpp` for `_snow` |
| Python | `src_python/tin_engine/snow/`, `src_python/tin_engine/radiation/` | `src/rasputin_snow/snow/`, `src/rasputin_snow/radiation/` |
| Python I/O used only by snow | `src_python/tin_engine/io/ensight.py` (its header line names `rasputin snow`) | `src/rasputin_snow/ensight.py` |
| CLI | the `synthetic-terrain` and `snow` commands and their helpers at the end of `cli.py` (`_snow_folder`, `_SNOW_STATIC`, `_Sun`, `_Shadows`, `_shadow_dem`, `_contains`, `_per_triangle`, and any the private increments add) | `src/rasputin_snow/cli.py`, a `typer.Typer` |
| Data | package data the private increments add under `src_python/tin_engine/data/` for snow | `src/rasputin_snow/data/` |
| Stubs | the snow and radiation entries of `_core.pyi` | `src/rasputin_snow/_snow.pyi` |
| C++ tests | `tests/cpp/unit/test_hbv_snow.cpp`, `test_hbv_snow_linear.cpp`, `test_radiation_{sun,receivers,horizon,evaluate}.cpp` and their six `add_terrain_test` lines | `tests/cpp/`, own `CMakeLists.txt` |
| Python tests | `test_cli_snow`, `test_core_snow`, `test_core_horizon`, `test_radiation_{evaluate,frame,search,sun}`, `test_snow_{forcing,geometry,job,params,run,terrain}`, `test_io_ensight`, plus the private increments' suites; helpers `ensightread.py` moved, `cli_helpers.py`, `gil_probe.py`, `vtkread.py` copied | `tests/python/` |
| Test-only dependency | `pvlib` in `dev` | the private `dev` extra |
| Docs | `docs/increments/46-hbv-snow.md`, `46-probes/`, `terrain-radiation.md`, `docs/benchmarks/2026-10-14/46-*`; the private increments' records and probes | `docs/` |

What stays public, though only snow uses it today: `io/geotiff_write.py`
(a GeoTIFF encoder, the reader's counterpart and a GIS output Ola wants) and
`io/vtk_read.py` (the mesh reader; its only production caller is `snow`,
`cli.py:3165`, but reading a mesh back is the core's own format, and
`tests/python/test_io_vtk_read.py` keeps it covered after P2).

Catchment hydrology (flood fill, accumulation, upstream) stays public: it
delineates catchments for the mesh, so it is meshing, not hydrological
modelling. Future runoff modelling goes private.

The private headers keep their paths and the `terrain::snow` and
`terrain::radiation` namespaces, so the move is a rename and not an edit. The
public core therefore never reuses `include/terrain/snow/` or
`include/terrain/radiation/`; T10 enforces it.

### 3.2 The three places the core points the wrong way today

1. `bindings/core.cpp` declares and calls `bind_snow` and `bind_radiation`
   inside `_core` (lines 54-55, 1576-1577), and `CMakeLists.txt` compiles
   both files into `_core` (line 89). After the split the private package
   builds its own extension module, `rasputin_snow._snow`.
2. `bindings/radiation.cpp` takes `BoundRasterView`, a type `core.cpp` binds,
   as `horizon_search`'s `near` and `far`. Resolved in section 4.3.
3. `cli.py` imports snow and radiation at module level (line 135, the EnSight
   writer, and lines 176-206).
   Replaced by the plugin hook (section 4.1).

## 4. What the public core exposes

Three narrow interfaces, each with a public test that fails when it changes
(section 9.1). Everything else in `tin_engine` stays free to change; a private
build that reaches past these is the private side's defect.

### 4.1 The plugin hook (Python)

New module `src_python/tin_engine/plugins.py`.

- **Group:** `rasputin.commands`. An entry point's value is a `typer.Typer`;
  each of its commands becomes a top-level `rasputin` command. The private
  package declares
  `[project.entry-points."rasputin.commands"] snow = "rasputin_snow.cli:app"`.
- **`discover() -> tuple[Plugin, ...]`:** one frozen record per entry point in
  the group: distribution name and version, entry-point name, and either the
  loaded `typer.Typer` or the reason it could not be loaded. Loading a broken
  plugin never raises.
- **`application(plugins: Sequence[Plugin] | None = None) -> typer.Typer`:**
  a **new** app holding the core's commands and then each plugin's; `None`
  means `discover()`. `cli.app` itself is never changed, so a test that wants
  the core alone uses `cli.app`, and a test of the hook passes its own plugins.
  Pure: the same plugins give the same app.
- **Refusals, each one line on stderr naming the distribution, and the core
  runs on:**
  - a plugin that fails to import, or whose value is not a `typer.Typer`:
    skipped, with the error;
  - a plugin command with the name of a core command: the core command is
    kept, the plugin's skipped;
  - two plugins with the same command name: the first in sorted
    distribution-name order is kept.
- **Who calls it:** the console script (`[project.scripts] rasputin =
  "tin_engine.plugins:main"`, where `main()` runs `application()`),
  `__main__.py` (the MCP server's child process), and the `mcp` command,
  which passes `application()` to `build_server`, so a plugin's commands are
  MCP tools too (question 5).
- **`rasputin version`** adds one line per loaded plugin, distribution name
  and version, so a run record can name what produced it.
- **Startup cost: discovery runs before every command, eagerly.** Measured
  2026-10-15 on an Apple M1 Max, Python 3.14.7, a venv of 69 installed
  distributions and no plugin, six fresh processes each:
  `entry_points(group="rasputin.commands")` takes 2.4 to 2.6 ms once
  `tin_engine.cli` is imported. Importing `importlib.metadata` itself costs
  about 18 ms in a bare interpreter, but `tin_engine.cli` already imports it
  (`installed_version`), so the core pays that today; a bare
  `import importlib.metadata` plus the scan is the roughly 22 ms the review
  saw. For scale: importing `tin_engine.cli` takes about 200 ms and
  `python -m tin_engine --help` about 310 ms wall time. The scan grows with the
  number of installed distributions (checked at 69 only). **The plugin's own
  import cost is unmeasured** until the private package exists (Q1); it is
  paid by every command, `mesh` included, because a Typer app must be loaded
  for its command names to appear. Lazy discovery is not needed for the scan:
  `rasputin mcp` pays it once per server start and a `mesh` run lasts seconds
  or more. Q1 measures `rasputin mesh --help` with and without the package,
  and if the package adds more than 50 ms (a sixth of today's `--help`), the
  private `cli.py` moves its heavy imports into the command bodies; that is a
  private change, and the public hook stays eager.

### 4.2 The Python API a plugin may import (`tin_engine.plugin_api`)

New module `src_python/tin_engine/plugin_api.py`: re-exports only, with an
`__all__`. The private package imports from `tin_engine` through this module
and nothing else (private test QT2). The list is what the private code at
tip of the private increments' branch (`<private-tip>`) imports or calls,
derived with:

```
for f in $(git ls-tree -r --name-only <private-tip> -- \
    src_python/tin_engine/snow src_python/tin_engine/radiation); do
  git show <private-tip>:$f | grep -E '^ *from tin_engine'; done
```

plus the names the moved CLI section calls.

| From | Names |
|---|---|
| `crs` | `frame_factors`, `parse_crs`, `reprojector`, `same_crs`, `site_of` |
| `dem_input` | `DemRequest`, `open_dem`, `repository_for` |
| `io.models` | `Bounds`, `DemTile`, `RasterMeta` |
| `io.vtk_read` | `read_mesh` |
| `io.geotiff_write` | `encode_geotiff` |
| `raster` | `core_args` (new, section 4.3) |
| `resources` | `Resources` |
| `run_record` | `plural` |
| `cli_kit` (new) | `destination`, `cli_resources`, `words`, `THREADS` |

`cli_kit.py` is today's `_destination`, `_resources`, `_words` and `THREADS`
moved out of `cli.py` unchanged; `cli.py` imports them back under their old
names (`from tin_engine.cli_kit import destination as _destination`, ...), so
no call site in `cli.py` changes and a private branch that edits the snow
command merges without conflict.

The private tests also use `features.DEFAULT_VOCABULARY`,
`io.vtk_legacy.write_vtk` and `io.geotiff.decode_dem` to build fixtures. They
are not promised; a change there breaks only private tests, which the nightly
run against public master (section 6) reports.

### 4.3 The raster across two extension modules

The one real design question in the split. The radiation kernel reads a
float32 DEM through `terrain::raster::RasterView<float>`; today Python builds
that view with `_core.raster_view` and hands the bound object to
`_core.horizon_search`. After the split `horizon_search` lives in `_snow`,
built separately, so `_core`'s bound type cannot be its argument (section
2.1).

- **Python: `tin_engine.raster.core_args(tile) -> RasterArgs`.** A frozen
  record of exactly what crosses into C++: the read-only array, `x_min`,
  `y_max`, `delta_x`, `delta_y`, `nodata`, after today's checks (a writeable
  array refused, a geographic DEM refused). `to_core(tile)` becomes
  `raster_view(**core_args(tile))`. Validation stays in one place.
- **C++: `bindings/numpy_raster.hpp`**, public, pybind11-dependent,
  header-only: `template <class Tag> struct TaggedRasterView`, today's
  `BoundRasterView` (from `bindings/bound_raster_view.hpp`, which it replaces)
  with a tag parameter that is never used in the body; today's `try_view` from
  `core.cpp` (lines 287-308) and the body of `_core.raster_view` (lines
  833-843), each made a template over `Tag` and otherwise unchanged. `_core`
  declares `struct CoreRasterTag;` in `core.cpp` and
  `using BoundRasterView = TaggedRasterView<CoreRasterTag>;`, so no call site
  in `core.cpp` changes. It is installed beside the
  headers (section 4.4) and lives under `bindings/`, never under
  `include/terrain/`, so the core headers stay free of pybind11 (the
  `modern-cxx` skill's binding isolation).
- **The private module** declares its own tag in its own source
  (`struct SnowRasterTag;` in the private `bindings/module.cpp`, never in a
  public header), binds `TaggedRasterView<SnowRasterTag>` as `RasterView`
  with `py::module_local()`, and its own `raster_view` function from that
  header. The private Python builds the view as
  `_snow.raster_view(**dataclasses.asdict(core_args(tile)))`.
- **Why the tag, read in pybind11 3.1.0** (the version in all 35 worktree venvs; `pyproject.toml` asks for `pybind11>=2.12`). pybind11 keys
  its type registry on `std::type_index`, and under libc++ compares by the
  `typeid` name (`detail/internals.h`, `same_type`, `type_hash` and
  `type_equal_to`, lines 182-206; under libstdc++ `type_info` equality is name-based too). A
  module-local caster that finds no local match retries with the global
  registration of the same C++ type (`detail/type_caster_base.h`, lines
  1294-1299), and accepts another module's module-local instance of the same
  type (`try_load_foreign_module_local`, lines 1198-1218, `same_type`). So
  with one struct in both modules, `module_local` or not, a `_core.RasterView`
  would load as `_snow`'s argument whenever the two share pybind11 internals.
  With distinct tags the two `typeid` names differ, both lookups miss, and the
  argument is refused with a `TypeError`. When the two modules do not share
  internals (another pybind11 ABI), neither registry sees the other's type and
  the result is the same `TypeError`. No C++ type is registered twice, so the
  `generic_type: type "RasterView" is already registered!` failure
  (`pybind11.h`, line 1823) cannot occur either. `module_local` is kept so
  `_snow`'s class never enters the global registry, but the tag, not
  `module_local`, is what refuses.
- **Considered and not taken:** `_snow.horizon_search` taking the arrays and
  scalars of both rasters directly, with no bound view type. It registers
  nothing, but gives the call twelve raster arguments in place of two, and
  moves the `nodata` and dtype checks into each kernel binding instead of the
  one `raster_view` function.
- Only NumPy arrays and Python scalars cross between `_core` and `_snow`. A
  `_core.RasterView` passed to `_snow.horizon_search` is a `TypeError`
  (QT3 checks it).

### 4.4 The C++ surface a plugin may include

The public wheel installs, in P1:

- `include/terrain/**` to `tin_engine/include/terrain/`, and
  `bindings/numpy_raster.hpp` to `tin_engine/include/rasputin_bindings/`;
- `tin_engine/share/cmake/rasputin/rasputinConfig.cmake`, hand-written,
  defining the imported INTERFACE target `rasputin::terrain_headers`: the
  include directory relative to the config file, `cxx_std_20`, and the
  hardening definitions the core was built with (`RASPUTIN_HARDENING`), so a
  plugin compiles the shared headers the way `_core` does;
- the entry point `[project.entry-points."cmake.root"] rasputin =
  "tin_engine"`, so a scikit-build-core build with the core in its build
  environment finds it with `find_package(rasputin CONFIG REQUIRED)`.

Promised headers (what the private headers, bindings and tests include,
followed transitively): `terrain/core/point.hpp`,
`terrain/raster/{geometry,raster,view,sample}.hpp`,
`terrain/parallel_util/chunks.hpp`, `rasputin_bindings/numpy_raster.hpp`. The
rest are installed but not promised. No compiled core library is shipped:
the private module needs none (the promised headers are header-only).

## 5. The private package

### 5.1 Layout, on `rasputin-private` master after Q1

```
pyproject.toml          rasputin-snow; scikit-build-core>=0.11
CMakeLists.txt          find_package(rasputin), pybind11_add_module(_snow ...)
include/terrain/{snow,radiation}/
bindings/{module,snow,radiation}.cpp
src/rasputin_snow/      cli.py, ensight.py, raster.py, snow/, radiation/,
                        data/, _snow.pyi, py.typed
tests/cpp/              CMakeLists.txt (Catch2 by FetchContent, as public)
tests/python/
docs/increments/        46, terrain-radiation and the private increments
docs/benchmarks/        the 46 acceptance runs
ROADMAP.md, CLAUDE.md, .claude/, tools/   (section 5.3)
```

### 5.2 Build

- **The core is pinned once, by full commit**, as both a build requirement
  (for the headers and the CMake config) and a runtime dependency (for the
  Python): `rasputin @ git+https://github.com/expertanalytics/rasputin@<sha>`
  in `[build-system].requires` and in `dependencies`. The two must be the
  same string (QT5). A direct git reference also sidesteps PyPI, where the
  name `rasputin` belongs to someone else; it means the private package is
  never uploaded to PyPI, which it is not meant to be.
- The headers come from the installed core, never from a second checkout, so
  the headers and the Python the package runs against are the same commit.
- Bumping the pin is a private PR whose CI is the compatibility check.
- `CMakeLists.txt`: `find_package(rasputin CONFIG REQUIRED)`,
  `find_package(pybind11 CONFIG REQUIRED)`,
  `pybind11_add_module(_snow bindings/module.cpp bindings/snow.cpp
  bindings/radiation.cpp)`, linked to `rasputin::terrain_headers` and
  `Threads::Threads`, `-Wall -Wextra -Wpedantic -Werror`, installed to
  `rasputin_snow`. The private C++ tests link the same target.
- Python imports: `tin_engine.snow.*` becomes `rasputin_snow.snow.*`, and so
  on; `from tin_engine import _core` in snow and radiation code becomes
  `from rasputin_snow import _snow`; core imports go through
  `tin_engine.plugin_api`. `importlib.resources.files("tin_engine")` in
  package data becomes `files("rasputin_snow")`.
- Licence: proprietary, Expert Analytics, all rights reserved (the licensing
  decision's default), stated in the private `LICENSE` and `pyproject.toml`.
  The public `LICENSE` and `NOTICE.md` are untouched.

### 5.3 Working in the private repository

Agents run in the repository a session starts in (the spawn guard refuses a
spawn outside it), so snow work runs in a clone of `rasputin-private` with its
own sessions. Q1 therefore keeps the agent harness in the private repository
(`CLAUDE.md`, `.claude/`, the governance tools and their tests under
`tools/` and `tests/python/`), with the private `CLAUDE.md` §4 changed for the
private build and merge (below). A harness change is made in public first and
brought across by a private PR that checks out the harness paths from a named
public commit (`git checkout <public-sha> -- .claude tools ...`, the public
repository added as a remote in the private clone). `@orchestrator` owns
keeping the two in step (question 3).

Increment numbers stay one sequence across both repositories, so a number is
never reused; the public ROADMAP has no rows for private increments.

## 6. Private CI

On `rasputin-private`, GitHub Actions, Team plan: no merge queue.

- **Gating, every PR, Linux:** install the package (which builds the pinned
  core and `_snow`); `ruff check`, `ruff format --check`, `mypy` (strict over
  `src/rasputin_snow`); `pytest`; configure `tests/cpp` against the installed
  core and run `ctest`, reading `cmake --build`'s exit status first;
  `tools/check_prohibited_deps.py`; `tools/count_loc.py` against the base.
  Branch protection requires this one job.
- **macOS:** on demand (`workflow_dispatch`) and weekly, not on every PR
  (macOS minutes cost ten times Linux minutes on a private repository;
  question 4).
- **Nightly, not gating:** the same job with the pin replaced by public
  `master`'s head, so a public change that breaks the package shows within a
  day rather than at the next pin bump.
- **Merging:** `gh pr merge <pr> --merge` once the gating job is green; a
  merge commit, never a squash, as on public. Ola's yes per merge, as on
  public.

## 7. Where the private increments land

Two approved private increments are stacked on public master, as one chain
of branches; the tip is called `<private-tip>` here. Their branch names and
contents stay out of this public record. Checked with
`git diff --stat master...<private-tip>`: every file they touch is in a
moving path of section 3.1, or is a shared file whose snow part moves
(`cli.py`'s snow command, `_core.pyi`, `test_layering.py`,
`project_structure.md`, `ROADMAP.md`). Part of the chain is on no remote yet
(`git ls-remote private`).

**They land on private master first, as they are**, before Q1 (Q0 below):
private master is today a copy of public master (`private/master` =
`a058bb6b`), so each merges as it was reviewed, and Q1 then moves their files
with the rest, with their history. They never go to the public remote.

## 8. The PR sequence

Each PR is under the 700-line rule of `CLAUDE.md` §2, counted by
`tools/count_loc.py`. Estimates; the reviewer counts.

| # | Repo | What | Net counted lines |
|---|---|---|---|
| Q0 | private | the private increments (section 7) merged into private master, one PR each, in chain order | as each was reviewed |
| P1 | public | the seam: `plugins.py`, `plugin_api.py`, `cli_kit.py` (moved), `raster.core_args`, `bindings/numpy_raster.hpp` (moved), header and CMake-config install, `cmake.root` entry point, console script and `mcp` through `application()`, `rasputin version` lists plugins | about +190 |
| P2 | public | the removal: everything in section 3.1's left column, `bind_snow`/`bind_radiation` out of `core.cpp` and `CMakeLists.txt`, the six C++ tests out of `tests/cpp/CMakeLists.txt` and CI's sanitizer list, `pvlib` out of `dev`, the snow and radiation paths out of `[tool.check_code_history]`, layering rows, `project_structure.md` entries, ROADMAP row 46 reduced to a pointer | about −2,600 |
| Q1 | private | the package: core copy deleted from the private tree, snow files moved into the layout of 5.1, `_snow` module, CMake and pyproject, pin = P2's merge commit, CI workflow, private `CLAUDE.md` §4 | far below zero; about +150 added |

**Order, and why neither repository breaks:**

1. Q0 and P1 are independent and may run side by side.
2. P2 after P1, so the hook exists before anything relies on it. From P2 on,
   public master has no snow; snow still runs from private master, which is
   still a full copy.
3. Q1 after P2 merges, pinned to **P2's merge commit**, not P1's. Pinned to
   P1, the core would still carry its own `snow` command, the clash rule
   would keep it, and the private CLI tests would exercise the core's copy,
   not the package's.
4. Private master never merges public master after P2 (that would delete the
   snow code before Q1 moves it).
5. Before Q1, a tag on private master keeps the last full copy
   (`full-copy-<date>`); pushing it needs Ola's yes.

`@perf`'s acceptance run does not apply: no refine or mesh code changes.

## 9. Tests

`@tester` writes them red, per PR. Fixtures are tiny and built in the test.

### 9.1 Public (P1 unless marked)

The plugin fixture: a directory under `tmp_path` holding one module
(`fakeplug.py`, one Typer command `hello` that prints `hello`, 6 lines) and
one `fakeplug-0.1.dist-info/` (`METADATA` with name and version,
`entry_points.txt` with one `rasputin.commands` entry), put on `sys.path` (or
`PYTHONPATH` for a subprocess), so real `importlib.metadata` discovery runs,
no monkeypatching. Variants change one line each.

- **T1** `application()` with the fixture lists `hello`, and
  `python -m tin_engine --help` in a subprocess lists it; `hello` runs.
- **T2** a plugin command named `mesh`: `mesh` is still the core's (its help
  text unchanged), one stderr line names `fakeplug`, exit 0 for `--help`.
- **T3** a plugin module that raises `ImportError`, and one whose value is
  not a `typer.Typer`: the core's commands all work; one stderr line each,
  naming the distribution.
- **T4** two plugins with the same command: the first in distribution-name
  order wins, one stderr line.
- **T5** `rasputin version` with the fixture prints `fakeplug 0.1`.
- **T6** `build_server(application(...))` derives a tool `hello`.
- **T7** `plugin_api.__all__` equals section 4.2's list, and each name
  resolves to the object in its home module (identity, not equality).
- **T8** `core_args` refuses a writeable array and a geographic tile with
  today's messages; `to_core(t)` gives the same view as
  `raster_view(**asdict(core_args(t)))` (same sampled values at 9 points of a
  3 × 3 float32 tile).
- **T9** `encode_geotiff` round trip: a 4 × 5 float32 grid, EPSG:25833,
  10 m, read back by `io.geotiff`: same array, origin, step and CRS. (Its only
  callers' tests leave in P2; this keeps it covered.)
- **T10 (P2)** the absence test: none of `include/terrain/snow`,
  `include/terrain/radiation`, `bindings/snow.cpp`, `bindings/radiation.cpp`,
  `src_python/tin_engine/{snow,radiation}`, `io/ensight.py` exists;
  `tin_engine.snow` and `tin_engine.radiation` do not import; `_core` has none
  of `SnowState`, `SnowDistribution`, `SunCorrection`, `horizon_search`,
  `sun_positions`, `solar_position`, `evaluate_instant`, `delta_t_seconds`;
  `pvlib` is not in `pyproject.toml`. It also guards against a private branch
  merged into public by mistake.
- **CI step (P1), consumer build:** in the Python job, after the wheel is
  installed, configure and build a 2-file CMake project under
  `tests/consumer/` with `find_package(rasputin CONFIG REQUIRED)`: one
  translation unit includes every promised header, builds a
  `RasterView<float>` over a 2 × 2 array and samples it; a second includes
  `rasputin_bindings/numpy_raster.hpp` with pybind11. The step fails if a
  promised header is not installed or the config is wrong; checked by
  deleting the header install rule once, locally, before the commit.
- **What proves the core still works with no snow code:** after P2, the full
  public suite, `ctest` and every gate green in CI, with the counts of
  collected Python tests and registered C++ tests recorded in P2's review
  next to master's, the difference equal to the moved tests.

### 9.2 Private (Q1)

- **QT1** every moved test passes, changed only in imports (`tin_engine.snow`
  to `rasputin_snow.snow`, `_core` to `_snow` for snow and radiation names,
  helpers' imports). The reviewer checks with `git diff -M` on `tests/` that
  no other line changed, and that the count of collected tests equals the
  count on private master before Q1.
- **QT2** an import scan over `src/rasputin_snow`: the only `tin_engine`
  module imported is `tin_engine.plugin_api`.
- **QT3** one process imports `tin_engine._core` and `rasputin_snow._snow`,
  builds a raster view in each, and runs `_snow.horizon_search` on its own
  view; a `_core.RasterView` passed to it is a `TypeError`. The check is
  planted once, locally, before Q1's red commit: `_snow` built with
  `CoreRasterTag` in place of `SnowRasterTag` must make QT3 fail (the
  `_core` view is then accepted), and the result is recorded in Q1's review.
- **QT4** the installed console script, in a subprocess: `rasputin --help`
  lists `snow` and `synthetic-terrain`; `rasputin version` lists
  `rasputin-snow`.
- **QT5** the core pin in `[build-system].requires` equals the one in
  `dependencies`.
- **QT6, the split changes no result:** a probe, run once and recorded in
  Q1's review: `synthetic-terrain`, `mesh` with `--tolerance 5`, then
  `snow --days 10` with the 1 April default start, once with the last full
  copy (the tag of section 8) and once with the package; every `.var` and
  `snow.geo` file byte-identical. Its figures (triangles, frames, bytes) are
  measured when it runs; none is assumed here.

No mutation round: the loader's rules are small and each has its own test.

## 10. Pinned or assumed beyond Ola's words

1. The private package's working name `rasputin-snow`, import name
   `rasputin_snow` (question 1).
2. `rasputin-private` is turned into the package's repository; the full copy
   is kept as a tag (question 2).
3. `io/geotiff_write.py` and `io/vtk_read.py` stay public; `io/ensight.py`
   moves.
4. The private C++ keeps `include/terrain/{snow,radiation}/` and the
   `terrain::` namespaces; public never reuses those paths.
5. `scikit-build-core>=0.11` honours the `cmake.root` entry point (from a
   packager's changelog, not upstream). If Q1's build shows otherwise, the
   private `CMakeLists.txt` sets `rasputin_ROOT` from
   `importlib.util.find_spec("tin_engine")` instead; the public side is the
   same either way.
6. A broken plugin is skipped with a warning rather than stopping the CLI.
7. Increment numbers are one sequence across both repositories.
8. The design records (46, terrain-radiation, the 46 benchmarks) leave public
   master in P2; the retrospectives that mention them stay as written.
9. This public record names no private increment's content; section 1 quotes
   the question Ola answered verbatim, which names wind modelling as a
   subject.

## 11. Questions for Ola

1. **The private package's name?** Default: `rasputin-snow`, imported as
   `rasputin_snow`, even though it will also hold radiation and later
   runoff. A broader name such as `rasputin-hydro` is cheap now and costly
   later.
2. **Turn the existing private repository `rasputin-private` into the
   package's repository**, after tagging its last full copy? Default: yes.
3. **Keep a copy of the agent harness (`CLAUDE.md`, `.claude/`, the
   governance tools) in the private repository**, synced from public by a
   private PR when it changes? Default: yes; without it no agent can work
   there.
4. **Private CI on Linux for every PR, macOS weekly and on demand?** Default:
   yes; macOS on every PR costs about ten times as much per minute.
5. **Should the snow commands also appear as tools in `rasputin mcp`** when
   the private package is installed? Default: yes.

Ruled by Ola, 2026-10-15: "Yes to all." All five defaults hold.

## 12. Ideas file

The ideas file is private; the entries this design takes up, overlaps and
leaves out are named in `@architect`'s handback, not here. One rule for them:
a private idea that needs more from the core than section 4.2 lists adds to
`plugin_api` in a public PR of its own, which says what it exposes and not
what uses it.

## Review

**Design review, round 1, 2026-10-15.** Range `master..ac9fef5a`, 0 counted lines, @reviewer: CHANGES REQUESTED. (1) Section 4.3 and QT3: pybind11 3.1.0 finds types by name (pybind11 header detail/internals.h, line 202), and `_snow`'s module-local caster falls back to `_core`'s global registration (pybind11 header detail/type_caster_base.h, lines 1294-1299), so a `_core.RasterView` is accepted by `_snow.horizon_search`, not refused; bind a differently named C++ type on the private side (e.g. a tag template parameter). (2) No startup cost stated for `discover()` before every command (entry-point scan about 22 ms warm with no plugin); state it, plugin import cost unmeasured until Q1. Wordings: `io/vtk_read.py`'s only production caller is `snow` (`cli.py:3165`); add `_SNOW_STATIC` (`cli.py:3069`) to the moving helpers; section 3.2 cites line 135 too; the status line and ROADMAP row 51 say the five questions were ruled at their defaults on 2026-10-15.
