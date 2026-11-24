# Increment 51: the plugin seam

**Status: P1 merged** (#270, `3663448f`); **P2 merged** (#271, `9a7c312d`).

This is the public record of the seam the core offers a plugin: the plugin
hook, the Python API a plugin may import, the raster across two extension
modules, the installed C++ headers and the consumer build that tests them. A
private package, `rasputin-lsm` in the private repository
`expertanalytics/rasputin-private`, builds on it. Its own record is private.

## 1. Prior art

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
  core command, and a broken plugin cannot stop the core CLI (section 2.1).
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
- **Two pybind11 modules.** pybind11 makes a class binding global across
  modules by default; `py::module_local()` keeps a binding to one module
  (pybind11.readthedocs.io/en/stable/advanced/classes.html, *Module-local
  class bindings*). `module_local` alone does not keep two modules apart
  (section 2.3). So the seam shares **no** bound type between the core's
  module and a plugin's: only NumPy arrays and Python scalars cross.

No novelty is claimed. This is packaging practice, applied.

**Rule: a plugin depends on the core, never the reverse.** Nothing under the
public tree imports, includes, links or names a plugin's code.

## 2. What the public core exposes

Three narrow interfaces, each with a public test that fails when it changes
(section 3). Everything else in `tin_engine` stays free to change; a plugin
that reaches past these carries its own defect.

### 2.1 The plugin hook (Python)

Module `src_python/tin_engine/plugins.py`.

- **Group:** `rasputin.commands`. An entry point's value is a `typer.Typer`;
  each of its commands becomes a top-level `rasputin` command.
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
  MCP tools too.
- **`rasputin version`** adds one line per loaded plugin, distribution name
  and version, so a run record can name what produced it.
- **Startup cost: discovery runs before every command, eagerly.** Measured
  2026-10-15 on an Apple M1 Max, Python 3.14.7, a venv of 69 installed
  distributions and no plugin, six fresh processes each:
  `entry_points(group="rasputin.commands")` takes 2.4 to 2.6 ms once
  `tin_engine.cli` is imported, which already imports `importlib.metadata`.
  For scale: importing `tin_engine.cli` takes about 200 ms and
  `python -m tin_engine --help` about 310 ms wall time.

### 2.2 The Python API a plugin may import (`tin_engine.plugin_api`)

Module `src_python/tin_engine/plugin_api.py`: re-exports only, with an
`__all__`. A plugin imports from `tin_engine` through this module and nothing
else.

| From | Names |
|---|---|
| `crs` | `frame_factors`, `parse_crs`, `reprojector`, `same_crs`, `site_of` |
| `dem_input` | `DemRequest`, `open_dem`, `repository_for` |
| `io.models` | `Bounds`, `DemTile`, `RasterMeta` |
| `io.vtk_read` | `read_mesh` |
| `io.geotiff_write` | `encode_geotiff` |
| `raster` | `core_args` (section 2.3) |
| `resources` | `Resources` |
| `run_record` | `plural` |
| `cli_kit` | `destination`, `cli_resources`, `words`, `THREADS` |

`cli_kit.py` holds `_destination`, `_resources`, `_words` and `THREADS`
moved out of `cli.py` unchanged; `cli.py` imports them back under their old
names (`from tin_engine.cli_kit import destination as _destination`, ...).
Increment 64 adds `accumulate`.

A plugin that needs more from the core adds to `plugin_api` in a public PR
of its own, which says what it exposes and not what uses it.

### 2.3 The raster across two extension modules

A plugin's C++ kernel may read a float32 DEM through
`terrain::raster::RasterView<float>`. Python builds that view with
`_core.raster_view`, but a plugin's module is built separately, so `_core`'s
bound type cannot be its argument.

- **Python: `tin_engine.raster.core_args(tile) -> RasterArgs`.** A frozen
  record of exactly what crosses into C++: the read-only array, `x_min`,
  `y_max`, `delta_x`, `delta_y`, `nodata`, after the checks (a writeable
  array refused, a geographic DEM refused). `to_core(tile)` is
  `raster_view(**core_args(tile))`. Validation stays in one place.
- **C++: `bindings/numpy_raster.hpp`**, public, pybind11-dependent,
  header-only: `template <class Tag> struct TaggedRasterView`, with a tag
  parameter that is never used in the body, and `try_view` and the body of
  `raster_view`, each a template over `Tag`. `_core` declares
  `struct CoreRasterTag;` in `core.cpp` and
  `using BoundRasterView = TaggedRasterView<CoreRasterTag>;`. It is installed
  beside the headers (section 2.4) and lives under `bindings/`, never under
  `include/terrain/`, so the core headers stay free of pybind11 (the
  `modern-cxx` skill's binding isolation).
- **A plugin's module** declares its own tag in its own source, never in a
  public header, binds `TaggedRasterView<ItsTag>` as `RasterView` with
  `py::module_local()`, and its own `raster_view` function from that header.
  Its Python builds the view as
  `module.raster_view(**dataclasses.asdict(core_args(tile)))`.
- **Why the tag, read in pybind11 3.1.0.** pybind11 keys its type registry on
  `std::type_index`, and under libc++ compares by the `typeid` name
  (`detail/internals.h`, `same_type`, `type_hash` and `type_equal_to`, lines
  182-206; under libstdc++ `type_info` equality is name-based too). A
  module-local caster that finds no local match retries with the global
  registration of the same C++ type (`detail/type_caster_base.h`, lines
  1294-1299), and accepts another module's module-local instance of the same
  type (`try_load_foreign_module_local`, lines 1198-1218). So with one struct
  in both modules, `module_local` or not, a `_core.RasterView` would load as
  the plugin's argument whenever the two share pybind11 internals. With
  distinct tags the two `typeid` names differ, both lookups miss, and the
  argument is refused with a `TypeError`; pybind11's C++ conduit between
  modules that do not share internals (`detail/type_caster_base.h`, lines
  1317 and 1329-1352) looks up the same name and misses the same way. No C++
  type is registered twice, so `generic_type: type "RasterView" is already
  registered!` (`pybind11.h`, line 1823) cannot occur. `module_local` keeps
  the plugin's class out of the global registry, but the tag is what refuses.

### 2.4 The C++ surface a plugin may include

The public wheel installs:

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

Promised headers: `terrain/core/point.hpp`,
`terrain/raster/{geometry,raster,view,sample}.hpp`,
`terrain/parallel_util/chunks.hpp`, `rasputin_bindings/numpy_raster.hpp`. The
rest are installed but not promised. No compiled core library is shipped.

## 3. Tests

### 3.1 Public tests of the seam

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
- **T7** `plugin_api.__all__` equals section 2.2's list, and each name
  resolves to the object in its home module (identity, not equality).
- **T8** `core_args` refuses a writeable array and a geographic tile with
  the core's messages; `to_core(t)` gives the same view as
  `raster_view(**asdict(core_args(t)))` (same sampled values at 9 points of a
  3 × 3 float32 tile).
- **T9** `encode_geotiff` round trip: a 4 × 5 float32 grid, EPSG:25833,
  10 m, read back by `io.geotiff`: same array, origin, step and CRS.
- **CI step, consumer build:** in the Python job, after the wheel is
  installed, configure and build a 2-file CMake project under
  `tests/consumer/` with `find_package(rasputin CONFIG REQUIRED)`: one
  translation unit includes every promised header, builds a
  `RasterView<float>` over a 2 × 2 array and samples it; a second includes
  `rasputin_bindings/numpy_raster.hpp` with pybind11. The step fails if a
  promised header is not installed or the config is wrong.

`io/geotiff_write.py` and `io/vtk_read.py` are public: a GeoTIFF encoder,
the reader's counterpart, and the reader of the core's own mesh format;
`tests/python/test_io_vtk_read.py` covers the second.

## 4. As built

- **P1** (#270, `3663448f`): `plugins.py`, `plugin_api.py`, `cli_kit.py`,
  `raster.core_args`, `bindings/numpy_raster.hpp`, the header and
  CMake-config install, the `cmake.root` entry point, the console script and
  `mcp` through `application()`, `rasputin version` listing plugins. Code
  review: 157 net production lines, approved with wordings.
- **P2** (#271, `9a7c312d`): the code that moved to the private package left
  the public tree; −2,317 net production lines, approved with wordings.
