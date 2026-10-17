# Writing a rasputin plugin

A plugin is a Python package, installed beside rasputin, that adds commands to
`rasputin`. It may also ship a C++ extension module built against rasputin's
installed headers. This page says what a plugin declares, what it may import,
and what it may include.

## Add a command

Declare an entry point in the group `rasputin.commands` whose value is a
`typer.Typer`. Each command of that app becomes a top-level `rasputin`
command, and a tool of the MCP server that `rasputin mcp` runs.

A package `hello_plugin`, with one module:

```python
# hello_plugin/__init__.py
import typer

app = typer.Typer()


@app.command()
def hello() -> None:
    """Print hello."""
    typer.echo("hello")
```

and in its `pyproject.toml`:

```toml
[project.entry-points."rasputin.commands"]
hello = "hello_plugin:app"
```

Once the package is installed in rasputin's environment, `rasputin hello`
prints `hello`, and `rasputin version` lists the plugin's distribution and
version.

## Refusals

A plugin never replaces a core command, and a broken plugin never stops the
core. Each refusal is one line on stderr, starting `rasputin: plugin <name>`,
and the core runs on:

- a plugin that fails to import, or whose entry point's value is not a
  `typer.Typer`, is skipped;
- a command with the name of a core command is skipped; the core keeps its own;
- of two plugins with a command of the same name, the first in
  distribution-name order keeps it, and the other's is skipped.

`rasputin version` lists only the plugins that loaded.

## Output paths

`rasputin mcp` lets a tool call write only inside its output root
(`--out-root`, default `~/rasputin-output`). It tells an output path from an
input path by Click's own marker: a path option or argument declared
`writable=True`. Mark every option your command writes to:

```python
out: Annotated[Path, typer.Option("--out", writable=True, help="Where to write.")]
```

A value of such an option that resolves outside the root, through `..` or a
symbolic link too, is refused before your command runs. An option left
unmarked is not fenced. Typer ignores `writable` on a `str` option, so
declare it as `Path`; add `allow_dash=True` where `-` means stdout. Pass the
path to `destination` (below) before writing: it refuses a symbolic link as
the file itself and makes the missing folders.

## Ship datasets

A plugin may add datasets to the catalogue that `rasputin sources` lists
and `--dem` reads by id. Declare an entry point in the group
`rasputin.sources` whose value names a package, not an attribute; the
package holds a `sources.toml` in the format of your own sources file
(`[source.<id>]` tables; see the README, *Connect an agent*):

```toml
[project.entry-points."rasputin.sources"]
mydata = "my_plugin.data"
```

A relative `path` in it is taken from the package's directory. A plugin
needs no command to ship data. The catalogue is the package's sources, then
the user's file, then each plugin's in distribution-name order, and a plugin
never replaces an entry: its whole file is skipped, with one line on stderr
starting `rasputin: plugin <name>`, when it cannot be read or validated or
when one of its ids is already taken. `rasputin sources` lists the skipped
files under `plugin_problems`, and each source's `origin` says where it
came from (`plugin <distribution>`).

## Where it appears in `--help`

A command that sets `rich_help_panel` is shown under that panel. Any other
plugin command is shown under *Plugins*, after the core's panels. A command
with `hidden=True` is left out of `--help` and is not offered as an MCP tool.

## What you may import

Import the core through `tin_engine.plugin_api` and nothing else. A name in
it is promised: it is removed or changed only in a release that says so. Any
other module of `tin_engine` may change in any commit. The list is what
plugins use today, and grows on request.

| Name | What it is |
|---|---|
| `THREADS` | The `--threads` option, as a Typer annotation |
| `Resources` | What a run may use of the machine: threads and a memory budget |
| `cli_resources` | `--threads` and `--memory-budget` resolved against the environment, as `Resources` |
| `destination` | `destination(out, hint="--out")`: the resolved output path, refusing a symbolic link as the file itself and making the missing folders; a refusal names the option `hint`. It no longer takes `out_parent` or `name`: the output root of `rasputin mcp` replaced them |
| `words` | A refusal's own words from a `ValueError`, without Pydantic's wrapper |
| `plural` | `"1 tile"`, `"2 tiles"` |
| `DemRequest` | Which DEM to read: a directory, files or a cached source, and an optional box or domain |
| `open_dem` | Lists, plans and decodes a `DemRequest` into one grid |
| `repository_for` | The tile repository over a directory, files or a cached source |
| `DemTile` | One decoded tile: its `RasterMeta` and a read-only 2-D array |
| `RasterMeta` | A tile's node grid: origin, spacing, size, CRS, nodata |
| `Bounds` | A box in a DEM's CRS |
| `core_args` | A tile's array and grid scalars, checked, for a `raster_view` function (see the raster tag rule below) |
| `accumulate` | `accumulate(tile)`: each node's drainage over a `DemTile`, by the Priority-Flood that fills its pits, on the node grid alone (spacing and CRS unused). A node drains to the neighbour that flooded it: its lowest filled neighbour, ties first in, first out; not D8's steepest descent. `.count` (uint32, rows × cols): 0 on NoData, else the nodes draining through the node, itself included. `.flow_to` (uint8): that neighbour, `3*(dr+1)+(dc+1)`, with `dr`, `dc` its row and column minus the node's; 255 for an outlet or NoData. Both read-only. 2^32 nodes or more is a `ValueError`. Releases the GIL. The result's type and other attributes are not promised |
| `encode_geotiff` | A 2-D grid as the bytes of a float32 GeoTIFF |
| `read_mesh` | The bytes of a binary `.vtk` that `rasputin mesh` wrote, as arrays |
| `parse_crs` | Any CRS text pyproj reads, or a `ValueError` naming it |
| `same_crs` | Whether two CRSs name the same points |
| `reprojector` | A vectorised map of `(x, y)` points from one CRS to another |
| `frame_factors` | Meridian convergence and point scale of a projected CRS |
| `site_of` | Where a point of a projected CRS lies on the globe |

Where a name has a docstring, it is the full contract.

## C++

The wheel installs rasputin's headers and a CMake config:

- `include/terrain/**`, under `tin_engine/include/terrain/`;
- `bindings/numpy_raster.hpp`, under `tin_engine/include/rasputin_bindings/`;
- `tin_engine/share/cmake/rasputin/rasputinConfig.cmake`, which defines the
  target `rasputin::terrain_headers`: the include directory, C++20, and the
  bounds-check setting the core was built with.

Promised headers: `terrain/core/point.hpp`,
`terrain/raster/{geometry,raster,view,sample}.hpp`,
`terrain/parallel_util/chunks.hpp` and `rasputin_bindings/numpy_raster.hpp`.
The others are installed but may change in any commit. All are header-only;
there is no compiled core library to link.

rasputin declares the entry point `[project.entry-points."cmake.root"]
rasputin = "tin_engine"`, so a scikit-build-core build with rasputin in its
build environment finds the config with no path given:

```cmake
find_package(rasputin CONFIG REQUIRED)
target_link_libraries(my_module PRIVATE rasputin::terrain_headers)
```

Elsewhere, pass `-Drasputin_ROOT=<the installed tin_engine directory>`.
`tests/consumer/CMakeLists.txt` in this repository is a complete example: an
executable and a pybind11 module built against the installed core alone.

## The raster tag rule

To read a DEM in C++, bind your own raster view from
`rasputin_bindings/numpy_raster.hpp`:

```cpp
#include <rasputin_bindings/numpy_raster.hpp>
#include <pybind11/stl.h>  // std::optional, for nodata

struct MyRasterTag;  // in your own source, never in a header you share
using MyRasterView = TaggedRasterView<MyRasterTag>;

PYBIND11_MODULE(_mine, m) {
    pybind11::class_<MyRasterView>(m, "RasterView", pybind11::module_local());
    m.def("raster_view", &rasputin_bindings::raster_view<MyRasterTag>,
          pybind11::arg("array"), pybind11::kw_only(), pybind11::arg("x_min"),
          pybind11::arg("y_max"), pybind11::arg("delta_x"), pybind11::arg("delta_y"),
          pybind11::arg("nodata") = pybind11::none());
}
```

and build it in Python from `core_args`:

```python
import dataclasses
from tin_engine.plugin_api import core_args

view = _mine.raster_view(**dataclasses.asdict(core_args(tile)))
```

Never accept rasputin's own `_core.RasterView`: only NumPy arrays and Python
scalars cross between two extension modules. With a tag of your own, your
view is a different C++ type from rasputin's, so pybind11 refuses one where
the other is expected with a `TypeError`, instead of accepting an object whose
layout another module's build decided. The full reasoning is in
`docs/increments/51-open-core.md`, *4.3 The raster across two extension
modules*.
