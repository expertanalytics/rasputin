# Rasputin

Rasputin turns a digital elevation model (DEM) into a triangulated irregular
network (TIN): a triangle mesh of the terrain that uses few triangles where the
ground is smooth and many where it is rough, and never strays from the DEM by
more than a height tolerance you choose.

![The Bygdin catchment meshed at 10 m tolerance, coloured by CORINE land cover](docs/benchmarks/2026-09-29/bygdin-landcover/bygdin_landcover_oblique.png)

*The Bygdin catchment (305 km²) at a 10 m tolerance, each triangle coloured by
its CORINE land-cover class. Elevation © Kartverket (CC BY 4.0); land cover
© European Union, Copernicus Land Monitoring Service 2018, EEA. How it was made:
`docs/benchmarks/2026-09-29/bygdin-landcover.md`.*

## What it does

- **DEM in, mesh out.** One GeoTIFF, several, or a directory of tiles
  (Kartverket's DTM10, for example), stitched into one grid.
- **DEMs from the network.** `rasputin fetch glo30 --domain X` copies the
  Copernicus GLO-30 tiles a mesh of that domain reads into a tile cache, so
  `rasputin mesh --dem glo30` then runs offline; `anadem-v1` (South America)
  likewise.
- **A domain.** Mesh only the inside of a polygon, such as a catchment
  (GeoJSON or WKT, in any CRS that pyproj reads).
- **Features as constraints.** Lakes, land-cover borders, roads and rivers
  (GeoJSON, GeoPackage or GML) become edges of the mesh, so no triangle
  crosses them.
- **Land-cover labels, cleaned up.** With a CORINE class map, every triangle
  gets the class of the polygon it lies in. First the borders are repaired
  (`--features-repair`), merged between neighbours of one class, snapped onto
  the domain outline where they nearly touch it, and simplified to
  `--features-tolerance` metres with each class keeping its area.
- **Refinement to a tolerance.** Triangles are split at the worst DEM point
  until every triangle is within `--tolerance` metres of the DEM, then kept
  Delaunay, with a minimum-angle start so there are no slivers along the
  boundary.
- **A tolerance that varies.** `--tolerance-near FILE N` with
  `--tolerance-ramp START END` holds the mesh to N metres along named lines
  (a railway, a river) and ramps to the general tolerance with distance;
  `--tolerance-slope N START END` does the same on steep ground.
- **Catchments.** `rasputin catchment` computes the catchment of a lake
  (a `--seed` point in it) or of a gauge on a mapped river (`--rivers`) from
  the DEM, and writes a GeoJSON polygon ready for `--domain`.
  `rasputin fetch-stations nve-hrd` fetches NVE's reference gauging stations
  with their catchments and rivers; `rasputin station-catchments` runs them
  all and compares.

## Quick example

Install first ([INSTALL.md](INSTALL.md)); the example needs the `codecs`
extra, because the committed DEM tile is LZW-compressed. The install is
bounds-checked by default; INSTALL.md's "Bounds checks" section says how to
build an unchecked copy for heavy runs. Every run uses every core unless
`--threads N` or the environment variable `RASPUTIN_THREADS` says otherwise
(0 means every core). From the repository root:

```sh
# A 50 km DTM10 tile, meshed inside a quarter disc of radius 30 km,
# to within 1 m of the DEM.
rasputin mesh \
    --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif \
    --domain docs/benchmarks/2026-09-26/quarter.geojson \
    --tolerance 1 \
    --out quarter.vtk

# The same, with CORINE land cover as constraints and labels,
# at a 5 m tolerance.
rasputin mesh \
    --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif \
    --domain docs/benchmarks/2026-09-26/quarter.geojson \
    --features tests/fixtures/corine/clc2018_7908_3.gpkg --features-map corine \
    --tolerance 5 \
    --out quarter_landcover.vtk

# Natural colours for the land-cover classes, for ParaView.
rasputin palette corine --out corine.json
```

The first command writes about 430,000 triangles in under two seconds on a
laptop (Apple M1 Max). Add `--stats -` to print sizes, mesh quality and
timings, and `--record run.json` to write the run's record, for reproducing it.

`rasputin --help` lists the commands, and `rasputin mesh --help` explains
every option.

## Output and viewing

- **`.vtk`** (legacy VTK, one file) is for [ParaView](https://www.paraview.org).
  It holds the triangles, the constraint edges with their feature bits and,
  with a class map, the cell array `land_cover_code`. To colour by land cover,
  import `corine.json` in ParaView's Colour Map Editor (Choose Preset, Import)
  and colour by `land_cover_code`.
- **`.ply`** is for QGIS and other mesh tools. The surface goes to `--out`,
  and the constraint edges, if you want them, to a second file named by
  `--out-edges`.

Both formats are binary by default; `--ascii` writes text you can read with
`head`. Every mesh file records what it was made from, and stderr gets a
one-line summary of the run.

## Connect an agent

An AI agent (Claude Desktop, Codex, or any other client of the Model
Context Protocol, MCP) can drive rasputin through a small server built in.
Install the `mcp` extra (the `dev` extra includes it):

```sh
pip install "rasputin[mcp]"
```

`rasputin mcp` then serves the commands over stdin and stdout. Point the
client at the `rasputin` of the environment you installed into, or at
`python -m tin_engine mcp` with that environment's `python`. Claude Desktop
takes this in `claude_desktop_config.json` (on macOS under
`~/Library/Application Support/Claude/`):

```json
{"mcpServers": {"rasputin": {"command": "/path/to/venv/bin/rasputin", "args": ["mcp"]}}}
```

Codex takes the same two fields in `~/.codex/config.toml`:

```toml
[mcp_servers.rasputin]
command = "/path/to/venv/bin/rasputin"
args = ["mcp"]
```

Those two fields are all a client needs. An `env` table is optional, and only
says where the tile cache lives: with
`env = { RASPUTIN_DATA = "/path/to/data" }` the server reads the cache in
`/path/to/data/cache`; without it, in `~/.cache/rasputin`.

The agent sees one tool per rasputin command (`sources`, `mesh`, `catchment`,
`fetch`, `fetch_stations`, ...; the names use underscores), each with the command's
options as its input schema, help text included, so an option added to the
CLI is a new tool argument with no further work. A call returns the exit code,
stdout and stderr, each output file as a `file://` link and, for `mesh`, the
run record and the `--stats` report as data; the mesh itself never travels
through the agent. A refused command comes back as an error carrying the
CLI's own `Error:` words. During a long run the server reports every line the
command prints to stderr and a heartbeat every 5 s while it is silent; cancelling the
call kills the run. One call runs at a time unless `--jobs N` says more. One
caveat: a `mesh` cancelled while writing a `.vtk` or `.ply` can leave a
truncated file under the asked name, so rerun before reading it; a `.gpkg`
goes through a `.part` file and a rename, so it cannot.

`sources` is the one to call first, and the server tells the agent so: it lists
every elevation and land-cover source the machine knows, local and remote,
with its coverage as a longitude/latitude box, cell size in metres, CRS,
height unit, description and credit, read from the files' headers (with
`bbox`, the share of an area each covers). The agent passes the chosen
source's id as `dem` to `mesh`, `catchment` or `station_catchments`.

Data already on the machine is registered by name in your sources file,
`~/.config/rasputin/sources.toml` (or `$XDG_CONFIG_HOME/rasputin/sources.toml`;
a `$RASPUTIN_DATA/sources.toml`, if one exists, is read instead). One table per
dataset:

```toml
[source.my-terrain]
kind = "local"
path = "/data/terrain-10m"    # a directory of GeoTIFF tiles, or one GeoTIFF
credit = "© The data's owner"
licence_note = "CC BY 4.0"
description = "Bare-ground terrain model, 10 m cells."
```

A relative `path` is taken from the file's directory. The CRS, cell size and
coverage come from the GeoTIFFs themselves; `z_unit = "metre"` names the
height unit for files that do not state it. `rasputin sources` on the
command line shows the same listing, and says which sources file it read.

## How it is built

- **A C++20 core** (`include/terrain/`, `src/`): exact geometric predicates,
  a constrained Delaunay triangulation, a noder that makes input lines meet
  cleanly, and the parallel refinement. The core never opens a file and knows
  nothing about coordinate systems. The only third-party C++ code is
  [detria](https://github.com/Kimbatt/detria), vendored in `lib/detria/`.
- **A Python layer** (`src_python/tin_engine/`): reading GeoTIFFs (tifffile),
  vector files and CRS (shapely, pyproj), fetching remote data,
  configuration (Pydantic) and the CLI (Typer). The two meet through a
  pybind11 module, `tin_engine._core`.

`project_structure.md` describes the layout and the rules between the layers.

## Where to read more

- [INSTALL.md](INSTALL.md): installing, testing, and getting data.
- [ROADMAP.md](ROADMAP.md): what has been built and what comes next.
- [docs/increments/](docs/increments/): one design record per step of the
  work, with its measurements.

## Publications

Bhattarai, B. C., Silantyeva, O., Teweldebrhan, A. T., Helset, S., Skavhaug,
O., and Burkhart, J. F.: Impact of Catchment Discretization and Imputed
Radiation on Model Response: A Case Study from Central Himalayan Catchment,
Water, 12, 2020. https://doi.org/10.3390/w12092339

Silantyeva, O., Skavhaug, O., Bhattarai, B. C., Helset, S., Tallaksen, L. M.,
Nordaas, M., and Burkhart, J. F.: Shyft and Rasputin: a toolbox for hydrologic
simulations on triangular irregular networks. https://doi.org/10.31223/X5CS95

## Licence

MIT; see [LICENSE](LICENSE). Third-party code and data credits are in
[NOTICE.md](NOTICE.md). Rasputin is developed by Expert Analytics.
