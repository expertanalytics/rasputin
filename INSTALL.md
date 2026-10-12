# Installing rasputin

Rasputin is not on PyPI yet. You install it from a clone, and `pip` (or `uv`)
compiles the C++ core on the way. It takes a minute or less.

Steps marked **(not verified)** were not run when this file was written.
Every other install step was run on 2026-09-29 on macOS arm64 (AppleClang 21,
CMake 4.4, Python 3.13), in a fresh clone of the local repository (not of
GitHub), with both `uv` and `pip`; or it is what CI runs. The data sections
were checked again on 2026-10-11.

## Platforms

| platform | what CI checks (`.github/workflows/main.yaml`) |
|---|---|
| Linux x86_64 (`ubuntu-latest`) | C++ build and tests, sanitizers, and the Python suite on Python 3.12, 3.13 and 3.14 |
| macOS arm64 (`macos-latest`) | C++ build and tests only |

The Python suite is not run on macOS in CI, but it passes there locally (see
"Running the tests"). Windows is not supported.

## Prerequisites

- **A C++20 compiler.** CI builds with GCC 13.3 (Ubuntu) and AppleClang 21
  (macOS); the compiler must provide the C++20 `<chrono>` calendar types.
- **CMake 3.24 or newer** (`cmake_minimum_required` in `CMakeLists.txt`).
- **Python 3.12 or newer** (`requires-python` in `pyproject.toml`), with its
  development headers.
- **git**, and network access during the build: pip fetches the build tools
  (scikit-build-core, pybind11), and the C++ test build fetches Catch2.

That is the whole list: everything else the install needs comes from PyPI
during `pip install`, and the Python packages bring their own native libraries.

On **macOS**, the Xcode Command Line Tools give you the compiler
**(not verified: this machine already had them)**:

```sh
xcode-select --install
brew install cmake
```

On **Ubuntu or Debian** **(not verified)**:

```sh
sudo apt-get install build-essential cmake git python3-dev python3-venv
```

If your distribution's Python is older than 3.12, use `uv` (below), which
downloads a Python for you.

## Install

Clone, then install into a virtual environment. With
[uv](https://docs.astral.sh/uv/) (`brew install uv`, or see its site):

```sh
git clone https://github.com/expertanalytics/rasputin.git
cd rasputin
uv venv --python 3.13
uv pip install ".[codecs]"
source .venv/bin/activate
rasputin version
```

Or with plain `pip`:

```sh
git clone https://github.com/expertanalytics/rasputin.git
cd rasputin
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[codecs]"
rasputin version
```

`rasputin version` prints the version (`0.2.0.dev0` at the time of writing).
Then try the quick example in [README.md](README.md).

The build is driven by scikit-build-core, which runs the project's
`CMakeLists.txt` with the Python extension on and the C++ tests off, and
installs the compiled module as `tin_engine._core`.

### Extras

| extra | what it adds | when you need it |
|---|---|---|
| `codecs` | `imagecodecs` | To read LZW-compressed GeoTIFFs, or those with the floating-point predictor. Kartverket's DTM10 tiles and the committed tile `tests/fixtures/dem_archive/7908_3_10m_z33.tif` need it. Without it, rasputin reads uncompressed and Deflate files and refuses the rest with a message naming the extra. |
| `viewer` | `vtk` (100-140 MB) | Only for the tests that read `.vtk` output back through VTK. You do not need it to view meshes; install ParaView for that. |
| `dev` | pytest, hypothesis, mypy, ruff | To run the tests and the lint gates. |

Combine them as `".[dev,codecs]"`.

### For development

An editable install lets you change the Python code without reinstalling:

```sh
uv pip install -e ".[dev,codecs]"      # or: python -m pip install -e ".[dev,codecs]"
```

Changes to the C++ code are **not** picked up automatically: reinstall
(`uv pip install -e ".[dev,codecs]"` again) after editing anything under
`include/`, `src/` or `bindings/`.

### Bounds checks (on by default)

A normal install builds the C++ core with the standard library's own bounds
checks on: `_LIBCPP_HARDENING_MODE_FAST` with libc++ (macOS),
`_GLIBCXX_ASSERTIONS` with libstdc++ (Linux). An out-of-range index into a
`std::vector` or `std::span` then stops the process at the faulting call (a
crash, not a Python exception) instead of reading or writing the memory next to
it and carrying on. The cost, measured: about 5 % on refinement (4 to 9 % across
thread counts) and 2 to 3 % on a whole `rasputin mesh` run, with identical
meshes ([the measurement](docs/benchmarks/2026-10-03/release-hardening/README.md)).

`rasputin version` says which build is installed, on its second line:
`bounds checks: on (libc++ fast)` (or `on (libstdc++ assertions)`), or
`bounds checks: off`. `rasputin mesh --stats` writes the same line into its
report.

For a heavy run where those percent matter, build an unchecked copy **in a
second virtual environment**, so the everyday one never changes mode:

```sh
uv venv .venv-unchecked
uv pip install --python .venv-unchecked ".[codecs]" -C cmake.define.RASPUTIN_HARDENING=OFF
.venv-unchecked/bin/rasputin version     # bounds checks: off
```

With plain `pip`, the install line is
`python -m pip install ".[codecs]" -C cmake.define.RASPUTIN_HARDENING=OFF`.
For the C++ tests, or any plain CMake build:
`cmake -S . -B build-unchecked -DCMAKE_BUILD_TYPE=Release -DRASPUTIN_HARDENING=OFF`.

The unchecked environment is a non-editable install, so it too needs
reinstalling after C++ changes (see Troubleshooting).

Another libc++ mode (`extensive`, `debug`) is chosen by setting the option OFF
and passing the mode yourself, e.g.
`-DCMAKE_CXX_FLAGS=-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE`;
with the option ON as well, the two definitions collide. This works only for a
build without the C++ tests (`pip install`, or `-DRASPUTIN_BUILD_TESTS=OFF`),
whose compile-time guard refuses an OFF option beside a hardening mode.

## Running the tests

Python, from the repository root with the `dev` extra installed:

```sh
pytest
```

This runs `tests/python/` and fails if line coverage drops below 85 %. On
macOS arm64 with Python 3.13 and `".[dev,codecs]"` it takes under two
minutes. Every test should pass or be skipped, with none failing. The skips
are expected: a few tests need the `viewer`
extra, a few test the refusal when `codecs` is absent, and a few need the full
national datasets in `../rasputin_data/` (a directory next to the clone) and
skip when it is not there.

C++ (no Python needed; CMake fetches Catch2 v3.6.0):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

The static checks CI runs are listed in `CLAUDE.md`, section 4.

## Getting data

### Committed test data

`tests/fixtures/` holds small real data you can mesh at once:

- `dem_archive/7908_3_10m_z33.tif`: one Kartverket DTM10 tile, 50 km square,
  10 m cells, EPSG:25833 (needs `codecs`).
- `dtm10/seam/` and `dtm10/lattices/`: small windows of neighbouring DTM10
  tiles, for `--dem DIR`.
- `corine/clc2018_7908_3.gpkg`: CORINE Land Cover 2018 polygons over that
  tile, for `--features`. Its credits are in `corine/NOTICE`.

### Remote DEMs: `rasputin fetch`

Two DEMs rasputin fetches itself, into a tile cache (`$RASPUTIN_DATA/cache`
by default, or `--cache DIR`):

```sh
rasputin fetch glo30 --domain catchment.geojson   # --dry-run prints the plan
rasputin mesh --dem glo30 --domain catchment.geojson --out-crs EPSG:32723 \
    --tolerance 5 --out catchment.vtk             # offline from here
```

`glo30` is the Copernicus GLO-30 (30 m, worldwide) and `anadem-v1` the ANADEM
terrain model of South America. Both are in geographic coordinates, so `mesh`
needs `--out-crs`; its refusal suggests one. The source's cache directory
holds a `NOTICE.txt` with the distributor's credit and licence, and every mesh
made from it carries the credit line.

### Norway: Kartverket DTM10

Kartverket's 10 m terrain model of mainland Norway is open data under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), credit
"© Kartverket" (its Geonorge record and Kartverket's terms of use, read
2026-10-11). It is on offer in two places, and no login is needed for either:

- **The current model, from [hoydedata.no](https://hoydedata.no/LaserInnsyn2/).**
  Open the "Nedlasting" (download) tab, then "Landsdekkende" (nationwide),
  "UTM sone 33", and the "DTM10" entry (the click path is from the project's own
  2026-09-25 download; the download endpoints, not the web page, were
  checked by script): one zip of about 19.6 GB named `DTM10_UTM33_<date>.zip`,
  re-exported every few weeks. The 2026-09-25 export unpacked to 254 tiles of
  the national 50 km grid as GeoTIFF (`7908_3_10m_z33.tif` and so on:
  5051 × 5051 float32 cells, LZW, EPSG:25833, NoData −32767; about 20 GB);
  a later export's tile list was not checked. The committed tile and the
  project's own runs come from here.
- **The 2013-based model, from Geonorge**, one zip per tile:
  `https://nedlasting.geonorge.no/geonorge/Basisdata/DTM10UTM33/TIFF/Basisdata_<cell>_Celle_25833_DTM10UTM33_TIFF.zip`,
  `<cell>` from `6400-1` to `7910-3`; all 254 are listed in the Atom feed
  `https://nedlasting.geonorge.no/geonorge/ATOM-feeds/DTM10UTM33_AtomFeedTIFF.xml`,
  7.5 GB in all. Its tiles are 5041 × 5041 and sit half a cell off
  hoydedata.no's lattice (checked on tile 6400-1), so keep the two sets in
  separate directories. Kartverket says this set is no longer updated.

Unpack the tiles into one directory and give that directory to `--dem`; use
`--bbox XMIN YMIN XMAX YMAX` or `--domain` to mesh part of it. The tiles need
the `codecs` extra. A `rasputin fetch dtm10` is possible, since both routes
are plain HTTP: `docs/research/2026-10-11-dtm10-fetcher.md`.

The project's own runs keep data in a directory next to the clone, for example
`../rasputin_data/DTM10_UTM33_20260925/`. Nothing reads a fixed location: every
path is given on the command line.

### Norway: NVE's reference gauging stations

```sh
rasputin fetch-stations nve-hrd --out-dir ../rasputin_data/nve-hrd
```

fetches the stations of NVE's Hydrological Reference Dataset, their catchment
polygons and the river lines round them, with a `NOTICE.txt` (licence NLOD,
"Kilde: NVE"); files already there are kept unless `--refresh`. They are the
inputs of `rasputin station-catchments`.

### Europe: CORINE Land Cover **(not verified)**

CORINE Land Cover 2018 comes from the
[Copernicus Land Monitoring Service](https://land.copernicus.eu/) (EEA), free
with a registered account. Download the vector GeoPackage
(`U2018_CLC2018_V2020_20u1.gpkg`, about 8 GB) and pass it as

```sh
--features U2018_CLC2018_V2020_20u1.gpkg \
--features-layer U2018_CLC2018_V2020_20u1 --features-map corine
```

Rasputin reads only the polygons near your domain, through the file's spatial
index. If you publish results, credit the data as `corine/NOTICE` shows.

## Troubleshooting

- **A GeoTIFF is refused with "cannot be decoded as installed; the `codecs`
  extra (imagecodecs) may decode it".** Install the extra:
  `pip install ".[codecs]"`.
- **The build fails on `<chrono>` calendar types, `concept` or other C++20
  features.** Your compiler is too old. On macOS, check
  `c++ --version`; if `/usr/bin/c++` points at an old toolchain, run
  `sudo xcode-select --switch /Library/Developer/CommandLineTools`. You can pick
  a compiler with `export CXX=...` before installing; for GCC, `CXX` must be
  `g++`, not `gcc`.
- **`CMake 3.24 or higher is required`.** Install a newer CMake, from your
  package manager or with `pip install cmake` **(not verified)**.
- **`Python.h: No such file or directory`** (Linux). Install your Python's
  development package (`python3-dev`), or use a Python from `uv`.
- **C++ changes have no effect.** An editable install does not rebuild the
  extension; reinstall. The same holds for an unchecked environment
  ("Bounds checks" above): reinstall it with the `-C` setting.
