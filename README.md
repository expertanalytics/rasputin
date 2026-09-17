# Rasputin

Rasputin can convert a point set of `(x, y, z)` coordinates to a triangulated
irregular network. Specifically, it has been developed to convert raster dems
(digital elevation models) into simplified triangulated surface meshes. The
`rasputin_store` program did this in the CGAL era and does not exist in the
post-CGAL tree; the installed CLI is `rasputin`, which currently exposes only
`version` while the core is rebuilt. See `docs/increments/` for what has landed.

It is also possible to compute the shade cast from a given, planar sun ray
vector. This shade is computed based on the cell center of the simplified
surface mesh.

## Implementation strategy

Rasputin is being rebuilt on a CGAL-free and GDAL-free stack: a C++20 core
with async Python bindings. The C++ dependencies are header-only.

### C++ core
 * [pybind11](https://pybind11.readthedocs.io/en/stable/) generates the Python
   wrappers.
 * A constrained Delaunay triangulation library, used exactly once to build the
   initial triangulation — everything downstream of it is owned in-tree.
   Detria (MIT, header-only, C++20) is the intended choice, with poly2tri
   (BSD-2) as fallback; see `parallel_refinement.md` for the comparison.
 * Exact `orient2d` and `incircle` predicates via a filtered kernel over a
   vendored backend (`lib/detria/`, MIT, pinned). Not header-only: the backend
   is confined to one translation unit, `src/predicates/detria_exact.cpp`.
 * [Catch2 v3](https://github.com/catchorg/Catch2) for unit tests. Property
   tests use Catch2 `GENERATE` over a seeded range rather than rapidcheck —
   each property suite records why at the top of the file.

Parallelism uses the standard library by default (`<thread>`, `<atomic>`,
`std::execution`); TBB or OpenMP can be opted in via a CMake flag.

### Python layer
 * [Pydantic V2](https://docs.pydantic.dev/) for models and validation.
 * [Typer](https://typer.tiangolo.com/) for the CLI.
 * [Shapely](https://shapely.readthedocs.io/) and
   [PyProj](https://pyproj4.github.io/pyproj/) for geometry and CRS handling.
 * [NumPy](https://numpy.org/) for array interchange across the binding layer.

Raster reading deliberately avoids Rasterio, which wraps GDAL. A pure-Python
reader will be introduced under `src_python/tin_engine/io/`.

### Removed
CGAL, GMP and MPFR (triangulation and simplification move to the CDT wrapper
plus in-tree refinement), Armadillo, Pillow, Meshio, and Howard Hinnant's
`date` library (superseded by C++20 `<chrono>`). The pre-migration tree is
preserved under `legacy/` and still references several of these.


## Installation

The C++ dependencies are header-only, and Catch2 is fetched automatically by
CMake via `FetchContent`, so no manual checkout is needed.

Create the Python environment:

```
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"
```

Building the Python extension needs no extra step: `pybind11` is declared in
`[build-system].requires`, and scikit-build-core invokes CMake for you.
Rasputin does not aim at being backwards compatible with older compilers.
The build requires C++20. Date and time handling uses the C++20 `<chrono>`
calendar types (`sys_days`, `year`/`month`/`day`) directly, so the compiler must
provide them. The following are verified to build and pass the test suite:
 * AppleClang 21.0.0
 * g++ 16.2.0

All date handling is UTC; no timezone database is required, so libc++ is fine
despite not yet shipping `std::chrono::zoned_time`.

On macOS, note that `/usr/bin/c++` dispatches through `xcode-select`. If it
resolves to an old toolchain, point it at the Command Line Tools:
```
sudo xcode-select --switch /Library/Developer/CommandLineTools
```
You can ensure that the right compiler is used for building Rasputin by setting
the `CXX` environment variable, for example:
```
export CXX=/opt/homebrew/bin/g++-16
```
If you are using gcc, make sure that `CXX` points to `g++` and not `gcc`.

Rasputin is build using [CMake](https://cmake.org). On Ubuntu, CMake can be installed with the command
```
sudu apt-get install cmake
```
or on Arch,
```
sudo pacman -S cmake
```
A relatively recent version of CMake is needed; the build declares a minimum of
`3.24`.

Additionally, you need Python 3.11 or newer.
Then, to install Rasputin, change to the Rasputin root source directory and run
```
pip install .
```
or, for a development install with the test and lint tooling, use the
virtualenv shown above.


## Minimal Example

The post-CGAL core is under construction: the meshing pipeline is not wired up
yet, so the installable surface is currently the geometry primitives and the
CLI. To check that the compiled extension imported correctly, run this in
ipython:

```
from tin_engine import Point2, Point3, dot, cross

a = Point2(3.0, 4.0)
b = Point2(1.0, 2.0)

print("dot(a, b)  =", dot(a, b))
print("cross(a, b) =", cross(a, b))
print("cross(x, y) =", cross(Point3(1.0, 0.0, 0.0), Point3(0.0, 1.0, 0.0)))
```

This should print out:
```
dot(a, b)  = 11.0
cross(a, b) = 2.0
cross(x, y) = Point3(0, 0, 1)
```

The CLI is installed as `rasputin`:
```
rasputin version
```

The legacy CGAL-based pipeline that used to be demonstrated here is archived
under `legacy/` for reference during the port. It is not packaged and not
importable from an installed rasputin.

## Data

High quality DTM data for Norway can be downloaded from free [here](https://hoydedata.no/LaserInnsyn/).
Choose "Nedlasting" from the left hand side of the map, and choose "Landsdekkende", check "UTM-sone 33"
and finally click DTM10. Download and unpack in, for instance, `$HOME/rasputin_data/dem_archive`, and
`export RASPUTIN_DATA_DIR=$HOME/rasputin_data`.

It is possible to include land cover types in your triangulation, through the 
[GlobCover dataset](http://due.esrin.esa.int/page_globcover.php) from ESA. It is a raster based 
300m (approx) resolution data set that contains 23 different land cover types. 
Download the data set and unpack it in `$RASPUTIN_DATA_DIR/globcov` to access the land types using
the `rasputin.globcov_repository.GlobCovRepository` class.

## Acknowledges

The original layout of this project followed the recommendation from an
[excellent blog post by Benjamin R.
Jack](http://www.benjack.io/2018/02/02/python-cpp-revisited.html), and both the
`CMakeExtension` and the `CMakeBuild` classes were taken from his blog as well.
They lived in the `setup.py` that the post-CGAL migration replaced with
`pyproject.toml`, so they are no longer in the tree -- but the debt stands.
Thanks!

## Use cases

Bhattarai, B. C., Silantyeva, O., Teweldebrhan, A. T., Helset, S., Skavhaug, O., and Burkhart, J. F.: Impact of Catchment Discretization and
Imputed Radiation on Model Response: A Case Study from Central Himalayan Catchment, Water, 12, 2020b; https://doi.org/10.3390/w12092339


Silantyeva, O.,  Skavhaug, O., Bhattarai, B.C., Helset, S., Tallaksen, L.M.,  Nordaas, M.,  and Burkhart, J.F.: Shyft and Rasputin: a toolbox for hydrologic simulations on triangular irregular networks. https://doi.org/10.31223/X5CS95
