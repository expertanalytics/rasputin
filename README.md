# Rasputin

Rasputin can convert a point set of `(x, y, z)` coordinates to a triangulated
irregular network. Specifically, it has been developed to convert raster dems
(digital elevation models) into simplified triangulated surface meshes. The
`rasputin_store` program can read `GeoTIFF` files and construct surface
meshes in various formats. Run `rasputin_store --help` for more info.

It is also possible to compute the shade cast from a given, planar sun ray
vector. This shade is computed based on the cell center of the simplified
surface mesh.

## Implementation strategy

The heavy lifting in Rasputin is done by external software:
 * [CGAL](https://www.cgal.org/) is used for triangulation and simplification
   routines.
 * [pybind11](https://pybind11.readthedocs.io/en/stable/) is used to generate
   the Python wrappers.
 * [Pillow](https://python-pillow.org/) is used to read
   [GeoTIFF](https://en.wikipedia.org/wiki/GeoTIFF) files.
 * [Meshio](https://github.com/nschloe/meshio) is used to write results.
 * [Armadillo](http://arma.sourceforge.net/) for speedy arithmetics.
 * [date](https://github.com/HowardHinnant/date) for date and time on top of `chrono`.
 * [Catch2](https://github.com/catchorg/Catch2) for unit testing of the c++ code.


## Installation

Installing Rasputin is easy, as the C++ dependencies are header only. Simply, either install the dependencies using your system's package manager, or clone the source repository for each dependency locally on your computer and install them in such a way that CMake will find them.
For examaple, in some location where you have your source code, type:

```
git clone https://github.com/pybind/pybind11.git 
git clone https://gitlab.com/conradsnicta/armadillo-code.git 
git clone https://github.com/boostorg/geometry.git 
git clone https://github.com/catchorg/Catch2.git 
git clone https://github.com/CGAL/cgal.git 
```
For Howard Hinnant's `date` library to work, enter the rasputin source root directory and checkout the source under the `lib` folder:
```
cd lib
git clone https://github.com/HowardHinnant/date.git 
```

Rasputin does not aim at being backwards compatible with older compilers.
Hence, you will need something quite new. The following compilers are known to
work:
 * g++ 8.3
 * clang 11.0.0

Note that g++ 7 no longer works, due to the use of `<chrono>` from `stl`.
You can ensure that the right compiler is used for building Rasputin by setting the `CXX` environment variable.
For example, to use `g++` 8.x write the following in the terminal window:
```
export CXX=/usr/bin/g++-8
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
A relatively recent version of CMake will be needed, and `3.15.6` or newer is
known to work.

CGAL requires the two libraries [GMP](http://gmplib.org/) and
[MPFR](http://www.mpfr.org/) to be installed in order to work satisfactory. On
Ubuntu, these libraries can be installed the usual way by typing
```
sudo apt-get install libgmp-dev libmpfr-dev
```
or for Arch:
```
sudo pacman -Syy gmp mpfr
```
in a terminal window. Also, CGAL depends on [Boost](https://www.boost.org/),
see [here](https://doc.cgal.org/latest/Manual/installation.html#title21).

Additionally, you need Python 3.
Then, to install Rasputin, change to the Rasputin root source directory and run
```
pip3 install .
```
Or, if you prefer the old style:
```
python3 setup.py install
```


## Docker build
Take a look at the [Dockerfile](Dockerfile) to see how to setup required dependencies for a Debian system.

You can build rasputin and run tests by building the Docker image: `docker build . -t rasputin-test`


## Minimal Example
To test the installation run this for example in ipython:

```
import numpy as np
import pyproj
from rasputin.reader import Rasterdata
from rasputin.mesh import Mesh

def construct_rasterdata():
    raster = np.array([0, 0, 0, 
                       0, 1, 0, 
                       0, 0, 0], dtype=np.float32).reshape(3,3)
    cs = pyproj.CRS.from_epsg(32633)
    return Rasterdata(shape=(raster.shape[1], raster.shape[0]), x_min=0, 
                      y_max=20, delta_x=10, delta_y=10, array=raster,
                      coordinate_system=cs.to_proj4(), info={})

if __name__ == "__main__":
    rd = construct_rasterdata()
    mesh = Mesh.from_raster(data=rd)
    pts = mesh.points
    for face in mesh.faces:
        print("Face:", *[f'{fc:2d}' for fc in face])
        print(f"pts[{face[0]}]:", *[f'{pt:4.1f}' for pt in pts[face[0]]])
        print(f"pts[{face[1]}]:", *[f'{pt:4.1f}' for pt in pts[face[1]]])
        print(f"pts[{face[2]}]:", *[f'{pt:4.1f}' for pt in pts[face[2]]])
        print()
```

This should print out:
```
Face:  0  1  2
pts[0]: 10.0 10.0  1.0
pts[1]: 10.0 20.0  0.0
pts[2]:  0.0 20.0  0.0

Face:  0  2  3
pts[0]: 10.0 10.0  1.0
pts[2]:  0.0 20.0  0.0
pts[3]:  0.0 10.0  0.0

Face:  0  4  1
pts[0]: 10.0 10.0  1.0
pts[4]: 20.0 10.0  0.0
pts[1]: 10.0 20.0  0.0

Face:  4  5  1
pts[4]: 20.0 10.0  0.0
pts[5]: 20.0 20.0  0.0
pts[1]: 10.0 20.0  0.0

Face:  3  6  0
pts[3]:  0.0 10.0  0.0
pts[6]: 10.0  0.0  0.0
pts[0]: 10.0 10.0  1.0

Face:  3  7  6
pts[3]:  0.0 10.0  0.0
pts[7]:  0.0  0.0  0.0
pts[6]: 10.0  0.0  0.0

Face:  6  8  0
pts[6]: 10.0  0.0  0.0
pts[8]: 20.0  0.0  0.0
pts[0]: 10.0 10.0  1.0

Face:  8  4  0
pts[8]: 20.0  0.0  0.0
pts[4]: 20.0 10.0  0.0
pts[0]: 10.0 10.0  1.0
```
Congratulations! You just triangulated a small mountain.

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

The layout of this project follows the recommentation from an [excellent blog
post by Benjamin R.
Jack](http://www.benjack.io/2018/02/02/python-cpp-revisited.html). Both the
`CMakeExtension` and the `CMakeBuild` classes in `setup.py` are are taken from
his blog as well. Thanks!

## Use cases

Bhattarai, B. C., Silantyeva, O., Teweldebrhan, A. T., Helset, S., Skavhaug, O., and Burkhart, J. F.: Impact of Catchment Discretization and
Imputed Radiation on Model Response: A Case Study from Central Himalayan Catchment, Water, 12, 2020b; https://doi.org/10.3390/w12092339


Silantyeva, O.,  Skavhaug, O., Bhattarai, B.C., Helset, S., Tallaksen, L.M.,  Nordaas, M.,  and Burkhart, J.F.: Shyft and Rasputin: a toolbox for hydrologic simulations on triangular irregular networks. https://doi.org/10.31223/X5CS95
