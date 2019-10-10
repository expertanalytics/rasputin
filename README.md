# Rasputin

Rasputin can convert a point set of `(x, y, z)` coordinates to a triangulated
irregular network. Specifically, it has been developed to convert raster dems
(digital elevation models) into simplified triangulated surface meshes. The
`rasputin_triangulate` program can read `GeoTIFF` files and construct surface
meshes in various formats.

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

Rasputin does not aim at being backwards compatible with older compilers.
Hence, you will need something quite new. The following compilers are known to
work:
 * g++ 8.3 
 * clang 11.0.0

Note that g++ 7 no longer works, due to the use of `<chrono>` from `stl`.

## Installation

Installing Rasputin is easy, as the C++ dependencies are header only. Simply
download and unpack the source code for `pybind11` and `CGAL` and place them
under the `lib` directory using the names `pybind11` and `CGAL`, respectively.

As of November 2018, the latest release of `pybind11` generates quite a lot of deprecation warnings for Python3.7. In this case, pulling the source directly from the master branch at github is recommended:
```
cd <rasputin_directory>/lib
git clone git@github.com:pybind/pybind11.git
```

Alternatively, use the latest released sources:
```
cd <rasputin_directory>/lib
wget https://github.com/pybind/pybind11/archive/v2.2.3.tar.gz
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13/CGAL-4.13.tar.xz
wget http://sourceforge.net/projects/arma/files/armadillo-9.200.7.tar.xz
git clone git@github.com:HowardHinnant/date.git
git clone git@github.com:catchorg/Catch2.git
tar xf v2.2.3.tar.gz && mv pybind11-2.2.3 pybind11
tar xf CGAL-4.13.tar.xz && mv CGAL-4.13 CGAL
tar xf armadillo-9.200.7.tar.xz && mv mv armadullo-9.200.7 armadillo
```

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

Additionally, you need Python 3, a modern compiler supporting C++17, and CMake.
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
from rasputin import triangulate_dem

vs = triangulate_dem.point3_vector([[1,0,0], [0,1,0], [0,0,0], [0.25,0.25,1]])
points, faces = triangulate_dem.lindstrom_turk_by_ratio(vs, 2.0)
for f in faces:
   print(f)
```

This should print out:
```
>> (3, 1, 2)
>> (0, 1, 3)
>> (0, 3, 2)
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
