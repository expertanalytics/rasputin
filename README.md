# Rasputin

Rasputin can convert a point set of `(x, y, z)` coordinates to a triangulated irregular network. Specifically, it has been developed to convert raster dems (digital elevation models) into simplified triangulated surface meshes. The `rasputin_triangulate` program can read `GeoTIFF` files and construct surface meshes in the `off` format.

## Implementation strategy

All the heavy lifting in Rasputin is done by external software. The triangulation and simplification routines are done in `CGAL` and the wrapper code is built by `pybind11`.  Finally, the processing of the GeoTIFF files is handled by the Python module `Pillow`.

## Installation

Installing Rasputin is easy, as the C++ dependencies are header only. Simply download and unpack the source code for `pybind11` and `CGAL` and place them under the `lib` directory using the names `pybind11` and `CGAL`, respectively:
```
cd <rasputin_directory>/lib
wget https://github.com/pybind/pybind11/archive/v2.2.3.tar.gz
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13/CGAL-4.13.tar.xz
tar zxf v2.2.3.tar.gz && mv pybind11-2.2.3 pybind11
tar xjf CGAL-4.13.tar.xz && mv CGAL-4.13 CGAL
```
An alternative for `pybind11`, currently resulting in a lot fewer compiler warnings, is to clone the master branch from git:
```
cd <rasputin_directory>/lib
git clone git@github.com:pybind/pybind11.git
```

Additionally, you need Python 3, a modern compiler supporting C++11, and CMake. Then, to install Rasputin, change to the Rasputin root source directory and run
```
python3 setup.py install
```

## Data

The Rasputin `GeoTIFF` data reader has only been used on one file, and hence, should not be expected to work on a wide variety of files. However, this data set can be downloaded here: http://blog.mastermaps.com/2016/09/creating-tin-from-raster-dem.html

## Acknowledges

The layout of this project follows the recommentation from an excellent blog post by Benjamin R. Jack: http://www.benjack.io/2018/02/02/python-cpp-revisited.html. Also, in `setup.py` both the `CMakeExtension` and the `CMakeBuild` are taken from this blog.
