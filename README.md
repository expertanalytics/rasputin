# Rasputin

Rasputin can convert a point set of (x,y,z) coordinates to a triangulatd irregular network. Specifically, it has been developed to convert raster dems (digital elevation models) into a simplified triangulated surface mesh. The `rasputin_triangulate` program can read `GeoTIFF` files and construct surface meshes in the `off` format.

## Implementation strategy

All the heavy lifting in Rasputin is made by external software. The triangulation and simplication routines are done in `CGAL`, the processing of the GeoTIFF files in handled by `Pillow`. Finally, the wrapper code is built by `pybind11`.

## Installation

Installing Rasputin is easy, as the C++ dependencies are header only. Simply download and unpack the source code for `pybind11` and `CGAL` and place them under the `lib` directory using the names `pybind11` and `CGAL`, respectively:
```
cd <rasputin_directory>/lib
wget https://github.com/pybind/pybind11/archive/v2.2.3.tar.gz
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13/CGAL-4.13.tar.xz
tar zxf v2.2.3.tar.gz && mv pybind11-2.2.3 pybind11
tar xjf CGAL-4.13.tar.xz && mv CGAL-4.13 CGAL
```

Additionally, you need Python 3, a modern compiler supporting C++11, and CMake. Then, to install Rasputin, change to the Rasputin root source directory and run
```
python3 setup.py install
```
