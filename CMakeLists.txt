cmake_minimum_required (VERSION 3.9.1)
project (Rasputin)

add_definitions(-DARMA_DONT_USE_WRAPPER)

set(CMAKE_PREFIX_PATH  lib/CGAL lib/armadillo)
add_definitions(-DHEADER_ONLY="ON")

add_subdirectory(lib/pybind11)

find_package(CGAL REQUIRED)
if ( CGAL_FOUND )
    include( ${CGAL_USE_FILE} )
endif()

add_definitions("-std=c++17")
set(SOURCE_DIR "src/rasputin")
include_directories(${SOURCE_DIR})
include_directories(lib/armadillo/include)
pybind11_add_module(triangulate_dem "${SOURCE_DIR}/bindings.cpp")
