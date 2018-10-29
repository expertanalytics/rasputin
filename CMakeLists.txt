cmake_minimum_required (VERSION 3.11.3)
project (Rasputin)

set(CMAKE_PREFIX_PATH  lib/CGAL)
add_definitions(-DHEADER_ONLY="ON")

add_subdirectory(lib/pybind11)

find_package(CGAL)
if ( CGAL_FOUND )
    include( ${CGAL_USE_FILE} )
endif()

add_definitions("-std=c++17")
set(SOURCE_DIR "src/rasputin")
include_directories(${SOURCE_DIR})
pybind11_add_module(triangulate_dem "${SOURCE_DIR}/bindings.cpp")
