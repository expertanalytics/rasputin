cmake_minimum_required (VERSION 3.9.1)
project (Rasputin)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_VERBOSE_MAKEFILE ON)

set(ENABLE_TESTS ON)

# Default location for dependencies
set(RASPUTIN_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# Make cmake search for installations in 'lib' before looking elsewhere
# (used to find CGAL)


# Find Rasputin dependencies
# --------------------------
# All dependencies (except BLAS) are header only and can be installed under lib
# with git clone. Dependencies installed here should be prefererred over
# installations in other locations, e.g. /usr/local/include

set(RASPUTIN_DEPENDENCIES) # Append to list as they are found


# Pybind11
# --------
if (EXISTS "${RASPUTIN_LIB_DIR}/pybind11")
    add_subdirectory("${RASPUTIN_LIB_DIR}/pybind11")
else()
    find_package(pybind11 REQUIRED)
endif()


# Boost::geometry
# ---------------
if (EXISTS "${RASPUTIN_LIB_DIR}/geometry")
    add_library(geometry INTERFACE)
    target_include_directories(geometry PRIVATE INTERFACE
        "${RASPUTIN_LIB_DIR}/geometry/include")
    list(APPEND RASPUTIN_DEPENDENCIES geometry)
    add_definitions(-DHAS_BOOST)
else()
    find_package(Boost REQUIRED)
    # Import into scope
    add_library(Boost::geometry INTERFACE IMPORTED)
    list(APPEND RASPUTIN_DEPENDENCIES Boost::geometry)
    add_definitions(-DHAS_BOOST)
endif()


# date
# ----
# TODO: Check CXX version if this is really needed?
if (EXISTS "${RASPUTIN_LIB_DIR}/date")
    add_library(date INTERFACE)
    target_include_directories(date INTERFACE
        "${RASPUTIN_LIB_DIR}/date/include")
else()
    find_package(date REQUIRED)
    add_library(date INTERFACE IMPORTED)
endif()
list(APPEND RASPUTIN_DEPENDENCIES date)


# CGAL
# ----
# Required for now but should be optional in the future
set(CMAKE_PREFIX_PATH  "${RASPUTIN_LIB_DIR}/cgal" ${CMAKE_PREFIX_PATH})
find_package(CGAL REQUIRED)
if (CGAL_FOUND)
    add_definitions(-DHAS_CGAL) # For future use when CGAL is optional
    list(APPEND RASPUTIN_DEPENDENCIES CGAL::CGAL)
endif()
# NOTE: Calling`include(${CGAL_USE_FILE})` will modify INCLUDE_DIRECTORIES,
# with the consquence that directories added by CGAL will be searched before
# anything else, potentiall messing up the include order.


# Armadillo
# ---------
if (EXISTS "${RASPUTIN_LIB_DIR}/armadillo")
    set(ARMADILLO_INCLUDE_DIRS "${RASPUTIN_LIB_DIR}/armadillo/include")
else()
    # NOTE: This fails if the wrapper library is not found
    find_package(Armadillo REQUIRED)
endif()
add_library(Armadillo INTERFACE)
target_include_directories(Armadillo PRIVATE INTERFACE
    ${ARMADILLO_INCLUDE_DIRS})
add_definitions(-DARMA_DONT_USE_WRAPPER) # Header-only usage


# BLAS
# ----
# Need to link to BLAS libraries when not using armadillo wrappers
find_package(BLAS REQUIRED)
find_package(LAPACK)


# Header-only Rasputin library
# ----------------------------
set(RASPUTIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/rasputin")
add_library(rasputin INTERFACE)
target_include_directories(rasputin PRIVATE INTERFACE ${RASPUTIN_SOURCE})

if (ENABLE_TESTS)
    add_subdirectory(cpp_test)
endif()


# Python bindnigs
# ---------------
# Use pybind11 as a system library to avoid problems with include order
pybind11_add_module(triangulate_dem MODULE SYSTEM
    "${RASPUTIN_SOURCE}/bindings.cpp")
# Need to link BLAS for header-only Aramdillo
target_link_libraries(triangulate_dem PRIVATE
    ${RASPUTIN_DEPENDENCIES}
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
