cmake_minimum_required (VERSION 3.9.1)
project (Rasputin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make cmake search for installations in 'lib' before looking elsewhere
set(CMAKE_PREFIX_PATH  lib/CGAL ${CMAKE_PREFIX_PATH})

# Include pybind11 CMakeLists.txt in 'lib/pybind11' if it exists
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/pybind11")
    add_subdirectory(lib/pybind11)
else()
    find_package(pybind11 REQUIRED)
endif()

# CGAL repository root has config file for header-only use
find_package(CGAL REQUIRED)
if (CGAL_FOUND)
    include(${CGAL_USE_FILE})
endif()

# Include Armadillo source in 'lib/armadillo' if it exists
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/armadillo")
    set(ARMADILLO_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/armadillo/include")
else()
    find_package(Armadillo REQUIRED)  # This fails if the wrapper library is not found
endif()

# Use Armadillo as a header-only library
add_definitions(-DARMA_DONT_USE_WRAPPER)

# Need to link to BLAS libraries when not using armadillo wrappers
find_package(BLAS REQUIRED)
find_package(LAPACK)

set(SOURCE_DIR "src/rasputin")

include_directories(${SOURCE_DIR} ${ARMADILLO_INCLUDE_DIRS})
pybind11_add_module(triangulate_dem "${SOURCE_DIR}/bindings.cpp")
target_link_libraries(triangulate_dem PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
