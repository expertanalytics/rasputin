cmake_minimum_required (VERSION 3.9.1)
project (Rasputin)

set(CMAKE_CXX_STANDARD 20)
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

# Check if date library is installed in 'lib/date'
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/date")
    set(DATE_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/date/include")
else()
    find_package(date REQUIRED)
    get_target_property(DATE_INCLUDE_DIR date::date INTERFACE_INCLUDE_DIRECTORIES)
endif()

add_subdirectory(cpp_test)

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

pybind11_add_module(triangulate_dem "${SOURCE_DIR}/bindings.cpp")
target_include_directories(triangulate_dem PRIVATE
    ${SOURCE_DIR} ${ARMADILLO_INCLUDE_DIRS} ${DATE_INCLUDE_DIR})
target_link_libraries(triangulate_dem PRIVATE
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
