cmake_minimum_required(VERSION 3.24)
project(rasputin
    VERSION 0.2.0
    DESCRIPTION "Parallel TIN engine for terrain meshing"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

add_library(terrain_headers INTERFACE)
target_include_directories(terrain_headers INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_compile_features(terrain_headers INTERFACE cxx_std_20)

# The Python extension is built by scikit-build-core when installing the
# wheel; it is off by default so a plain `cmake ..` still configures only the
# header library and the C++ tests.
option(RASPUTIN_BUILD_PYTHON "Build the Python extension module" OFF)
if(RASPUTIN_BUILD_PYTHON)
    find_package(pybind11 CONFIG REQUIRED)

    pybind11_add_module(_core bindings/core.cpp)
    target_link_libraries(_core PRIVATE terrain_headers)
    target_compile_options(_core PRIVATE -Wall -Wextra -Wpedantic -Werror)

    install(TARGETS _core LIBRARY DESTINATION tin_engine)
endif()

option(RASPUTIN_BUILD_TESTS "Build C++ unit tests" ON)
if(RASPUTIN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests/cpp)
endif()
