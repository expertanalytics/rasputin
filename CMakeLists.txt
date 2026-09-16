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

option(RASPUTIN_BUILD_TESTS "Build C++ unit tests" ON)
if(RASPUTIN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests/cpp)
endif()
