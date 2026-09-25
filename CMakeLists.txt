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
# The vendored detria predicates (lib/detria, MIT by election -- see its
# README). SYSTEM is load-bearing rather than tidy: the tree builds at
# -Wall -Wextra -Wpedantic -Werror and detria is unmodified C++17-era
# third-party code we have a written policy against patching, so its diagnostics
# must not become pressure to either patch it or weaken our own flags.
add_library(detria INTERFACE)
target_include_directories(detria SYSTEM INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/detria"
)

# The predicates backend: the only compiled target in the core, and the only
# place detria.hpp is ever included.
#
# `detria` is linked PRIVATE on purpose. That keeps lib/detria off this
# library's INTERFACE include path, so no consumer -- test target or future
# triangulation -- can resolve `#include <detria.hpp>` even by accident. The
# "detria.hpp lives in one translation unit" rule is therefore structural here,
# with tools/check_detria_boundary.py and the suite's preprocessor guard as the
# second and third lines of defence.
add_library(terrain_predicates STATIC
    src/predicates/detria_exact.cpp
)
target_link_libraries(terrain_predicates
    PUBLIC terrain_headers
    PRIVATE detria
)
target_compile_options(terrain_predicates PRIVATE
    -Wall -Wextra -Wpedantic -Werror
)

# The CDT backend: the second compiled target in the core, and the second and
# last place detria.hpp is ever included.
#
# `detria` is linked PRIVATE for the same reason it is on terrain_predicates:
# that keeps lib/detria off this library's INTERFACE include path, so no
# consumer can resolve `#include <detria.hpp>` even by accident, and the
# "detria.hpp lives in one translation unit per backend" rule stays structural.
add_library(terrain_cdt STATIC
    src/cdt/detria_backend.cpp
)
target_link_libraries(terrain_cdt
    PUBLIC terrain_headers
    PRIVATE detria
)
target_compile_options(terrain_cdt PRIVATE
    -Wall -Wextra -Wpedantic -Werror
)

option(RASPUTIN_BUILD_PYTHON "Build the Python extension module" OFF)
if(RASPUTIN_BUILD_PYTHON)
    find_package(pybind11 CONFIG REQUIRED)

    pybind11_add_module(_core bindings/core.cpp)
    # terrain_cdt for DetriaBackend::triangulate, and terrain_predicates
    # because build_pslg validates under DefaultKernel -- the same pair, for
    # the same two reasons, that tests/cpp/CMakeLists.txt's add_terrain_cdt_test
    # links. Omitting the second links on macOS (a MODULE is allowed undefined
    # symbols) and fails at import, so the C++ build is not where it is caught.
    target_link_libraries(_core PRIVATE terrain_headers terrain_cdt terrain_predicates)
    # refine starts threads (include/terrain/parallel_util/chunks.hpp).
    find_package(Threads REQUIRED)
    target_link_libraries(_core PRIVATE Threads::Threads)

    # Both static libraries are linked into a SHARED module here, so their
    # objects must be position-independent. On ELF the linker refuses outright
    # -- "relocation R_X86_64_PC32 ... can not be used when making a shared
    # object; recompile with -fPIC" -- while Mach-O does not care, so this is
    # invisible on macOS and fails every Linux leg. Set on the libraries rather
    # than globally: it is a property of being linked into _core, and the ctest
    # executables that also consume them need nothing.
    set_property(TARGET terrain_cdt terrain_predicates
                 PROPERTY POSITION_INDEPENDENT_CODE ON)
    target_compile_options(_core PRIVATE -Wall -Wextra -Wpedantic -Werror)

    install(TARGETS _core LIBRARY DESTINATION tin_engine)
endif()

option(RASPUTIN_BUILD_TESTS "Build C++ unit tests" ON)
if(RASPUTIN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests/cpp)
endif()
