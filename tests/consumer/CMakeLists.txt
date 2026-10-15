# A plugin's build against the installed core, and nothing else: no path into
# this repository's include/ or bindings/. It fails when a promised header is
# not installed or the CMake config is wrong.
# Spec: docs/increments/51-open-core.md, section 4.4 and the consumer-build
# step of section 9.1.
#
#   cmake -S tests/consumer -B build-consumer \
#         -Drasputin_ROOT=<the installed tin_engine directory> \
#         -Dpybind11_DIR="$(python -m pybind11 --cmakedir)"
#   cmake --build build-consumer && ./build-consumer/consumer_headers
#
# EXPECT_HARDENING: the core was built with RASPUTIN_HARDENING ON (the wheel's
# default), so the imported target must carry the standard library's checks.
cmake_minimum_required(VERSION 3.24)
project(rasputin_consumer LANGUAGES CXX)

option(EXPECT_HARDENING "The installed core was built with RASPUTIN_HARDENING ON" ON)

find_package(rasputin CONFIG REQUIRED)
find_package(Threads REQUIRED)

# The config alone must give C++20: nothing here asks for it.
add_executable(consumer_headers headers.cpp)
target_link_libraries(consumer_headers PRIVATE rasputin::terrain_headers Threads::Threads)
target_compile_options(consumer_headers PRIVATE -Wall -Wextra -Wpedantic -Werror)
if(EXPECT_HARDENING)
    target_compile_definitions(consumer_headers PRIVATE CONSUMER_EXPECT_HARDENING)
endif()

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(consumer_ext bindings.cpp)
target_link_libraries(consumer_ext PRIVATE rasputin::terrain_headers)
target_compile_options(consumer_ext PRIVATE -Wall -Wextra -Wpedantic -Werror)
