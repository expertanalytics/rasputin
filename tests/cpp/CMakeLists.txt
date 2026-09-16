include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.6.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

add_executable(test_point test_point.cpp)
target_link_libraries(test_point PRIVATE
    terrain_headers
    Catch2::Catch2WithMain
)
target_compile_options(test_point PRIVATE
    -Wall -Wextra -Wpedantic -Werror
)

catch_discover_tests(test_point)
