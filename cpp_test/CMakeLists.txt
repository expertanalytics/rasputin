enable_testing()

# Find Catch2
# -----------

# First look under lib/Catch2
if (EXISTS "${RASPUTIN_LIB_DIR}/catch2")
    include("${RASPUTIN_LIB_DIR}/catch2/contrib/Catch.cmake")
    include(CTest)
    add_library(Catch2 INTERFACE)
    target_include_directories(Catch2 INTERFACE
        "${RASPUTIN_LIB_DIR}/catch2/single_include")

# Otherwise look for catch2 installed on system
else()
    find_package(Catch2 REQUIRED)
    # Import into scope
    add_library(Catch2 INTERFACE IMPORTED)
    include(Catch)
    include(CTest)
endif()

# Build tests
# -----------
add_executable(foo test_sun_position.cpp)
target_link_libraries(foo rasputin ${Catch2} ${RASPUTIN_DEPENDENCIES})

catch_discover_tests(foo)
