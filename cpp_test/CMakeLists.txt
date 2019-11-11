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
    set(catchlib Catch2)

# Otherwise look for catch2 installed on system
else()
    find_package(Catch2 REQUIRED)
    set(catchlib Catch2::Catch2)
    include(Catch)
    include(CTest)
endif()

# Build tests
# -----------
add_executable(rasputin_test test_sun_position.cpp)
target_link_libraries(rasputin_test rasputin ${catchlib} ${RASPUTIN_DEPENDENCIES})

catch_discover_tests(rasputin_test)
