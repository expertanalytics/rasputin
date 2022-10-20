enable_testing()

# Find Catch2
# -----------
find_package(Catch2 3 REQUIRED)

# Build tests
# -----------
add_executable(test_sun_position test_sun_position.cpp)
target_link_libraries(test_sun_position PRIVATE rasputin Catch2::Catch2WithMain ${RASPUTIN_DEPENDENCIES})
add_executable(test_triangulate test_triangulate.cpp)
target_link_libraries(test_triangulate PRIVATE rasputin Catch2::Catch2WithMain ${RASPUTIN_DEPENDENCIES})
