enable_testing()
# CMake configuration for tests

# Find Catch2
find_package(Catch2 3 REQUIRED)

# Add test sources
file(GLOB TEST_SOURCES "*.cpp")
list(REMOVE_ITEM TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bench_seqan2_vs_igor.cpp")

add_executable(igor_tests ${TEST_SOURCES})

# Link against the core library and Catch2
# (Assuming igor::Core is the main library target)
target_link_libraries(igor_tests PRIVATE igor::Core Catch2::Catch2WithMain)

# Provide the project source directory so tests can locate data files
target_compile_definitions(igor_tests PRIVATE IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# Set include directories if needed
# target_include_directories(igor_tests PRIVATE ../src/igor/Core)

# Register the tests with CTest
include(Catch)
catch_discover_tests(igor_tests
    PROPERTIES TIMEOUT 3600  # generation test with 1M sequences needs headroom
)

add_executable(bench_seqan2_vs_igor bench_seqan2_vs_igor.cpp)
target_link_libraries(bench_seqan2_vs_igor PRIVATE igor::Core)
target_compile_definitions(bench_seqan2_vs_igor PRIVATE IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
