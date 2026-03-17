cmake_minimum_required(VERSION 3.14)

include(FetchContent)

# Download doctest header-only library. We use FetchContent_Populate (not
# MakeAvailable) to avoid running doctest's own CMakeLists.txt, which declares
# a cmake_minimum_required version that modern CMake no longer accepts.
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG        v2.4.11
  GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(doctest)
if(NOT doctest_POPULATED)
  FetchContent_Populate(doctest)
endif()

# doctest_SOURCE_DIR is set by FetchContent_Populate; fall back to the
# conventional _deps path in case variable scoping is tricky.
if(NOT doctest_SOURCE_DIR)
  set(doctest_SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/doctest-src)
endif()

# Include the parent of doctest/doctest.h so that #include <doctest/doctest.h> works
set(TEST_INCLUDES
  ${doctest_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/core/include
  ${CMAKE_SOURCE_DIR}/Quark_Photon_Vertex/include
)

add_executable(test_core test_core.cpp)
target_include_directories(test_core PRIVATE ${TEST_INCLUDES})
target_compile_features(test_core PRIVATE cxx_std_20)
target_compile_options(test_core PRIVATE -Wall)
add_test(NAME test_core COMMAND test_core)

add_executable(test_physics test_physics.cpp)
target_include_directories(test_physics PRIVATE ${TEST_INCLUDES})
target_compile_features(test_physics PRIVATE cxx_std_20)
target_compile_options(test_physics PRIVATE -Wall)
add_test(NAME test_physics COMMAND test_physics)
