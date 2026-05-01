cmake_minimum_required(VERSION 3.16)
project(Simulations)

# sometimes we run the program on a server that has ancient standard compilers. In these cases, the following lines should be set to ensure the used std headers are up to date with the c++ 14 standard
#set(CMAKE_CXX_COMPILER /usr/local/bin/g++-11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fopenmp -ffast-math -march=native -funroll-loops -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")

# Enable link-time / interprocedural optimisation for Release if the toolchain
# supports it. Lets the linker inline across translation units, which matters
# here because most hot code lives in headers anyway but the few .cpp boundaries
# (Simulation.cpp ↔ tests) still benefit.
include(CheckIPOSupported)
check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_err)
if(_ipo_ok)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  message(STATUS "IPO/LTO enabled for Release builds")
else()
  message(STATUS "IPO/LTO not supported: ${_ipo_err}")
endif()

# Locate OpenMP centrally so subprojects can link OpenMP::OpenMP_CXX directly.
find_package(OpenMP REQUIRED)

# Reduced grid sizes for performance benchmarking. See bench/run_bench.sh.
option(BENCHMARK_PROFILE "Compile with reduced grid sizes for benchmarking" OFF)
if(BENCHMARK_PROFILE)
  add_compile_definitions(BENCHMARK_PROFILE)
  message(STATUS "BENCHMARK_PROFILE enabled — using parameters_bench.hh grid sizes")
endif()

# copy this to the project dir whenever a change to file tree structure has been made.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core/include)

# the subprojects. core is the most important, as it provides the library used for all the other projects that contains basic tasks always needed
add_subdirectory(Quark_Photon_Vertex)

enable_testing()
add_subdirectory(tests)
