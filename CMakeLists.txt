cmake_minimum_required(VERSION 3.16)
project(Simulations)

# sometimes we run the program on a server that has ancient standard compilers. In these cases, the following lines should be set to ensure the used std headers are up to date with the c++ 14 standard
#set(CMAKE_CXX_COMPILER /usr/local/bin/g++-11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# GPU (OpenACC) offload path. Requires the NVIDIA HPC SDK (nvc++); the CPU
# OpenMP build is unaffected and remains the default. Configure with e.g.
#   source /opt/nvidia/hpc_sdk/activate   # puts nvc++ on PATH + LD_LIBRARY_PATH
#   CXX=nvc++ cmake -DENABLE_OPENACC=ON -DCMAKE_BUILD_TYPE=Release ..
option(ENABLE_OPENACC "Build the GPU (OpenACC) offload path with nvc++" OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
  # NVIDIA HPC SDK. OpenMP is -mp (not -fopenmp); -Wall/-march/-funroll-loops
  # are GNU/Clang-only. -fast bundles the aggressive opt + relaxed-math flags.
  # The -mp/-acc flags must appear on the LINK line too so nvc pulls in the
  # OpenMP/OpenACC runtimes; CMake's FindOpenMP can't derive them for nvc, so
  # we wire them by hand and skip find_package(OpenMP) below.
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fast -mp -DNDEBUG")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -mp")
  # The HDF5 imported target exports /usr/include as a SYSTEM include, so CMake
  # would pass `-isystem /usr/include`. With nvc++ that reorders /usr/include
  # ahead of the C++ header dir and breaks libstdc++'s `#include_next <stdlib.h>`
  # (catastrophic build error). Demoting imported-target SYSTEM includes to plain
  # -I preserves the search order (verified: -I/usr/include compiles, -isystem
  # does not) and fixes the build.
  set(CMAKE_NO_SYSTEM_FROM_IMPORTED ON)
  if(ENABLE_OPENACC)
    # cc89 = RTX 4070 (Ada). -Minfo=accel prints the offload/kernel report.
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -acc=gpu -gpu=cc89 -Minfo=accel")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -acc=gpu -gpu=cc89")
    message(STATUS "OpenACC GPU offload ENABLED (nvc++, target cc89)")
  else()
    message(STATUS "nvc++ detected; building CPU/OpenMP path (OpenACC OFF)")
  endif()
else()
  # GCC / Clang: the established CPU/OpenMP flags.
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fopenmp -ffast-math -march=native -funroll-loops -DNDEBUG")
  if(ENABLE_OPENACC)
    message(FATAL_ERROR
      "ENABLE_OPENACC requires the NVIDIA HPC SDK. Re-run with a login shell "
      "that has nvc++ on PATH and -DCMAKE_CXX_COMPILER=nvc++ (or CXX=nvc++).")
  endif()
endif()
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
# nvc++ handles OpenMP via -mp (wired above); CMake's FindOpenMP does not know
# how to probe it, so we skip the module for NVHPC.
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
  find_package(OpenMP REQUIRED)
endif()

# HDF5 (C library) for the .h5 output. Module mode is portable and yields
# HDF5_INCLUDE_DIRS / HDF5_LIBRARIES. The vendored, header-only C++ wrapper
# lives in third_party/hdf5lib (included via <hdf5lib/...>).
find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/third_party ${HDF5_INCLUDE_DIRS})

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
