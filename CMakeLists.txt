cmake_minimum_required(VERSION 2.6)
project(Simulations)

# sometimes we run the program on a server that has ancient standard compilers. In these cases, the following lines should be set to ensure the used std headers are up to date with the c++ 14 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CXX_COMPILER_LIBS "" CACHE FILEPATH "Path to cpp compiler libs")
link_directories(${CXX_COMPILER_LIBS})

set(CXX_STANDARD_INCLUDE_DIRECTORIES "" CACHE FILEPATH "Path to cpp compiler includes")
SET(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES "${CXX_STANDARD_INCLUDE_DIRECTORIES}")

# copy this to the project dir whenever a change to file tree structure has been made.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# this puts the all-purpose GPU stuff in top dir and creates a compiler macro to get the dir
add_compile_definitions(OPENCL_KERNELS_DIR="empty")

include_directories(core/include)

# the subprojects. core is the most important, as it provides the library used for all the other projects that contains basic tasks always needed
add_subdirectory(warm_up)
add_subdirectory(Quark_Photon_Vertex)
