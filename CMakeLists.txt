cmake_minimum_required(VERSION 2.6)
project(Simulations)

# sometimes we run the program on a server that has ancient standard compilers. In these cases, the following lines should be set to ensure the used std headers are up to date with the c++ 14 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")

set(CXX_COMPILER_LIBS "" CACHE FILEPATH "Path to cpp compiler libs")
link_directories(${CXX_COMPILER_LIBS})
link_directories(/usr/local/lib)

# copy this to the project dir whenever a change to file tree structure has been made.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(core/include)
include_directories(/usr/local/include)

# the subprojects. core is the most important, as it provides the library used for all the other projects that contains basic tasks always needed
add_subdirectory(Quark_Photon_Vertex)
