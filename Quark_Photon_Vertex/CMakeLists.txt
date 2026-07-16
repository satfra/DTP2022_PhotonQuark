project(Quark_Photon_Vertex)

# add the files of our current project
include_directories(./include)

# build the current executable
add_executable(Quark_Photon_Vertex Simulation.cpp ${SRCFILES})

# find pthread library
find_package(Threads REQUIRED)

# link against the external libraries. For GCC/Clang, OpenMP::OpenMP_CXX (found
# by the parent CMakeLists) puts -fopenmp on the link line. For nvc++ the
# OpenMP/OpenACC runtimes come from the -mp/-acc link flags set by the parent,
# and no OpenMP CMake target exists — so link it only when present.
target_link_libraries(Quark_Photon_Vertex Threads::Threads ${HDF5_LIBRARIES})
if(TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(Quark_Photon_Vertex OpenMP::OpenMP_CXX)
endif()
