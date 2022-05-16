project(Quark_Photon_Vertex)

# add the files of our current project
include_directories(./include)
FILE(GLOB SRCFILES ./source/*.cpp)

# build the current executable
add_executable(Quark_Photon_Vertex Simulation.cpp ${SRCFILES})

# link against the external libraries
target_link_libraries(Quark_Photon_Vertex pthread)

# lastly, we copy the config and the openCL kernels so that the program can be used.
add_compile_definitions(KPM_OPENCL_KERNELS_DIR="empty")
