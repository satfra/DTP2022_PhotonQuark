project(Quark_Photon_Vertex)

# add the files of our current project
include_directories(./include)

# build the current executable
add_executable(Quark_Photon_Vertex Simulation.cpp ${SRCFILES})

# find pthread library
find_package(Threads REQUIRED)

# link against the external libraries
target_link_libraries(Quark_Photon_Vertex Threads::Threads)
