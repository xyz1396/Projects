set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx")

add_executable(SiprosV3mpi ${SRC_LIST} ${PROJECT_SOURCE_DIR}/mpi/mpimain.cpp)
