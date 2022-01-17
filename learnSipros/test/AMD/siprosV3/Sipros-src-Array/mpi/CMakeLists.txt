set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(SiprosV3mpi ${SRC_LIST} ${PROJECT_SOURCE_DIR}/mpi/mpimain.cpp)
