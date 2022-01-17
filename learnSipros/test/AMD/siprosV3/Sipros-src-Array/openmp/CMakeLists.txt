set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

add_executable(SiprosV3omp ${SRC_LIST} ${PROJECT_SOURCE_DIR}/openmp/main.cpp)

# target_link_libraries(SiprosV3omp -fopenmp -lprofiler)
target_link_libraries(SiprosV3omp -fopenmp)
