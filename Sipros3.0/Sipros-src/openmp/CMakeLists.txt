set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

add_executable(SiprosV3omp ${SRC_LIST} ${PROJECT_SOURCE_DIR}/openmp/main.cpp)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_link_libraries(SiprosV3omp -fopenmp -lprofiler)
else()
    target_link_libraries(SiprosV3omp -fopenmp -static)
endif()

