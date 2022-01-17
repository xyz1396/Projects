project(siprosV3)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_CXX_FLAGS "-fopenmp -std=gnu++11 -fpermissive -O3 -lprofiler")
set(CMAKE_CXX_FLAGS "-fopenmp -std=gnu++11 -fpermissive -O3")

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
include_directories(${PROJECT_SOURCE_DIR}/include)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

add_subdirectory(openmp)
add_subdirectory(mpi)

