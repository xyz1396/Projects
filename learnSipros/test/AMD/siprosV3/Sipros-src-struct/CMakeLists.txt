project(siprosV3)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message("${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 11)
find_package(OpenMP REQUIRED)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
include_directories(${PROJECT_SOURCE_DIR}/include)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

add_subdirectory(openmp)
add_subdirectory(mpi)
add_subdirectory(test)

