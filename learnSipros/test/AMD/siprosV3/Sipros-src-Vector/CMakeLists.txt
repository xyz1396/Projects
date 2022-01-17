project(siprosV3)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    message("Debug")
    set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "-fopenmp -std=gnu++11 -fpermissive -g -lprofiler -DGper")
else()
    message("Release")
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_CXX_FLAGS "-fopenmp -std=gnu++11 -fpermissive -O3")
endif()

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
include_directories(${PROJECT_SOURCE_DIR}/include)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

add_subdirectory(openmp)
add_subdirectory(mpi)
add_subdirectory(test)

