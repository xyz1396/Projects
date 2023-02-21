cmake_minimum_required(VERSION 3.1)
project(Sipros)


set(CMAKE_C_COMP "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("${CMAKE_BUILD_TYPE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -DGCC -D_FILE_OFFSET_BITS=64 -g2 -fpermissive -O0")
    add_compile_options(-fopenmp)
    add_compile_options(-static)
else()
    message("${CMAKE_BUILD_TYPE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -DGCC -D_FILE_OFFSET_BITS=64 -fpermissive -static -O3")
    add_compile_options(-fopenmp)
    add_compile_options(-static)
endif()
 
# 3rdParty library
# Build with original makefile
set(MS_LIBRARY ${PROJECT_SOURCE_DIR}/3rdParty/MSToolkit/libmstoolkitlite.a)
 
add_custom_target(build_MS ALL
        COMMAND ${CMAKE_MAKE_PROGRAM} lite -j6
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/3rdParty/MSToolkit
        COMMENT "Original mstoolkit makefile target")
 
add_library(MS STATIC IMPORTED)
 
set_property(TARGET MS APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(MS PROPERTIES IMPORTED_LOCATION_NOCONFIG "${MS_LIBRARY}")
 
add_dependencies(MS build_MS)

# enable_testing()
add_subdirectory(test)
 
# source file
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/3rdParty/MSToolkit/include)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")
add_executable(Sipros ${SRC_LIST} ${PROJECT_SOURCE_DIR}/main.cpp)

target_link_libraries(Sipros MS)
target_link_libraries(Sipros -fopenmp -static)