cmake_minimum_required(VERSION 3.1)



set(CMAKE_C_COMP "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx")
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("${CMAKE_BUILD_TYPE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=gnu++11 -DGCC -D_FILE_OFFSET_BITS=64 -g -fpermissive")
    add_compile_options(-fopenmp)
    add_compile_options(-static)
else()
    message("${CMAKE_BUILD_TYPE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=gnu++11 -DGCC -D_FILE_OFFSET_BITS=64 -fpermissive -static -O3")
    add_compile_options(-fopenmp)
    add_compile_options(-static)
endif()

# # add_compile_options(-g) 
# add_compile_options(-std=gnu++11)
# add_compile_options(-DGCC)
# add_compile_options(-D_FILE_OFFSET_BITS=64)
# add_compile_options(-fpermissive)
# add_compile_options(-fopenmp)
# # add_compile_options(-O3)
# add_compile_options(-static)
# # add_compile_options(-Wall -Wextra)

project(Sipros)
 
# 3rdParty library
# Build with original makefile
set(MS_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/MSToolkit/libmstoolkitlite.a)
 
add_custom_target(build_MS ALL
        COMMAND ${CMAKE_MAKE_PROGRAM} lite
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/MSToolkit
        COMMENT "Original mstoolkit makefile target")
 
add_library(MS STATIC IMPORTED)
 
set_property(TARGET MS APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(MS PROPERTIES IMPORTED_LOCATION_NOCONFIG "${MS_LIBRARY}")
 
add_dependencies(MS build_MS)
 
# source file
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src SRC_LIST)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/MSToolkit/include)

add_executable(Sipros ${SRC_LIST})

target_link_libraries(Sipros MS)
target_link_libraries(Sipros -fopenmp -static)