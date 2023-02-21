cmake_minimum_required(VERSION 3.12)
project(Raxport LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../bin")

add_executable(Raxport raxportOMP.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(Raxport PUBLIC OpenMP::OpenMP_CXX)
