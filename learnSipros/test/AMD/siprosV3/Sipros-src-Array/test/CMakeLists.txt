set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

add_executable(SiprosV3test ${SRC_LIST} ${PROJECT_SOURCE_DIR}/test/main.cpp)

target_link_libraries(SiprosV3test -fopenmp -static)