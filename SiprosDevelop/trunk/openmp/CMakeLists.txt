set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

add_executable(SiprosV3omp ${SRC_LIST} ${PROJECT_SOURCE_DIR}/openmp/main.cpp)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_link_libraries(SiprosV3omp PUBLIC OpenMP::OpenMP_CXX profiler)
    target_compile_options(SiprosV3omp PUBLIC -ffast-math -march=native)
    target_compile_definitions(SiprosV3omp PUBLIC Gper)
else()
    target_link_libraries(SiprosV3omp PUBLIC OpenMP::OpenMP_CXX -static)
    target_compile_options(SiprosV3omp PUBLIC -ffast-math -march=native)
endif()
