aux_source_directory(src PRatioSrc)
include_directories(include ${PROJECT_SOURCE_DIR}/Common/include ${PROJECT_SOURCE_DIR}/SicForma/include)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")
add_executable(ProRata ${PRatioSrc} DriverProteomeInfo.cpp)
target_link_libraries(ProRata PUBLIC Common SicForma -static)