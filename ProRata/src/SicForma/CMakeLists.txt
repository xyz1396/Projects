aux_source_directory(src SicFormaSrc)
include_directories(include ${PROJECT_SOURCE_DIR}/Common/include ${PROJECT_SOURCE_DIR}/PRatio/include)
add_library(SicForma ${SicFormaSrc} ${PROJECT_SOURCE_DIR}/PRatio/src/peptideInfo.cpp)
target_link_libraries(SicForma PUBLIC Common)