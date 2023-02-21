AUX_SOURCE_DIRECTORY("${PROJECT_SOURCE_DIR}/test" DIR_TESTS)
AUX_SOURCE_DIRECTORY("${PROJECT_SOURCE_DIR}/src" SRC)

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include" 
    "${PROJECT_SOURCE_DIR}/3rdParty/MSToolkit/include")

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin/test")

add_executable(SiprosTest
  ${SRC}
  ${DIR_TESTS}
)

target_link_libraries(
  SiprosTest
  MS
  -fopenmp
  -static
)
