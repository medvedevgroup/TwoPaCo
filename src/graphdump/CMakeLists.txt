project(graphdump CXX)

set(CMAKE_PROJECT_NAME graphdump)
set(CMAKE_CXX_STANDARD 14)

add_executable(graphdump graphdump.cpp ../common/dnachar.cpp ../common/streamfastaparser.cpp)
link_directories(${TBB_LIB_DIR})
include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIR} "../common")
target_link_libraries(graphdump  "tbb" "pthread")
install(TARGETS graphdump RUNTIME DESTINATION bin)

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
