project(graphdump CXX)

set(CMAKE_PROJECT_NAME graphdump)
set(CMAKE_CXX_STANDARD 14)
add_executable(graphdump graphdump.cpp ../common/dnachar.cpp ../common/streamfastaparser.cpp)
include_directories(${twopaco_SOURCE_DIR} "../common")
target_link_libraries(graphdump "pthread")
install(TARGETS graphdump RUNTIME DESTINATION bin)

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
