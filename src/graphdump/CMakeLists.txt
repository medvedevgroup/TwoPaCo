project(graphdump CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME graphdump)
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x")
endif()

add_executable(graphdump graphdump.cpp ../common/dnachar.cpp ../common/streamfastaparser.cpp)
link_directories(${TBB_LIB_DIR})
include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIR} "../common")
target_link_libraries(graphdump  "tbb")


set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "3")
