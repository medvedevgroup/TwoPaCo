project(twopaco CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME twopaco)
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x")
endif()

add_executable(twopaco ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp common.cpp)
link_directories(${TBB_LIB_DIR})
include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIR} "../common") 
target_link_libraries(twopaco  "tbb")

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "3")
