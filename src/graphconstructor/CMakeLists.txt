project(graphconstructor CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME graphconstructor)
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x")
endif()

add_executable(graphconstructor ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp)
include_directories(${graphconstructor_SOURCE_DIR} "../common") 

include(CPack)
