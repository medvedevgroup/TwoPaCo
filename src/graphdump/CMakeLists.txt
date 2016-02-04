project(graphdump CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME graphdump)
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x")
endif()

add_executable(graphdump graphdump.cpp)
include_directories(${graphdump_SOURCE_DIR} "../common") 

