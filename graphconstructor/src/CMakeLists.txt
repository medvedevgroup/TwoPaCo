cmake_minimum_required(VERSION 2.8)
project(graphconstructor CXX)
set(CMAKE_BUILD_TYPE Release)
list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x -fopenmp")

add_executable(graphconstructor constructor.cpp dnastring.cpp streamfastaparser.cpp test.cpp vertexenumerator.cpp lib/SpookyV2.cpp)
set(CMAKE_PROJECT_NAME graphconstructor)
set(ROOT_DIR "${CMAKE_SOURCE_DIR}/../")

include(CPack)
