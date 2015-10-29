project(graphconstructor CXX)
cmake_minimum_required(VERSION 2.8)

find_package(Boost 1.54 COMPONENTS locale thread filesystem regex system REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")

set(CMAKE_BUILD_TYPE Release)
list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x -fopenmp")

add_executable(graphconstructor dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp streamfastaparser.cpp test.cpp vertexenumerator.cpp SpookyV2.cpp)
include_directories(${Boost_INCLUDE_DIRS} ${graphconstructor_SOURCE_DIR} ${TBB_INCLUDE_DIR}) 
target_link_libraries(graphconstructor ${Boost_LIBRARIES} "${TBB_LINK_DIR}/libtbb.so")

set(CMAKE_PROJECT_NAME graphconstructor)
set(ROOT_DIR "${CMAKE_SOURCE_DIR}/../")

include(CPack)
