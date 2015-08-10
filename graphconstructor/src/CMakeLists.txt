project(graphconstructor CXX)
cmake_minimum_required(VERSION 2.8)

find_package(Boost 1.57 COMPONENTS thread filesystem regex REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")

set(CMAKE_BUILD_TYPE Release)
list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x -fopenmp")


add_executable(graphconstructor constructor.cpp concurrentbitvector.cpp compressedstring.cpp streamfastaparser.cpp test.cpp vertexenumerator.cpp)
include_directories(${Boost_INCLUDE_DIRS} ${graphconstructor_SOURCE_DIR} ~/include) 
link_directories(${Boost_LIBRARY_DIR} "~/lib" "~/lib/tbb")
target_link_libraries(graphconstructor ${Boost_LIBRARIES} "libtbb.so" "/gpfs/home/ium125/lib/libtpie.a")

set(CMAKE_PROJECT_NAME graphconstructor)
set(ROOT_DIR "${CMAKE_SOURCE_DIR}/../")

include(CPack)
