project(graphconstructor CXX)
cmake_minimum_required(VERSION 2.8)

find_package(Boost 1.57 COMPONENTS thread REQUIRED)
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")

set(CMAKE_BUILD_TYPE Release)
list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x")

include_directories(${graphconstructor_SOURCE_DIR}/tclap)
add_executable(graphconstructor constructor.cpp bloomfilter.cpp dnastring.cpp streamfastaparser.cpp test.cpp vertexenumerator.cpp lib/SpookyV2.cpp)
include_directories(${Boost_INCLUDE_DIRS}) 
link_directories(${Boost_LIBRARY_DIR})
target_link_libraries(graphconstructor ${Boost_LIBRARIES})

set(CMAKE_PROJECT_NAME graphconstructor)
set(ROOT_DIR "${CMAKE_SOURCE_DIR}/../")

include(CPack)
