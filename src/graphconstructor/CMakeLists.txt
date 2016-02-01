project(graphconstructor CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME graphconstructor)
find_package(Boost 1.57 COMPONENTS locale thread filesystem regex system REQUIRED)
if(WIN32)
	set(extlibrary "${TBB_LIBRARY_DIR}/tbb.lib")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x -fopenmp")
	set(extlibrary "libtbb")
endif()

add_executable(graphconstructor ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp)
include_directories(${Boost_INCLUDE_DIR} ${graphconstructor_SOURCE_DIR} ${TBB_INCLUDE_DIR} "../common") 
target_link_libraries(graphconstructor ${Boost_LIBRARIES} ${extlibrary})



include(CPack)
