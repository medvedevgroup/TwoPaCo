project(graphconstructor CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME graphconstructor)
if(WIN32)
	set(extlibrary "${TBB_LIBRARY_DIR}/tbb.lib")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x -fopenmp")
	set(extlibrary "libtbb")
endif()

add_executable(graphconstructor ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp)
include_directories(${graphconstructor_SOURCE_DIR} ${TBB_INCLUDE_DIR} "../common") 
target_link_libraries(graphconstructor ${extlibrary})

include(CPack)
