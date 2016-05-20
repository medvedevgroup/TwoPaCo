project(twopaco CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME twopaco)
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x")
endif()

<<<<<<< HEAD
add_executable(twopaco ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp common.cpp)
link_directories(${TBB_LIB_DIR})
include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIR} "../common") 
target_link_libraries(twopaco  "tbb")
=======
add_executable(graphconstructor ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp common.cpp)
include_directories(${graphconstructor_SOURCE_DIR} "../common") 
target_link_libraries(graphconstructor  "tbb" "cmph")
>>>>>>> 9a39ef882f3bbabab346bad34ac73f2e5efc0173

