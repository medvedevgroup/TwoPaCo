project(twopaco CXX)

set(CMAKE_PROJECT_NAME twopaco)
set(CMAKE_CXX_STANDARD 14)

add_executable(twopaco ../common/dnachar.cpp constructor.cpp concurrentbitvector.cpp compressedstring.cpp ../common/streamfastaparser.cpp test.cpp vertexenumerator.cpp ../common/spooky/SpookyV2.cpp common.cpp)
link_directories(${TBB_LIB_DIR})
include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIR} "../common") 
target_link_libraries(twopaco  "tbb" "tbbmalloc" "tbbmalloc_proxy" "pthread" )
install(TARGETS twopaco RUNTIME DESTINATION bin)

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
