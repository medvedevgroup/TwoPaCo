cmake_minimum_required(VERSION 3.15)
project(root)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE RELEASE)
add_subdirectory(graphdump)
add_subdirectory(graphconstructor)
include(CPack)
set(THREADS_PREFER_PTHREAD_FLAG ON)
