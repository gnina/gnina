cmake_minimum_required(VERSION 2.8)
project (gnina C CXX)

set(USE_CUDNN 0)
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

include(caffe/cmake/Utils.cmake)
include(caffe/cmake/Cuda.cmake)
add_subdirectory(caffe)
add_subdirectory(src)
