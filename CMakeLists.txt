cmake_minimum_required(VERSION 2.8)
project (gnina C CXX)

#disable parts of caffe we don't need
set(USE_CUDNN 0)
set(USE_LMDB 0)
set(USE_LEVELDB 0)
set(USE_OPENCV 0)
set(BUILD_python 0)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

include(caffe/cmake/Utils.cmake)
include(caffe/cmake/Cuda.cmake)
#our version of caffe needs to pull in headers from gnina
include_directories(.) 
add_subdirectory(caffe)
include_directories(${CMAKE_BINARY_DIR}/caffe/include/)
include_directories(caffe/include)
add_subdirectory(gninasrc)
