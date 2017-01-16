cmake_minimum_required(VERSION 2.8)
project (gnina C CXX)

#disable parts of caffe we don't need
set(USE_CUDNN 0)
set(USE_LMDB 0)
set(USE_LEVELDB 0)
set(USE_OPENCV 0)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

include(caffe/cmake/Utils.cmake)
include(caffe/cmake/Cuda.cmake)
#our version of caffe needs to pull in headers from gnina and we'll use c++11
include_directories(.) 
add_compile_options(-std=c++11)
list(APPEND CUDA_NVCC_FLAGS -std=c++11)

add_subdirectory(caffe)
include_directories(${CMAKE_BINARY_DIR}/caffe/include/)
include_directories(caffe/include)
add_subdirectory(gninasrc)
