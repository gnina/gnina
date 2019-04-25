cmake_minimum_required(VERSION 3.13)
project (gnina C CXX CUDA)

#disable parts of caffe we don't need
set(USE_LMDB 0)
set(USE_LEVELDB 0)
set(USE_OPENCV 0)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA 9.0 REQUIRED)
find_package(libmolgrid REQUIRED)
include_directories($LIBMOLGRID_INCLUDE_DIR)


include(caffe/cmake/Utils.cmake)
include(caffe/cmake/Cuda.cmake)

#our version of caffe needs to pull in headers from gnina and we'll use c++11
include_directories(.) 

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)


add_subdirectory(caffe)

include_directories(${CMAKE_BINARY_DIR}/caffe/include/)
include_directories(caffe/include)
add_subdirectory(gninasrc)

# enable testing
include(CTest)
enable_testing()
# define tests

add_subdirectory(test)