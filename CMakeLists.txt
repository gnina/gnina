cmake_minimum_required(VERSION 3.12)
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

if(NOT ${CUDA_cublas_device_LIBRARY}) #workaround for oldercmake/newer cuda issue
 set(CUDA_cublas_device_LIBRARY "")
endif()

find_package(libmolgrid)
if(NOT LIBMOLGRID_LIBRARY)
 message("libmolgrid will be fetched from git")
 include(ExternalProject)
 find_package(Git REQUIRED)
 set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)
 ExternalProject_Add(libmolgrid
    GIT_REPOSITORY https://github.com/gnina/libmolgrid 
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
  )
  set(LIBMOLGRID_INCLUDE_DIR ${EXTERNAL_INSTALL_LOCATION}/include)
  set(LIBMOLGRID_LIBRARY ${EXTERNAL_INSTALL_LOCATION}/lib/libmolgrid.so)
else()
  add_custom_target(libmolgrid)	
endif()

include_directories(${LIBMOLGRID_INCLUDE_DIR})

find_package(OpenBabel3 REQUIRED)
#get_target_property(OPENBABEL_LIBRARIES openbabel LOCATION)
set(OPENBABEL_LIBRARIES $<TARGET_FILE:openbabel>)
set(OPENBABEL_INCLUDE_DIR ${OpenBabel3_INCLUDE_DIRS})
include_directories(${OPENBABEL_INCLUDE_DIR})

include(caffe/cmake/Utils.cmake)
include(caffe/cmake/Cuda.cmake)

#our version of caffe needs to pull in headers from gnina and we'll use c++11
include_directories(.) 

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)


add_subdirectory(caffe)
add_dependencies(caffe libmolgrid)

include_directories(${CMAKE_BINARY_DIR}/caffe/include/)
include_directories(caffe/include)
add_subdirectory(gninasrc)

# enable testing
include(CTest)
enable_testing()
# define tests

add_subdirectory(test)
