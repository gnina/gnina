cmake_minimum_required(VERSION 3.12)
project (gnina C CXX CUDA)


if(POLICY CMP0104)
  cmake_policy(SET CMP0104 OLD) #TODO: clean up Cuda.cmake to be modern to avoid this workaround
endif()
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON) # for vs code

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_COVERAGE "Build with code coverage" OFF)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -Wno-unknown-pragmas")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (BUILD_COVERAGE)
  set(CMAKE_CXX_FLAGS_RELEASE "-fprofile-arcs -ftest-coverage")
endif()

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA 12.0 REQUIRED)
set(USE_CUDNN 1)


#zlib
find_package(ZLIB REQUIRED)

#openbabel
find_package(OpenBabel3 REQUIRED)
set(OpenBabel3_LIBRARIES $<TARGET_FILE:openbabel>)
include_directories($<TARGET_PROPERTY:openbabel,INTERFACE_INCLUDE_DIRECTORIES>)
include_directories(${OpenBabel3_INCLUDE_DIRS})
  

#torch
include(FetchContent)
# https://discuss.pytorch.org/t/using-fetchcontent-declare-in-cmake-to-pull-libtorch/177449
FetchContent_Declare(
        libtorch
        URL "https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-2.2.0%2Bcu121.zip"
)

set(TORCH_CUDA_ARCH_LIST "5.0 5.2 6.0 7.0 8.0 8.6 8.9 9.0") #All includes 3.5 which is problematic
FetchContent_MakeAvailable(libtorch)
list(APPEND CMAKE_PREFIX_PATH "${libtorch_SOURCE_DIR}")
find_package(Torch REQUIRED)

#libmolgrid  
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
  set(LIBMOLGRID_LIBRARY ${EXTERNAL_INSTALL_LOCATION}/lib/libmolgrid.a)
  include_directories(${LIBMOLGRID_INCLUDE_DIR})
  find_package (Python COMPONENTS Interpreter Development)
  set(PYTHONPATH $ENV{PYTHONPATH})
  set(ENV{PYTHONPATH} "${PYTHONPATH}:${EXTERNAL_INSTALL_LOCATION}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/")
else()
  add_custom_target(libmolgrid)
endif()

include_directories(${LIBMOLGRID_INCLUDE_DIR})


#caffe
#disable parts of caffe we don't need
set(USE_LMDB 0)
set(USE_LEVELDB 0)
set(USE_OPENCV 0)
set(USE_CUDNN 0) #problematic and conflicts with torch version

set(CUDA_ARCH_NAME "All")
include(caffe/cmake/Utils.cmake)
include(caffe/cmake/Cuda.cmake)
#our version of caffe needs to pull in headers from gnina and we'll use c++11
include_directories(.)
add_subdirectory(caffe)

unset(CUDA_ARCH_NAME CACHE) #remove conflict with torch

add_dependencies(caffe libmolgrid)

include_directories(${CMAKE_BINARY_DIR}/caffe/include/)
include_directories(caffe/include)
add_subdirectory(gninasrc)

include_directories(${LIBMOLGRID_INCLUDE_DIR})


# enable testing
include(CTest)
enable_testing()
# define tests

add_subdirectory(test)
