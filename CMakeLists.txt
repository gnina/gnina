cmake_minimum_required(VERSION 3.12)
project (gnina C CXX CUDA)

#disable parts of caffe we don't need
set(USE_LMDB 0)
set(USE_LEVELDB 0)
set(USE_OPENCV 0)

if(POLICY CMP0104)
  cmake_policy(SET CMP0104 OLD) #TODO: clean up Cuda.cmake to be modern to avoid this workaround
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
find_package(CUDA 11.0 REQUIRED)

if(NOT ${CUDA_cublas_device_LIBRARY}) #workaround for oldercmake/newer cuda issue
 set(CUDA_cublas_device_LIBRARY "")
endif()

execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

# Convert command output into a CMake list
STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

function(print_properties)
    message ("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
endfunction(print_properties)

function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
      message("There is no target named '${tgt}'")
      return()
    endif()

    foreach (prop ${CMAKE_PROPERTY_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
    # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
    if(prop STREQUAL "LOCATION" OR prop MATCHES "^LOCATION_" OR prop MATCHES "_LOCATION$")
        continue()
    endif()
        # message ("Checking ${prop}")
        get_property(propval TARGET ${tgt} PROPERTY ${prop} SET)
        if (propval)
            get_target_property(propval ${tgt} ${prop})
            message ("${tgt} ${prop} = ${propval}")
        endif()
    endforeach(prop)
endfunction(print_target_properties)

find_package(OpenBabel3 REQUIRED)
print_target_properties(openbabel)
set(OpenBabel3_LIBRARIES $<TARGET_FILE:openbabel>)
include_directories($<TARGET_PROPERTY:openbabel,INTERFACE_INCLUDE_DIRECTORIES>)
include_directories(${OpenBabel3_INCLUDE_DIRS})
  
find_package(ZLIB REQUIRED)
  
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

include(caffe/cmake/Utils.cmake)
include(caffe/cmake/Cuda.cmake)

#our version of caffe needs to pull in headers from gnina and we'll use c++11
include_directories(.)

add_subdirectory(caffe)
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
