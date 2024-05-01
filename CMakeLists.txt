cmake_minimum_required(VERSION 3.25)
project (gnina C CXX CUDA)


if(POLICY CMP0104)
  cmake_policy(SET CMP0104 OLD) #TODO: clean up Cuda.cmake to be modern to avoid this workaround
endif()
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

#RPATH instead of RUNPATH, needed for transitive torch libraries
#TODO: get static build working so we don't need this
SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags")

add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON) # for vs code

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_COVERAGE "Build with code coverage" OFF)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_INSTALL_PREFIX}/lib/gnina")
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -Wno-unknown-pragmas")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (BUILD_COVERAGE)
  set(CMAKE_CXX_FLAGS_RELEASE "-fprofile-arcs -ftest-coverage")
endif()

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA 12.0 REQUIRED)
set(USE_CUDNN 1)

function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if ((NOT DEFINED ARGV0) OR _variableName MATCHES ${ARGV0})
            message(STATUS "${_variableName}=${${_variableName}}")
        endif()
    endforeach()
endfunction()

#zlib
find_package(ZLIB REQUIRED)

#openbabel
find_package(OpenBabel3 REQUIRED)
set(OpenBabel3_LIBRARIES $<TARGET_FILE:openbabel>)
include_directories($<TARGET_PROPERTY:openbabel,INTERFACE_INCLUDE_DIRECTORIES>)
include_directories(${OpenBabel3_INCLUDE_DIRS})

#make torch use system protobuf to be compat with caffe
find_package(Protobuf REQUIRED)

#torch
include(FetchContent)
# https://discuss.pytorch.org/t/using-fetchcontent-declare-in-cmake-to-pull-libtorch/177449
FetchContent_Declare(
        libtorch
        URL "https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-2.2.2%2Bcu121.zip"
)

#set(TORCH_CUDA_ARCH_LIST "5.0 5.2 6.0 7.0 8.0 8.6 8.9 9.0") #All includes 3.5 which is problematic
set(CMAKE_CUDA_ARCHITECTURES "all-major")
FetchContent_MakeAvailable(libtorch)
list(APPEND CMAKE_PREFIX_PATH "${libtorch_SOURCE_DIR}")
find_package(Torch REQUIRED)

cmake_path(GET TORCH_LIBRARY PARENT_PATH TORCH_LIBDIR)
set(TORCH_NEEDED
libgomp-98b21ff3.so.1 
libcudart-9335f6a2.so.12 
libnvToolsExt-847d78f2.so.1 
libcublas-37d11411.so.12 
libcublasLt-f97bfc2c.so.12 
libtorch_cpu.so
libtorch_cuda.so
libc10_cuda.so
libc10.so
libtorch.so
libcudnn.so.8
libcudnn_adv_infer.so.8
libcudnn_adv_train.so.8
libcudnn_cnn_infer.so.8
libcudnn_cnn_train.so.8
libcudnn_ops_infer.so.8
libcudnn_ops_train.so.8
libjitbackend_test.so
libnnapi_backend.so
libnvrtc-b51b459d.so.12
libnvrtc-builtins-6c5639ce.so.12.1
libshm.so
libtorchbind_test.so
libtorch_cuda_linalg.so
libtorch_global_deps.so
libtorch_python.so
)
list(TRANSFORM TORCH_NEEDED PREPEND "${TORCH_LIBDIR}/")

#install in subdirectory to avoid messing up a system torch install
install(FILES ${TORCH_NEEDED} DESTINATION lib/gnina)

# eventually build static pytorch libraries

# FetchContent_Declare(
#         pytorch
#         GIT_REPOSITORY "https://github.com/pytorch/pytorch"
#         GIT_TAG "v2.2.2"
# )
# Configure for minimal, static, openblas-enabled CPU library.
# set(BUILD_SHARED_LIBS OFF)
# set(BUILD_CAFFE2_OPS OFF)
# set(BUILD_CAFFE2 OFF)
# set(USE_MKL 0)

# set(CMAKE_CUDA_ARCHITECTURES "all-major")
# FetchContent_MakeAvailable(pytorch)
# list(APPEND CMAKE_PREFIX_PATH "${pytorch_BINARY_DIR}")

message(STATUS ${CMAKE_PREFIX_PATH})
#

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
include_directories(.)

add_subdirectory(gninasrc)


# enable testing
include(CTest)
enable_testing()
# define tests

add_subdirectory(test)
