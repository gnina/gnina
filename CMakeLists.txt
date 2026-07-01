cmake_minimum_required(VERSION 3.25)

# ─────────────────────────────────────────────────────────────────────────────
# Backend auto-detection
# On Apple platforms: Metal backend (M1/M2/M3, no CUDA).
# On Linux/Windows:   CUDA backend (original behaviour).
# ─────────────────────────────────────────────────────────────────────────────
if(APPLE)
  option(USE_METAL "Build with Metal GPU backend (Apple Silicon)" ON)
else()
  option(USE_METAL "Build with Metal GPU backend (Apple Silicon)" OFF)
endif()

if(USE_METAL)
  message(STATUS "Backend: Metal (Apple Silicon) — CUDA disabled")
  project(gnina C CXX OBJCXX)   # OBJCXX for metal_context.mm
else()
  message(STATUS "Backend: CUDA (NVIDIA)")
  project(gnina C CXX CUDA)
endif()

cmake_policy(SET CMP0104 NEW)
cmake_policy(SET CMP0135 NEW)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
cmake_policy(SET CMP0094 NEW)   # location-first python finding
cmake_policy(SET CMP0146 OLD)   # too baked in

set(CMAKE_CXX_STANDARD 17)
if(NOT USE_METAL)
  set(CMAKE_CUDA_STANDARD 17)
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Global compile definitions
# ─────────────────────────────────────────────────────────────────────────────
add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
if(USE_METAL)
  add_definitions(-DUSE_METAL)
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_COVERAGE "Build with code coverage" OFF)
option(GNINA_FORCE_EXTERNAL_LIBS "Use external libraries only" OFF)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_INSTALL_PREFIX}/lib/gnina")

if(USE_METAL)
  set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -Wno-unknown-pragmas -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION")
  set(CMAKE_OBJCXX_FLAGS "-Wno-deprecated-declarations")
else()
  set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -Wno-unknown-pragmas")
  # RPATH trick needed for transitive torch libraries on Linux
  SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(BUILD_COVERAGE)
  set(CMAKE_CXX_FLAGS_RELEASE "-fprofile-arcs -ftest-coverage")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# CUDA (Linux/Windows only)
# ─────────────────────────────────────────────────────────────────────────────
if(NOT USE_METAL)
  set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
  find_package(CUDA 12.0 REQUIRED)
  set(BUILD_CAFFE2 OFF)
  set(CAFFE2_USE_CUDNN 1)

  # Set default CUDA architectures depending on CUDA version
  if(CUDA_VERSION_MAJOR LESS 11)
    set(TORCH_CUDA_ARCH_LIST "3.7;5.0;6.0;6.1;7.0;7.5" CACHE STRING "")
  elseif(CUDA_VERSION_MAJOR EQUAL 11)
    set(TORCH_CUDA_ARCH_LIST "5.0;6.0;6.1;7.0;7.5;8.0;8.6;8.9" CACHE STRING "")
  elseif(CUDA_VERSION_MAJOR EQUAL 12)
    set(TORCH_CUDA_ARCH_LIST "6.0;7.0;7.5;8.0;8.6;8.9;9.0" CACHE STRING "")
  elseif(CUDA_VERSION_MAJOR GREATER_EQUAL 13)
    set(TORCH_CUDA_ARCH_LIST "7.5;8.0;8.6;8.9;9.0;9.0+PTX" CACHE STRING "")
  else()
    message(WARNING "Unrecognized CUDA version: ${CUDA_VERSION}, using generic arch list")
    set(TORCH_CUDA_ARCH_LIST "7.5;8.0;8.6")
  endif()

  set(CMAKE_CUDA_ARCHITECTURES "all-major")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Python
# ─────────────────────────────────────────────────────────────────────────────
find_package(Python3 COMPONENTS Interpreter Development)

# ─────────────────────────────────────────────────────────────────────────────
# zlib
# ─────────────────────────────────────────────────────────────────────────────
find_package(ZLIB REQUIRED)

# ─────────────────────────────────────────────────────────────────────────────
# OpenBabel
# ─────────────────────────────────────────────────────────────────────────────
find_package(OpenBabel3 3.1.0...<4 REQUIRED)
set(OpenBabel3_LIBRARIES $<TARGET_FILE:openbabel>)
include_directories($<TARGET_PROPERTY:openbabel,INTERFACE_INCLUDE_DIRECTORIES>)
include_directories(${OpenBabel3_INCLUDE_DIRS})

# ─────────────────────────────────────────────────────────────────────────────
# Protobuf
# ─────────────────────────────────────────────────────────────────────────────
find_package(Protobuf REQUIRED)

# ─────────────────────────────────────────────────────────────────────────────
# Torch / libtorch
#
# USE_METAL=ON  → macOS ARM build with MPS backend.
#   PyTorch 2.10 macOS ARM includes MPS support out of the box.
#   URL pattern: https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-X.Y.Z.zip
#
# USE_METAL=OFF → CUDA build, fetches cu121 build same as before.
# ─────────────────────────────────────────────────────────────────────────────
if(GNINA_FORCE_EXTERNAL_LIBS)
  find_package(Torch 2 REQUIRED)
elseif(NOT APPLE)
  # On Apple Silicon we always use the FetchContent libtorch (below) to avoid
  # protobuf target conflicts between the system-installed protobuf and whatever
  # torch version may be installed via Homebrew / conda.
  find_package(Torch)
endif()

if(NOT Torch_FOUND)
  message(STATUS "Torch not found, fetching via FetchContent")
  include(FetchContent)

  if(USE_METAL)
    # PyTorch 2.10 macOS ARM — MPS backend included in the CPU build.
    # The MPS backend is automatically available when torch::Device("mps") is used.
    FetchContent_Declare(
      libtorch
      URL "https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-2.10.0.zip"
    )
  else()
    # CUDA 12.1 build for Linux (unchanged from original)
    FetchContent_Declare(
      libtorch
      URL "https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-2.4.1%2Bcu121.zip"
    )
  endif()

  FetchContent_MakeAvailable(libtorch)
  # Pin Torch_DIR to the fetched libtorch to prevent CMake from finding a
  # system-installed Torch (e.g. Homebrew) whose bundled Caffe2/protobuf
  # targets conflict with the system protobuf on macOS.
  set(Torch_DIR "${libtorch_SOURCE_DIR}/share/cmake/Torch" CACHE PATH "" FORCE)
  list(APPEND CMAKE_PREFIX_PATH "${libtorch_SOURCE_DIR}")
  find_package(Torch REQUIRED)

  if(USE_METAL)
    # Install macOS torch dylibs so the installed binary can find them via RPATH
    cmake_path(GET TORCH_LIBRARY PARENT_PATH TORCH_LIBDIR)
    set(TORCH_NEEDED
      libc10.dylib
      libomp.dylib
      libshm.dylib
      libtorch.dylib
      libtorch_cpu.dylib
      libtorch_global_deps.dylib
      libtorch_python.dylib
    )
    list(TRANSFORM TORCH_NEEDED PREPEND "${TORCH_LIBDIR}/")
    install(FILES ${TORCH_NEEDED} DESTINATION lib/gnina)
  elseif(NOT USE_METAL)
    # Install CUDA-specific torch shared libraries (Linux only)
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
      libcudnn_adv.so.9
      libcudnn_cnn.so.9
      libcudnn_engines_precompiled.so.9
      libcudnn_engines_runtime_compiled.so.9
      libcudnn_graph.so.9
      libcudnn_heuristic.so.9
      libcudnn_ops.so.9
      libcudnn.so.9
      libjitbackend_test.so
      libnnapi_backend.so
      libnvrtc-b51b459d.so.12
      libnvrtc-builtins.so
      libshm.so
      libtorchbind_test.so
      libtorch_cuda_linalg.so
      libtorch_global_deps.so
      libtorch_python.so
    )
    list(TRANSFORM TORCH_NEEDED PREPEND "${TORCH_LIBDIR}/")
    install(FILES ${TORCH_NEEDED} DESTINATION lib/gnina)
  endif()
else()
  message(STATUS "Torch found")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Strip GNU-ld-only flags from torch interface on macOS
#
# PyTorch's cmake adds -Xlinker --exclude-libs -Xlinker ALL to its interface
# link options. Apple ld does not support --exclude-libs and fails with
# "unknown options: --exclude-libs". Strip it from all torch cmake targets.
# ─────────────────────────────────────────────────────────────────────────────
if(APPLE)
  foreach(_torch_tgt torch torch_cpu c10)
    if(TARGET ${_torch_tgt})
      get_target_property(_iface_link_opts ${_torch_tgt} INTERFACE_LINK_OPTIONS)
      if(_iface_link_opts)
        # Remove the four tokens that make up the flag
        list(REMOVE_ITEM _iface_link_opts "LINKER:--exclude-libs,ALL")
        list(REMOVE_ITEM _iface_link_opts "-Wl,--exclude-libs,ALL")
        list(REMOVE_ITEM _iface_link_opts "--exclude-libs")
        list(REMOVE_ITEM _iface_link_opts "ALL")
        # Also handle the -Xlinker form (stored as a single semicolon-joined string)
        string(REPLACE "-Xlinker;--exclude-libs;-Xlinker;ALL" "" _iface_link_opts "${_iface_link_opts}")
        string(REPLACE "-Xlinker --exclude-libs -Xlinker ALL" "" _iface_link_opts "${_iface_link_opts}")
        set_target_properties(${_torch_tgt} PROPERTIES INTERFACE_LINK_OPTIONS "${_iface_link_opts}")
      endif()
    endif()
  endforeach()
endif()

# ─────────────────────────────────────────────────────────────────────────────
# libmolgrid
#
# USE_METAL=ON  → build libmolgrid without CUDA (USE_CUDA=OFF).
# USE_METAL=OFF → original CUDA build.
# ─────────────────────────────────────────────────────────────────────────────
if(GNINA_FORCE_EXTERNAL_LIBS)
  find_package(libmolgrid REQUIRED)
else()
  find_package(libmolgrid)
endif()

if(NOT LIBMOLGRID_LIBRARY)
  message("libmolgrid will be fetched from git")
  include(ExternalProject)
  find_package(Git REQUIRED)
  set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

  if(USE_METAL)
    # Write a patch script that will be applied to the fnachon/libmolgrid source.
    # Using a file avoids shell quoting nightmares in PATCH_COMMAND.
    file(WRITE ${CMAKE_BINARY_DIR}/patch_libmolgrid.sh
      "#!/bin/bash\n"
      "SRC=$1\n"
      "set -e\n"
      "# Disable unused subdirectories (not guarded by options in the fork's CMakeLists.txt)\n"
      "sed -i '' 's/add_subdirectory(python)/# add_subdirectory(python)/' \"$SRC/CMakeLists.txt\"\n"
      "sed -i '' 's/add_subdirectory(test)/# add_subdirectory(test)/' \"$SRC/CMakeLists.txt\"\n"
      "sed -i '' 's/add_subdirectory(docs)/# add_subdirectory(docs)/' \"$SRC/CMakeLists.txt\"\n"
    )
    execute_process(COMMAND chmod +x ${CMAKE_BINARY_DIR}/patch_libmolgrid.sh)

    ExternalProject_Add(libmolgrid
      GIT_REPOSITORY https://github.com/fnachon/libmolgrid
      UPDATE_DISCONNECTED TRUE
      PATCH_COMMAND ${CMAKE_BINARY_DIR}/patch_libmolgrid.sh <SOURCE_DIR>
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        -DLIBMOLGRID_ENABLE_CUDA=OFF
        -DLIBMOLGRID_ENABLE_METAL=ON
        -DOPENBABEL3_INCLUDE_DIR=/opt/homebrew/include/openbabel3
        -DOPENBABEL3_LIBRARIES=/opt/homebrew/lib/libopenbabel.dylib
    )
  else()
    # TEMPORARY: point at the fnachon/libmolgrid fork instead of upstream
    # gnina/libmolgrid. gninasrc (molgridder.cpp/h, torch_model.cpp) was
    # updated to use libmolgrid's Vec3/UVec2 types, which only exist in the
    # fork (see gnina/libmolgrid#127, not yet merged upstream). Switch this
    # back to gnina/libmolgrid once that PR lands.
    ExternalProject_Add(libmolgrid
      GIT_REPOSITORY https://github.com/fnachon/libmolgrid
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        -DLIBMOLGRID_ENABLE_CUDA=ON
        -DLIBMOLGRID_ENABLE_METAL=OFF
    )
  endif()

  set(LIBMOLGRID_INCLUDE_DIR ${EXTERNAL_INSTALL_LOCATION}/include)
  set(LIBMOLGRID_LIBRARY ${EXTERNAL_INSTALL_LOCATION}/lib/libmolgrid.a)
  include_directories(${LIBMOLGRID_INCLUDE_DIR})
  find_package(Python COMPONENTS Interpreter Development)
  set(PYTHONPATH $ENV{PYTHONPATH})
  set(ENV{PYTHONPATH}
    "${PYTHONPATH}:${EXTERNAL_INSTALL_LOCATION}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/")
else()
  add_custom_target(libmolgrid)
endif()

include_directories(${LIBMOLGRID_INCLUDE_DIR})
include_directories(.)

add_subdirectory(gninasrc)

# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────
option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  include(CTest)
  enable_testing()
  add_subdirectory(test)
endif()
