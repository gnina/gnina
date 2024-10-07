cmake_minimum_required(VERSION 3.25)
project (gnina C CXX CUDA)


cmake_policy(SET CMP0104 NEW) 
cmake_policy(SET CMP0135 NEW)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
cmake_policy(SET CMP0094 NEW) #location first python finding
cmake_policy(SET CMP0146 OLD) #too baked in

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
option(GNINA_FORCE_EXTERNAL_LIBS "Use external libraries only" OFF)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_INSTALL_PREFIX}/lib/gnina")
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -Wno-unknown-pragmas")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (BUILD_COVERAGE)
  set(CMAKE_CXX_FLAGS_RELEASE "-fprofile-arcs -ftest-coverage")
endif()

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA 12.0 REQUIRED)
set(BUILD_CAFFE2 OFF)
set(CAFFE2_USE_CUDNN 1)

function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if ((NOT DEFINED ARGV0) OR _variableName MATCHES ${ARGV0})
            message(STATUS "${_variableName}=${${_variableName}}")
        endif()
    endforeach()
endfunction()

#python
find_package(Python3  COMPONENTS Interpreter Development)

#zlib
find_package(ZLIB REQUIRED)

#openbabel - in actuality need 3.1.1, but the 3-1-1 tag is incorrectly versioned
find_package(OpenBabel3 3.1.0...<4 REQUIRED)
set(OpenBabel3_LIBRARIES $<TARGET_FILE:openbabel>)
include_directories($<TARGET_PROPERTY:openbabel,INTERFACE_INCLUDE_DIRECTORIES>)
include_directories(${OpenBabel3_INCLUDE_DIRS})

#make torch use system protobuf to be compat with caffe
find_package(Protobuf REQUIRED)

#In order to build a mostly static executible, first build torch:
# cmake .. -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.0 -DCMAKE_CUDA_ARCHITECTURES="60;61;70;75;80;90" -DUSE_STATIC_MKL=ON -DUSE_MKLDNN=ON -DBLAS=OpenBLAS -DBUILD_SHARED_LIBS=OFF -DBUILD_PYTHON=OFF -DCAFFE2_STATIC_LINK_CUDA=ON -DUSE_TENSORPIPE=OFF -DUSE_STATIC_NCCL=ON -DUSE_STATIC_CUDNN=ON -DTORCH_CUDA_ARCH_LIST="6.0 6.1 7.0 7.5 8.0 9.0"
# need to build static mkl-dnn library with cmake ..  -DDNNL_LIBRARY_TYPE=STATIC
# configure gnina with cmake -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.0/ -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.0/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=all ..
# the commandline to build a mostly static (CUDA and standard system dependencies only) executable on my machine:
# /usr/bin/c++ -Wno-deprecated-declarations -Wno-unknown-pragmas -D_GLIBCXX_USE_CXX11_ABI=1 -O3 -DNDEBUG -Wl,--disable-new-dtags CMakeFiles/gnina.dir/main/main.cpp.o CMakeFiles/gninalibobj.dir/version.cpp.o CMakeFiles/gninalibobj.dir/lib/atom_constants.cpp.o CMakeFiles/gninalibobj.dir/lib/bfgs.cu.o CMakeFiles/gninalibobj.dir/lib/box.cpp.o CMakeFiles/gninalibobj.dir/lib/builtinscoring.cpp.o CMakeFiles/gninalibobj.dir/lib/cache.cpp.o CMakeFiles/gninalibobj.dir/lib/cache_gpu.cpp.o CMakeFiles/gninalibobj.dir/lib/cnn_torch_scorer.cpp.o CMakeFiles/gninalibobj.dir/lib/coords.cpp.o CMakeFiles/gninalibobj.dir/lib/covinfo.cpp.o CMakeFiles/gninalibobj.dir/lib/custom_terms.cpp.o CMakeFiles/gninalibobj.dir/lib/dl_scorer.cpp.o CMakeFiles/gninalibobj.dir/lib/device_buffer.cpp.o CMakeFiles/gninalibobj.dir/lib/everything.cpp.o CMakeFiles/gninalibobj.dir/lib/flexinfo.cpp.o CMakeFiles/gninalibobj.dir/lib/GninaConverter.cpp.o CMakeFiles/gninalibobj.dir/lib/grid.cpp.o CMakeFiles/gninalibobj.dir/lib/grid_gpu.cu.o CMakeFiles/gninalibobj.dir/lib/model.cpp.o CMakeFiles/gninalibobj.dir/lib/molgetter.cpp.o CMakeFiles/gninalibobj.dir/lib/monte_carlo.cpp.o CMakeFiles/gninalibobj.dir/lib/mutate.cpp.o CMakeFiles/gninalibobj.dir/lib/my_pid.cpp.o CMakeFiles/gninalibobj.dir/lib/naive_non_cache.cpp.o CMakeFiles/gninalibobj.dir/lib/non_cache.cpp.o CMakeFiles/gninalibobj.dir/lib/non_cache_cnn.cpp.o CMakeFiles/gninalibobj.dir/lib/obmolopener.cpp.o CMakeFiles/gninalibobj.dir/lib/parallel_mc.cpp.o CMakeFiles/gninalibobj.dir/lib/parallel_progress.cpp.o CMakeFiles/gninalibobj.dir/lib/parse_pdbqt.cpp.o CMakeFiles/gninalibobj.dir/lib/pdb.cpp.o CMakeFiles/gninalibobj.dir/lib/PDBQTUtilities.cpp.o CMakeFiles/gninalibobj.dir/lib/quasi_newton.cpp.o CMakeFiles/gninalibobj.dir/lib/quaternion.cu.o CMakeFiles/gninalibobj.dir/lib/random.cpp.o CMakeFiles/gninalibobj.dir/lib/result_info.cpp.o CMakeFiles/gninalibobj.dir/lib/ssd.cpp.o CMakeFiles/gninalibobj.dir/lib/szv_grid.cpp.o CMakeFiles/gninalibobj.dir/lib/terms.cpp.o CMakeFiles/gninalibobj.dir/lib/weighted_terms.cpp.o CMakeFiles/gninalibobj.dir/lib/conf.cpp.o CMakeFiles/gninalibobj.dir/lib/conf_gpu.cu.o CMakeFiles/gninalibobj.dir/lib/gpucode.cu.o CMakeFiles/gninalibobj.dir/lib/model.cu.o CMakeFiles/gninalibobj.dir/lib/non_cache_gpu.cu.o CMakeFiles/gninalibobj.dir/lib/precalculate_gpu.cu.o CMakeFiles/gninalibobj.dir/lib/torch_model.cpp.o CMakeFiles/gninalibobj.dir/lib/tree_gpu.cu.o CMakeFiles/gninalibobj.dir/lib/user_opts.cpp.o CMakeFiles/gninalibobj.dir/torch_models.cpp.o CMakeFiles/gnina.dir/cmake_device_link.o -o ../bin/gnina   -L/usr/local/cuda-12.0/targets/x86_64-linux/lib/stubs  -L/usr/local/cuda-12.0/targets/x86_64-linux/lib  -Wl,-rpath,/usr/local/lib:/usr/local/cuda-12.0/lib64: all_default_to_default_1.3_1.o all_default_to_default_1.3_2.o all_default_to_default_1.3_3.o crossdock_default2018_1.3_1.o crossdock_default2018_1.3_2.o crossdock_default2018_1.3_3.o crossdock_default2018_1.3_4.o crossdock_default2018_1.3.o crossdock_default2018_1.o crossdock_default2018_2.o crossdock_default2018_3.o crossdock_default2018_4.o crossdock_default2018_KD_1.o crossdock_default2018_KD_2.o crossdock_default2018_KD_3.o crossdock_default2018_KD_4.o crossdock_default2018_KD_5.o crossdock_default2018.o default2017.o dense_1.3_1.o dense_1.3_2.o dense_1.3_3.o dense_1.3_4.o dense_1.3.o dense_1.3_PT_KD_1.o dense_1.3_PT_KD_2.o dense_1.3_PT_KD_3.o dense_1.3_PT_KD_4.o dense_1.3_PT_KD_def2018_1.o dense_1.3_PT_KD_def2018_2.o dense_1.3_PT_KD_def2018_3.o dense_1.3_PT_KD_def2018_4.o dense_1.3_PT_KD_def2018.o dense_1.3_PT_KD.o dense_1.o dense_2.o dense_3.o dense_4.o dense.o general_default2018_1.o general_default2018_2.o general_default2018_3.o general_default2018_4.o general_default2018_KD_1.o general_default2018_KD_2.o general_default2018_KD_3.o general_default2018_KD_4.o general_default2018_KD_5.o general_default2018.o redock_default2018_1.3_1.o redock_default2018_1.3_2.o redock_default2018_1.3_3.o redock_default2018_1.3_4.o redock_default2018_1.3.o redock_default2018_1.o redock_default2018_2.o redock_default2018_3.o redock_default2018_4.o redock_default2018_KD_1.o redock_default2018_KD_2.o redock_default2018_KD_3.o redock_default2018_KD_4.o redock_default2018_KD_5.o redock_default2018.o /usr/lib/x86_64-linux-gnu/libboost_program_options.a /usr/lib/x86_64-linux-gnu/libboost_system.a /usr/lib/x86_64-linux-gnu/libboost_iostreams.a /usr/lib/x86_64-linux-gnu/libboost_timer.a /usr/lib/x86_64-linux-gnu/libboost_thread.a /usr/lib/x86_64-linux-gnu/libboost_serialization.a /usr/lib/x86_64-linux-gnu/libboost_filesystem.a /usr/lib/x86_64-linux-gnu/libboost_date_time.a /usr/lib/x86_64-linux-gnu/libboost_regex.a /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.a /usr/local/lib/libopenbabel.a /usr/local/lib/libmolgrid.a -Wl,--whole-archive /usr/local/lib/libtorch.a -Wl,--no-whole-archive -Wl,--whole-archive /usr/local/lib/libtorch_cpu.a /usr/lib/gcc/x86_64-linux-gnu/11/libgomp.a        -Wl,--no-whole-archive -Wl,--whole-archive /usr/local/lib/libtorch_cuda.a -Wl,--no-whole-archive -Wl,--whole-archive -L/usr/local/cuda-12.0/lib64 -lcudart -lcusparse -lcufft -lnvToolsExt -lcurand -lcublas -lcublasLt  -lcusolver /home/dkoes/git/pytorch/build-static/nccl/lib/libnccl_static.a /usr/local/lib/libc10_cuda.a -Wl,--no-whole-archive /usr/local/lib/libc10.a /usr/local/lib/libnnpack.a /lib/x86_64-linux-gnu/libnuma.a /usr/local/lib/libpytorch_qnnpack.a /usr/local/lib/libXNNPACK.a /usr/local/lib/libprotobuf-lite.a /usr/local/lib/libprotobuf.a /usr/local/lib/libprotoc.a /usr/local/lib/libfmt.a /usr/local/lib/libcpuinfo.a /usr/local/lib/libclog.a /usr/local/lib/libpthreadpool.a /usr/local/lib/libfbgemm.a /usr/local/lib/libdnnl.a /usr/local/lib/libdnnl.a /usr/local/lib/libsleef.a /usr/local/lib/libasmjit.a /usr/local/lib/libkineto.a /usr/local/cuda-12.0/lib64/libcupti_static.a -lnvToolsExt -lcudart_static -ldl /usr/lib/x86_64-linux-gnu/librt.a /usr/lib/x86_64-linux-gnu/libz.a /usr/lib/x86_64-linux-gnu/libboost_chrono.a /home/dkoes/git/pytorch/build-static/lib/libgloo.a /usr/local/lib/gnina/libgloo_cuda.a /usr/lib/x86_64-linux-gnu/libboost_atomic.a /usr/local/lib/libjsoncpp.a /home/dkoes/git/pytorch/build-static/lib/libonnx.a /home/dkoes/git/pytorch/build-static/lib/libonnx_proto.a /usr/lib/x86_64-linux-gnu/libopenblas.a /usr/lib/gcc/x86_64-linux-gnu/11/libgfortran.a /usr/lib/gcc/x86_64-linux-gnu/11/libquadmath.a -lcudadevrt -lcudart_static 

if(GNINA_FORCE_EXTERNAL_LIBS)
  find_package(Torch 2 REQUIRED)
else()
  find_package(Torch)
endif()

set(CMAKE_CUDA_ARCHITECTURES "all-major")

if (NOT Torch_FOUND)
    message(STATUS "Torch not found, fetching via FetchContent")

  include(FetchContent)
  #https://discuss.pytorch.org/t/using-fetchcontent-declare-in-cmake-to-pull-libtorch/177449
  FetchContent_Declare(
          libtorch
          URL "https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-2.4.1%2Bcu121.zip"
  )

  FetchContent_MakeAvailable(libtorch)
  list(APPEND CMAKE_PREFIX_PATH "${libtorch_SOURCE_DIR}")
  find_package(Torch REQUIRED)

  #the following may need to be updated if the version of torch is changed
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

  #install in subdirectory to avoid messing up a system torch install
  install(FILES ${TORCH_NEEDED} DESTINATION lib/gnina)
else()
  message(STATUS "Torch found")
endif()



#libmolgrid  
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
