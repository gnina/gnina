set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Boost REQUIRED COMPONENTS program_options iostreams timer
    thread serialization filesystem date_time regex unit_test_framework)
# boost_system is header-only since Boost 1.69 and no longer has a cmake config
find_package(OpenMP)
find_package(jsoncpp REQUIRED)

# ─────────────────────────────────────────────────────────────────────────────
# Compiler flags
# ─────────────────────────────────────────────────────────────────────────────
if(NOT USE_METAL)
  # CUDA NVCC flags (Linux/NVIDIA only)
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    list(APPEND CUDA_NVCC_FLAGS -O3 --default-stream per-thread -Xptxas -dlcm=ca)
  elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    list(APPEND CUDA_NVCC_FLAGS -O0 -G --default-stream per-thread)
  endif()
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    list(APPEND CUDA_NVCC_FLAGS -g)
  endif()
  if("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    list(APPEND CUDA_NVCC_FLAGS -lineinfo)
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${Boost_INCLUDE_DIRS})
if(NOT USE_METAL)
  include_directories(${CUDA_INCLUDE_DIRS})
endif()
include_directories(${TORCH_INCLUDE_DIRS})
if(USE_METAL)
  # On macOS, Homebrew also installs PyTorch headers under /opt/homebrew/include.
  # When both Homebrew (-isystem /opt/homebrew/include) and fetched libtorch are
  # in the search path, Clang deduplicates dirs that appear as both -I and -isystem
  # and may search Homebrew's torch/csrc/api/include/torch/ordered_dict.h before
  # the fetched one, causing a "redefinition of OrderedDict" error.
  #
  # Fix: generate an override header that provides an unconditional #pragma once
  # and then includes the fetched libtorch's ordered_dict.h by absolute path.
  # The override dir is added first, so any include of
  #   <torch/csrc/api/include/torch/ordered_dict.h>
  # hits the override first, which ensures the class is defined exactly once.
  # On macOS, Homebrew installs PyTorch headers under /opt/homebrew/include/torch/
  # which conflicts with the fetched libtorch (different version). Both end up as
  # -isystem dirs; Clang deduplicates user -I dirs that match existing -isystem dirs,
  # so the Homebrew version wins for torch/csrc/api/include/torch/*.h paths.
  #
  # Fix: generate an override directory whose headers provide #pragma once and
  # then include the fetched libtorch version by absolute path. The override dir
  # is prepended to CMAKE_CXX_FLAGS as a user -I dir (searched before -isystem).
  set(_override_base "${CMAKE_BINARY_DIR}/include_overrides")

  # Find which TORCH_INCLUDE_DIRS contains torch/csrc/api/include/torch/*.h
  set(_fetched_torch_api_dir "")
  foreach(_d ${TORCH_INCLUDE_DIRS})
    if(EXISTS "${_d}/torch/ordered_dict.h")
      set(_fetched_torch_api_dir "${_d}/torch")
      break()
    endif()
  endforeach()
  if(NOT _fetched_torch_api_dir)
    message(FATAL_ERROR "Cannot find fetched libtorch API dir in: ${TORCH_INCLUDE_DIRS}")
  endif()

  # Generate override headers for all .h files in torch/csrc/api/include/torch/
  set(_override_torch_dir "${_override_base}/torch/csrc/api/include/torch")
  file(MAKE_DIRECTORY "${_override_torch_dir}")
  file(GLOB _fetched_api_headers "${_fetched_torch_api_dir}/*.h")
  foreach(_hdr ${_fetched_api_headers})
    cmake_path(GET _hdr FILENAME _hdr_name)
    file(WRITE "${_override_torch_dir}/${_hdr_name}"
      "// Auto-generated — redirects to fetched libtorch to prevent Homebrew conflict.\n"
      "#pragma once\n"
      "#include \"${_fetched_torch_api_dir}/${_hdr_name}\"\n"
    )
  endforeach()

  set(CMAKE_CXX_FLAGS "-I${_override_base} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_OBJCXX_FLAGS "-I${_override_base} ${CMAKE_OBJCXX_FLAGS}")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Source files
#
# The 9 .cu files are compiled as:
#   USE_METAL=ON  → plain C++ (LANGUAGE CXX), using cuda_metal_compat.h stubs.
#                   Kernel <<<>>> launches are guarded by #ifndef USE_METAL.
#   USE_METAL=OFF → CUDA (compiled by NVCC as normal).
# ─────────────────────────────────────────────────────────────────────────────
set(LIB_SRCS
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
  lib/atom_constants.cpp
  lib/bfgs.cu
  lib/box.cpp
  lib/builtinscoring.cpp
  lib/cache.cpp
  lib/cache_gpu.cpp
  lib/cnn_torch_scorer.cpp
  lib/coords.cpp
  lib/covinfo.cpp
  lib/custom_terms.cpp
  lib/dl_scorer.cpp
  lib/device_buffer.cpp
  lib/everything.cpp
  lib/flexinfo.cpp
  lib/GninaConverter.cpp
  lib/grid.cpp
  lib/grid_gpu.cu
  lib/model.cpp
  lib/molgetter.cpp
  lib/monte_carlo.cpp
  lib/mutate.cpp
  lib/my_pid.cpp
  lib/naive_non_cache.cpp
  lib/non_cache.cpp
  lib/non_cache_cnn.cpp
  lib/obmolopener.cpp
  lib/parallel_mc.cpp
  lib/parallel_progress.cpp
  lib/parse_pdbqt.cpp
  lib/pdb.cpp
  lib/PDBQTUtilities.cpp
  lib/quasi_newton.cpp
  lib/quaternion.cu
  lib/random.cpp
  lib/result_info.cpp
  lib/ssd.cpp
  lib/szv_grid.cpp
  lib/terms.cpp
  lib/weighted_terms.cpp
  lib/conf.cpp
  lib/conf_gpu.cu
  lib/gpucode.cu
  lib/model.cu
  lib/non_cache_gpu.cu
  lib/precalculate_gpu.cu
  lib/torch_model.cpp
  lib/tree_gpu.cu
  lib/user_opts.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/torch_models.cpp"
)

find_package(Python3)   # needed for Metal shader embedding below

if(USE_METAL)
  # ── Metal shader compilation ────────────────────────────────────────────────
  # Step 1: compile gnina_kernels.metal → gnina_kernels.air
  add_custom_command(
    OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.air"
    COMMAND xcrun -sdk macosx metal
              -c "${CMAKE_CURRENT_SOURCE_DIR}/lib/gnina_kernels.metal"
              -o "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.air"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/lib/gnina_kernels.metal"
    COMMENT "Compiling Metal shaders → gnina_kernels.air"
  )

  # Step 2: link .air → gnina_kernels.metallib
  add_custom_command(
    OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.metallib"
    COMMAND xcrun -sdk macosx metallib
              "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.air"
              -o "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.metallib"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.air"
    COMMENT "Linking Metal library → gnina_kernels.metallib"
  )

  # Step 3: embed .metallib bytes as a C++ array using embed_binary.py
  add_custom_command(
    OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/gnina_metallib_embed.cpp"
    COMMAND ${Python3_EXECUTABLE}
              "${CMAKE_CURRENT_SOURCE_DIR}/embed_binary.py"
              "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.metallib"
              "${CMAKE_CURRENT_BINARY_DIR}/gnina_metallib_embed.cpp"
              "gnina_metallib"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gnina_kernels.metallib"
            "${CMAKE_CURRENT_SOURCE_DIR}/embed_binary.py"
    COMMENT "Embedding Metal library as C++ byte array"
  )

  # Compile metal_context.mm (Objective-C++ Metal runtime)
  list(APPEND LIB_SRCS
    lib/metal_context.mm
    "${CMAKE_CURRENT_BINARY_DIR}/gnina_metallib_embed.cpp"
  )

  # Tell CMake to treat all .cu files as C++ — cuda_metal_compat.h provides
  # the necessary stubs so Clang can compile them without NVCC.
  set(CUDA_SRCS
    lib/bfgs.cu
    lib/grid_gpu.cu
    lib/quaternion.cu
    lib/conf_gpu.cu
    lib/gpucode.cu
    lib/model.cu
    lib/non_cache_gpu.cu
    lib/precalculate_gpu.cu
    lib/tree_gpu.cu
  )
  foreach(_cu_file ${CUDA_SRCS})
    set_source_files_properties(${_cu_file} PROPERTIES LANGUAGE CXX)
  endforeach()

  # device_buffer.cpp uses #import <Metal/Metal.h> and ObjC bridge syntax —
  # must be compiled as Objective-C++.
  set_source_files_properties(lib/device_buffer.cpp PROPERTIES LANGUAGE OBJCXX)

  # grid_gpu.cu used to have -maxrregcount=64 for NVCC; not applicable to Clang.
  # The Metal shader equivalent will be set in Phase 4.
else()
  # CUDA: grid_gpu.cu must be callable from arbitrary kernels → limit registers.
  set_source_files_properties(lib/grid_gpu.cu PROPERTIES COMPILE_FLAGS -maxrregcount=64)
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Embedded PyTorch models (.pt → .o → torch_models.cpp)
# ─────────────────────────────────────────────────────────────────────────────
find_package(Python3)

function(model_to_obj out_var)
  set(result)
  foreach(in_f ${ARGN})
    cmake_path(GET in_f FILENAME out_f)
    string(REGEX REPLACE ".pt$" ".o" out_f ${out_f})
    set(out_f "${CMAKE_CURRENT_BINARY_DIR}/${out_f}")
    if(APPLE)
      # macOS ld does not support -b binary; use a Python script instead to
      # generate a .cpp with the same symbol names that GNU ld would produce,
      # then compile it to an object file.
      cmake_path(GET in_f STEM LAST_ONLY stem_f)  # e.g. "all_default_to_default_1.3_1"
      string(REGEX REPLACE "\\.o$" ".cpp" cpp_f ${out_f})
      add_custom_command(OUTPUT ${cpp_f}
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/embed_model.py"
                "${CMAKE_CURRENT_SOURCE_DIR}/${in_f}"
                "${cpp_f}"
                "${stem_f}"
        DEPENDS ${in_f} "${CMAKE_CURRENT_SOURCE_DIR}/embed_model.py"
        COMMENT "Embedding model ${in_f}"
        VERBATIM
      )
      add_custom_command(OUTPUT ${out_f}
        COMMAND ${CMAKE_CXX_COMPILER} -c ${cpp_f} -o ${out_f}
        DEPENDS ${cpp_f}
        COMMENT "Compiling embedded model ${out_f}"
        VERBATIM
      )
    else()
      add_custom_command(OUTPUT ${out_f}
        COMMAND ${CMAKE_LINKER} -r -b binary -o ${out_f} ${in_f}
        DEPENDS ${in_f}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Creating model object file ${out_f}"
        VERBATIM
      )
    endif()
    list(APPEND result ${out_f})
  endforeach()

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/torch_models.cpp"
    COMMAND ${Python3_EXECUTABLE} make_model_cpp.py "${result}" "${CMAKE_CURRENT_BINARY_DIR}/torch_models.cpp" $<IF:$<BOOL:${APPLE}>,apple,other>
    DEPENDS ${result} make_model_cpp.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Creating ${CMAKE_CURRENT_BINARY_DIR}/torch_models.cpp"
    VERBATIM
  )
  set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

set(MODELS
  lib/models/all_default_to_default_1.3_1.pt
  lib/models/all_default_to_default_1.3_2.pt
  lib/models/all_default_to_default_1.3_3.pt
  lib/models/crossdock_default2018_1.3_1.pt
  lib/models/crossdock_default2018_1.3_2.pt
  lib/models/crossdock_default2018_1.3_3.pt
  lib/models/crossdock_default2018_1.3_4.pt
  lib/models/crossdock_default2018_1.3.pt
  lib/models/crossdock_default2018_1.pt
  lib/models/crossdock_default2018_2.pt
  lib/models/crossdock_default2018_3.pt
  lib/models/crossdock_default2018_4.pt
  lib/models/crossdock_default2018_KD_1.pt
  lib/models/crossdock_default2018_KD_2.pt
  lib/models/crossdock_default2018_KD_3.pt
  lib/models/crossdock_default2018_KD_4.pt
  lib/models/crossdock_default2018_KD_5.pt
  lib/models/crossdock_default2018.pt
  lib/models/default2017.pt
  lib/models/dense_1.3_1.pt
  lib/models/dense_1.3_2.pt
  lib/models/dense_1.3_3.pt
  lib/models/dense_1.3_4.pt
  lib/models/dense_1.3.pt
  lib/models/dense_1.3_PT_KD_1.pt
  lib/models/dense_1.3_PT_KD_2.pt
  lib/models/dense_1.3_PT_KD_3.pt
  lib/models/dense_1.3_PT_KD_4.pt
  lib/models/dense_1.3_PT_KD_def2018_1.pt
  lib/models/dense_1.3_PT_KD_def2018_2.pt
  lib/models/dense_1.3_PT_KD_def2018_3.pt
  lib/models/dense_1.3_PT_KD_def2018_4.pt
  lib/models/dense_1.3_PT_KD_def2018.pt
  lib/models/dense_1.3_PT_KD.pt
  lib/models/dense_1.pt
  lib/models/dense_2.pt
  lib/models/dense_3.pt
  lib/models/dense_4.pt
  lib/models/dense.pt
  lib/models/general_default2018_1.pt
  lib/models/general_default2018_2.pt
  lib/models/general_default2018_3.pt
  lib/models/general_default2018_4.pt
  lib/models/general_default2018_KD_1.pt
  lib/models/general_default2018_KD_2.pt
  lib/models/general_default2018_KD_3.pt
  lib/models/general_default2018_KD_4.pt
  lib/models/general_default2018_KD_5.pt
  lib/models/general_default2018.pt
  lib/models/redock_default2018_1.3_1.pt
  lib/models/redock_default2018_1.3_2.pt
  lib/models/redock_default2018_1.3_3.pt
  lib/models/redock_default2018_1.3_4.pt
  lib/models/redock_default2018_1.3.pt
  lib/models/redock_default2018_1.pt
  lib/models/redock_default2018_2.pt
  lib/models/redock_default2018_3.pt
  lib/models/redock_default2018_4.pt
  lib/models/redock_default2018_KD_1.pt
  lib/models/redock_default2018_KD_2.pt
  lib/models/redock_default2018_KD_3.pt
  lib/models/redock_default2018_KD_4.pt
  lib/models/redock_default2018_KD_5.pt
  lib/models/redock_default2018.pt
)

model_to_obj(MODEL_OBJS "${MODELS}")

# ─────────────────────────────────────────────────────────────────────────────
# Version stamp
# ─────────────────────────────────────────────────────────────────────────────
add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/version.cmake)
add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
  PROPERTIES SYMBOLIC TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Library targets
# ─────────────────────────────────────────────────────────────────────────────
add_library(gninalibobj OBJECT ${LIB_SRCS})
set_property(TARGET gninalibobj PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(gninalibobj ${MODEL_OBJS} jsoncpp_lib)
add_dependencies(gninalibobj libmolgrid torch)

# Metal-specific: link Apple frameworks
if(USE_METAL)
  target_link_libraries(gninalibobj
    "-framework Metal"
    "-framework Foundation"
    "-framework CoreFoundation"
  )
  # Ensure -DUSE_METAL propagates to all consumers
  target_compile_definitions(gninalibobj PUBLIC USE_METAL)
endif()

add_library(gninalib_static STATIC $<TARGET_OBJECTS:gninalibobj>)
set_target_properties(gninalib_static PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib_static
  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${LIBMOLGRID_LIBRARY}
  ${TORCH_LIBRARIES} jsoncpp_lib ZLIB::ZLIB gninalibobj)

add_library(gninalib SHARED $<TARGET_OBJECTS:gninalibobj>)
set_target_properties(gninalib PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib
  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${LIBMOLGRID_LIBRARY}
  ${TORCH_LIBRARIES} jsoncpp_lib ZLIB::ZLIB gninalibobj)

if(USE_METAL)
  foreach(_tgt gninalib_static gninalib)
    target_link_libraries(${_tgt}
      "-framework Metal"
      "-framework Foundation"
      "-framework CoreFoundation"
    )
  endforeach()
endif()

install(TARGETS gninalib gninalib_static DESTINATION lib)

include_directories(lib)
include_directories(gninavis)
include_directories(${LIBMOLGRID_INCLUDE})

# ─────────────────────────────────────────────────────────────────────────────
# Executables
# ─────────────────────────────────────────────────────────────────────────────
add_executable(gnina main/main.cpp $<TARGET_OBJECTS:gninalibobj>)
target_link_libraries(gnina
  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${LIBMOLGRID_LIBRARY}
  ${TORCH_LIBRARIES} jsoncpp_lib ZLIB::ZLIB gninalibobj)
if(USE_METAL)
  target_link_libraries(gnina
    "-framework Metal" "-framework Foundation" "-framework CoreFoundation")
endif()

add_subdirectory(gninaserver)

option(BUILD_PYGNINA "Build the pygnina Python extension module" ON)
if(BUILD_PYGNINA)
  add_subdirectory(pygnina)
endif()

add_executable(fromgnina fromgnina/fromgnina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(fromgnina gninalib)

add_executable(gninagrid gninagrid/gninagrid.cpp gninagrid/molgridder.cpp lib/CommandLine2/CommandLine.cpp)
if(USE_METAL)
  target_link_libraries(gninagrid gninalib ${LIBMOLGRID_LIBRARY})
else()
  target_link_libraries(gninagrid gninalib ${CUDA_LIBRARIES} ${LIBMOLGRID_LIBRARY})
endif()

add_executable(gninatyper gninatyper/gninatyper.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(gninatyper
  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${LIBMOLGRID_LIBRARY}
  ${TORCH_LIBRARIES} jsoncpp_lib gninalib)

add_executable(tognina tognina/tognina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(tognina gninalib ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES})

install(TARGETS gnina gninagrid gninatyper fromgnina tognina RUNTIME DESTINATION bin)

# ─────────────────────────────────────────────────────────────────────────────
# gninavis (optional — requires RDKit)
# ─────────────────────────────────────────────────────────────────────────────
find_package(RDKit)
if(RDKIT_FOUND)
  include_directories(${RDKIT_INCLUDE_DIR})
  add_executable(gninavis gninavis/gninavis.cpp lib/CommandLine2/CommandLine.cpp
    gninavis/cnn_visualization.cpp gninavis/cnn_visualization.hpp)
  target_link_libraries(gninavis gninalib ${LIBMOLGRID_LIBRARY}
    ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${RDKIT_LIBRARIES})
  install(TARGETS gninavis RUNTIME DESTINATION bin)
else()
  message(WARNING "RDKIT not found. gninavis will NOT be built.")
endif()
