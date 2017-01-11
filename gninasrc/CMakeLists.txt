
find_package(CUDA REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options system iostreams timer thread serialization filesystem date_time regex)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  list(APPEND CUDA_NVCC_FLAGS -O3)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  list(APPEND CUDA_NVCC_FLAGS -O0)
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  list(APPEND CUDA_NVCC_FLAGS -g -lineinfo)
endif()

set(CUDA_SEPARABLE_COMPILATION "ON")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${Boost_INCLUDE_DIRS}) 
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${OPENBABEL2_INCLUDE_DIRS})
# need to use eigen 3.2 until CUDA fixes are released
include_directories(../Eigen)

#lib
set(LIB_SRCS
lib/atom_constants.cpp
lib/bfgs.cu
lib/box.cpp
lib/builtinscoring.cpp
lib/cache.cpp
lib/cnn_scorer.cpp
lib/coords.cpp
lib/custom_terms.cpp
lib/everything.cpp
lib/flexinfo.cpp
lib/GninaConverter.cpp
lib/grid.cpp
lib/molgetter.cpp
lib/monte_carlo.cpp
lib/mutate.cpp
lib/my_pid.cpp
lib/naive_non_cache.cpp
lib/nngridder.cpp
lib/non_cache.cpp
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
lib/conf_gpu.cu
lib/gpucode.cu
lib/gridmaker.cu
lib/model.cu
lib/non_cache_gpu.cu
lib/precalculate_gpu.cu
lib/tree_gpu.cu
)

cuda_add_library(gninalib ${LIB_SRCS} OPTIONS ${NVCC_FLAGS_EXTRA} --shared)
set_target_properties(gninalib PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})
install(TARGETS gninalib DESTINATION lib)

#cuda_add_library(gninalib_static STATIC ${LIB_SRCS} OPTIONS ${NVCC_FLAGS_EXTRA})
#set_target_properties(gninalib_static PROPERTIES OUTPUT_NAME gnina)
#target_link_libraries(gninalib_static caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

include_directories(lib)

# compile in libgnina - there are enough dependencies to deal with
cuda_add_executable(gnina main/main.cpp ${LIB_SRCS})
target_link_libraries(gnina caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

add_subdirectory(gninaserver)

cuda_add_executable(fromgnina fromgnina/fromgnina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(fromgnina  caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

cuda_add_executable(gninagrid gninagrid/gninagrid.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(gninagrid  caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

cuda_add_executable(gninatyper gninatyper/gninatyper.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(gninatyper  caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

cuda_add_executable(tognina tognina/tognina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(tognina  caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

install(TARGETS gnina gninagrid gninatyper fromgnina tognina RUNTIME DESTINATION bin)
