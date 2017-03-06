set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(CUDA REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options system iostreams timer thread serialization filesystem date_time regex)
find_package(RDKit REQUIRED)

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

set(CUDA_SEPARABLE_COMPILATION "ON")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${Boost_INCLUDE_DIRS}) 
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${OPENBABEL2_INCLUDE_DIRS})
include_directories(${RDKIT_INCLUDE_DIR})
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
gninavis/cnn_visualization.cpp
)

set(rdkit_lib_line "-lFileParsers -lSmilesParse -lDepictor -lGraphMol -lRDGeometryLib -lRDGeneral -lSubstructMatch -lSubgraphs -lDataStructs")

find_library(CUDA_DEV_RT NAMES cudadevrt PATHS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
cuda_add_library(gninalib ${LIB_SRCS} OPTIONS ${NVCC_FLAGS_EXTRA} --shared)
set_target_properties(gninalib PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES} ${RDKIT_LIBRARIES})
# MESSAGE(STATUS "variable is: " ${NVCC_FLAGS_EXTRA})

install(TARGETS gninalib DESTINATION lib)

#cuda_add_library(gninalib_static STATIC ${LIB_SRCS} OPTIONS ${NVCC_FLAGS_EXTRA})
#set_target_properties(gninalib_static PROPERTIES OUTPUT_NAME gnina)
#target_link_libraries(gninalib_static caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

include_directories(lib)
include_directories(gninavis)

# compile in libgnina - there are enough dependencies to deal with
cuda_add_executable(gnina main/main.cpp ${LIB_SRCS})
target_link_libraries(gnina caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES}
    ${RDKIT_LIBRARIES} ${CUDA_DEV_RT})

add_subdirectory(gninaserver)

cuda_add_executable(fromgnina fromgnina/fromgnina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(fromgnina  caffe gninalib ${Boost_LIBRARIES}
    ${OPENBABEL2_LIBRARIES} ${CUDA_DEV_RT})

cuda_add_executable(gninagrid gninagrid/gninagrid.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(gninagrid  caffe gninalib ${Boost_LIBRARIES}
    ${OPENBABEL2_LIBRARIES} ${RDKIT_LIBRARIES} ${CUDA_DEV_RT})

cuda_add_executable(gninatyper gninatyper/gninatyper.cpp lib/CommandLine2/CommandLine.cpp ${LIB_SRCS})
target_link_libraries(gninatyper  caffe ${Boost_LIBRARIES}
    ${OPENBABEL2_LIBRARIES} ${RDKIT_LIBRARIES} ${CUDA_DEV_RT})

cuda_add_executable(tognina tognina/tognina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(tognina  caffe gninalib ${Boost_LIBRARIES}
    ${OPENBABEL2_LIBRARIES} ${CUDA_DEV_RT})

cuda_add_executable(gninavis gninavis/gninavis.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(gninavis caffe gninalib ${Boost_LIBRARIES}
    ${OPENBABEL2_LIBRARIES} ${RDKIT_LIBRARIES} ${CUDA_DEV_RT})

install(TARGETS gnina gninagrid gninatyper fromgnina tognina gninavis RUNTIME DESTINATION bin)
