set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Boost REQUIRED COMPONENTS program_options system iostreams timer
    thread serialization filesystem date_time regex unit_test_framework)
find_package(OpenMP)

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

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${Boost_INCLUDE_DIRS}) 
include_directories(${CUDA_INCLUDE_DIRS})
# need to use eigen 3.2 until CUDA fixes are released
include_directories(../Eigen)

#lib
set(LIB_SRCS
lib/atom_constants.cpp
lib/bfgs.cu
lib/box.cpp
lib/builtinscoring.cpp
lib/cache.cpp
lib/cache_gpu.cpp
lib/cnn_scorer.cpp
lib/cnn_data.cpp
lib/coords.cpp
lib/custom_terms.cpp
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
lib/tree_gpu.cu
lib/user_opts.cpp
)

#grid member functions need to be called by arbitrary kernels
set_source_files_properties(lib/grid_gpu.cu PROPERTIES COMPILE_FLAGS -maxrregcount=64)

#test
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.c*)

add_library(gninalib ${LIB_SRCS})
set_target_properties(gninalib PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib caffe ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES})
# MESSAGE(STATUS "variable is: " ${NVCC_FLAGS_EXTRA})

add_library(gninalib_static STATIC ${LIB_SRCS})
set_target_properties(gninalib_static PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib_static caffe ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES})

install(TARGETS gninalib gninalib_static DESTINATION lib)

include_directories(lib)
include_directories(gninavis)
include_directories(${LIBMOLGRID_INCLUDE})

# compile in libgnina - there are enough dependencies to deal with
add_executable(gnina main/main.cpp ${LIB_SRCS})
target_link_libraries(gnina caffe ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES} ${LIBMOLGRID_LIBRARY})

add_subdirectory(gninaserver)

add_executable(fromgnina fromgnina/fromgnina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(fromgnina  caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES} )

add_executable(gninagrid gninagrid/gninagrid.cpp gninagrid/molgridder.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(gninagrid  caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES} ${CUDA_LIBRARIES})

add_executable(gninatyper gninatyper/gninatyper.cpp lib/CommandLine2/CommandLine.cpp ${LIB_SRCS})
target_link_libraries(gninatyper  caffe ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES} ${LIBMOLGRID_LIBRARY})

add_executable(tognina tognina/tognina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(tognina  caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES})

install(TARGETS gnina gninagrid gninatyper fromgnina tognina RUNTIME DESTINATION bin)

# gninavis uses rdkit, which can be a pain to install, so gracefully deal with its absence 
find_package(RDKit)

if(RDKIT_FOUND)
    include_directories(${RDKIT_INCLUDE_DIR})
    set(rdkit_lib_line "-lFileParsers -lSmilesParse -lDepictor -lGraphMol -lRDGeometryLib -lRDGeneral -lSubstructMatch -lSubgraphs -lDataStructs")
    
    add_executable(gninavis gninavis/gninavis.cpp lib/CommandLine2/CommandLine.cpp gninavis/cnn_visualization.cpp gninavis/cnn_visualization.hpp)
    target_link_libraries(gninavis caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES} ${RDKIT_LIBRARIES})
    install(TARGETS gninavis RUNTIME DESTINATION bin)
    
else()
    message(WARNING "RDKIT not found.  gninavis will NOT be built.")
endif()

