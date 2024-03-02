set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Boost REQUIRED COMPONENTS program_options system iostreams timer
    thread serialization filesystem date_time regex unit_test_framework)
find_package(OpenMP)

find_package(jsoncpp REQUIRED)

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${Boost_INCLUDE_DIRS}) 
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${TORCH_INCLUDE_DIRS})
include_directories(${JSON_INC_PATH})

# need to use eigen 3.2 until CUDA fixes are released
include_directories(../Eigen)

#lib
set(LIB_SRCS
${CMAKE_CURRENT_BINARY_DIR}/version.cpp
lib/atom_constants.cpp
lib/bfgs.cu
lib/box.cpp
lib/builtinscoring.cpp
lib/cache.cpp
lib/cache_gpu.cpp
lib/cnn_scorer.cpp
lib/cnn_torch_scorer.cpp
lib/cnn_data.cpp
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
  
#models
function(model_to_obj out_var)
  set(result)
  foreach(in_f ${ARGN})
    cmake_path(GET in_f FILENAME out_f)
    string(REGEX REPLACE ".pt$" ".o" out_f ${out_f})
    set(out_f "${CMAKE_CURRENT_BINARY_DIR}/${out_f}")
    add_custom_command(OUTPUT ${out_f}
      COMMAND ${CMAKE_LINKER} -r -b binary -o ${out_f} ${in_f}
      DEPENDS ${in_f} 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Creating model object file ${out_f}"
      VERBATIM
      )
    list(APPEND result ${out_f})
  endforeach()

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/torch_models.cpp"
    COMMAND ${PYTHON_EXECUTABLE} make_model_cpp.py "${result}" "${CMAKE_CURRENT_BINARY_DIR}/torch_models.cpp"
    DEPENDS ${result} make_model_cpp.py 
    WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Creating ${CMAKE_CURRENT_BINARY_DIR}/torch_models.cpp"
    VERBATIM
  )
  set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

set(MODELS
lib/models/crossdock_default2018_1.pt
lib/models/crossdock_default2018_2.pt
lib/models/crossdock_default2018_3.pt
lib/models/crossdock_default2018_4.pt
lib/models/crossdock_default2018.pt
lib/models/default2017.pt
lib/models/dense_1.pt
lib/models/dense_2.pt
lib/models/dense_3.pt
lib/models/dense_4.pt
lib/models/dense.pt
lib/models/general_default2018_1.pt
lib/models/general_default2018_2.pt
lib/models/general_default2018_3.pt
lib/models/general_default2018_4.pt
lib/models/general_default2018.pt
lib/models/redock_default2018_1.pt
lib/models/redock_default2018_2.pt
lib/models/redock_default2018_3.pt
lib/models/redock_default2018_4.pt
lib/models/redock_default2018.pt
)

model_to_obj(MODEL_OBJS "${MODELS}")

#grid member functions need to be called by arbitrary kernels
set_source_files_properties(lib/grid_gpu.cu PROPERTIES COMPILE_FLAGS -maxrregcount=64)

#test
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.c*)

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
# This is horribly contrived, but I can't figure out a better way -
# for some reason cmake insists on updating the timestamp of version.cpp
# if it is an output (but as a byproduct it doesn't register as a dependency)
# but we only want it to update if the version has actually changed, so
# we chain two custom commands because it seems to work this way. 
add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/version.cpp           
    COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_SOURCE_DIR}/cmake/version.cmake)
add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/_version.cpp PROPERTIES SYMBOLIC TRUE)
            
add_library(gninalibobj OBJECT ${LIB_SRCS})
set_property(TARGET gninalibobj PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(gninalibobj ${MODEL_OBJS})
add_dependencies(gninalibobj libmolgrid)


## Get all properties that cmake supports
execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
## Convert command output into a CMake list
STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)


add_library(gninalib_static STATIC $<TARGET_OBJECTS:gninalibobj>)
set_target_properties(gninalib_static PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib_static  ${Caffe_LINK}  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES}  ${TORCH_LIBRARIES} jsoncpp_lib ZLIB::ZLIB gninalibobj)

add_library(gninalib $<TARGET_OBJECTS:gninalibobj>)
set_target_properties(gninalib PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib  ${Caffe_LINK}  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES}  ${TORCH_LIBRARIES} jsoncpp_lib ZLIB::ZLIB gninalibobj)
# MESSAGE(STATUS "variable is: " ${NVCC_FLAGS_EXTRA})

install(TARGETS gninalib gninalib_static DESTINATION lib)

include_directories(lib)
include_directories(gninavis)
include_directories(${LIBMOLGRID_INCLUDE})

# compile in libgnina - there are enough dependencies to deal with
add_executable(gnina main/main.cpp $<TARGET_OBJECTS:gninalibobj> )
target_link_libraries(gnina  ${Caffe_LINK}  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${LIBMOLGRID_LIBRARY} ${TORCH_LIBRARIES} jsoncpp_lib ZLIB::ZLIB gninalibobj)

add_subdirectory(gninaserver)

add_executable(fromgnina fromgnina/fromgnina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(fromgnina  gninalib )

add_executable(dumpmodel dumpmodel/dumpmodel.cpp)
target_link_libraries(dumpmodel  gninalib )

add_executable(gninagrid gninagrid/gninagrid.cpp gninagrid/molgridder.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(gninagrid  gninalib  ${CUDA_LIBRARIES})

add_executable(gninatyper gninatyper/gninatyper.cpp lib/CommandLine2/CommandLine.cpp )
target_link_libraries(gninatyper   ${Caffe_LINK}  ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${LIBMOLGRID_LIBRARY} ${TORCH_LIBRARIES} jsoncpp_lib gninalib)

add_executable(tognina tognina/tognina.cpp lib/CommandLine2/CommandLine.cpp)
target_link_libraries(tognina   ${Caffe_LINK}  gninalib ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES})

install(TARGETS gnina gninagrid gninatyper fromgnina tognina RUNTIME DESTINATION bin)

# gninavis uses rdkit, which can be a pain to install, so gracefully deal with its absence 
find_package(RDKit)

if(RDKIT_FOUND)
    include_directories(${RDKIT_INCLUDE_DIR})
    set(rdkit_lib_line "-lFileParsers -lSmilesParse -lDepictor -lGraphMol -lRDGeometryLib -lRDGeneral -lSubstructMatch -lSubgraphs -lDataStructs")
    
    add_executable(gninavis gninavis/gninavis.cpp lib/CommandLine2/CommandLine.cpp gninavis/cnn_visualization.cpp gninavis/cnn_visualization.hpp)
    target_link_libraries(gninavis  ${Caffe_LINK}  gninalib ${Boost_LIBRARIES} ${OpenBabel3_LIBRARIES} ${RDKIT_LIBRARIES})
    install(TARGETS gninavis RUNTIME DESTINATION bin)
    
else()
    message(WARNING "RDKIT not found.  gninavis will NOT be built.")
endif()


