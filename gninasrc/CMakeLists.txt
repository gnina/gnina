
find_package(CUDA REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options system iostreams timer thread serialization filesystem date_time regex)

add_compile_options(-std=c++11)
list(APPEND CUDA_NVCC_FLAGS -std=c++11)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  list(APPEND CUDA_NVCC_FLAGS -O3)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  list(APPEND CUDA_NVCC_FLAGS -O0)
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  list(APPEND CUDA_NVCC_FLAGS -g)
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
#lib gets its own CMakeLists.txt
add_subdirectory(lib)

#everything else dealt with here

include_directories(lib)

cuda_add_executable(gnina main/main.cpp)
target_link_libraries(gnina caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

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
