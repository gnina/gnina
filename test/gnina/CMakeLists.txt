


#get all cpp files
set( TEST_SRCS
 test_cache.cu
 test_cache.h
 test_cnn.cpp
 test_cnn.h
 test_gpucode.cpp
 test_gpucode.h
 test_runner.cpp
 test_tree.h
 test_tree.cu
 test_utils.h
)

find_package(OpenMP)

#todo - factor into individual tests?
add_executable(gninacheck ${TEST_SRCS})
target_link_libraries(gninacheck gninalib_static caffe ${Boost_LIBRARIES} ${OPENBABEL_LIBRARIES} ${RDKIT_LIBRARIES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(gninacheck OpenMP::OpenMP_CXX)
endif()

if(Boost_MINOR_VERSION  VERSION_GREATER 64)
    add_test(NAME gninacheck COMMAND gninacheck -- --n_iters 1)
else()
    add_test(NAME gninacheck COMMAND gninacheck --n_iters=1)
endif()

add_test(NAME gninamin COMMAND ./test_min.py $<TARGET_FILE:gnina> WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME gninaflex COMMAND ./test_flex.py $<TARGET_FILE:gnina> WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})