
set(LIB_SRCS
atom_constants.cpp
box.cpp
builtinscoring.cpp
cache.cpp
cnn_scorer.cpp
coords.cpp
custom_terms.cpp
everything.cpp
flexinfo.cpp
GninaConverter.cpp
grid.cpp
molgetter.cpp
monte_carlo.cpp
mutate.cpp
my_pid.cpp
naive_non_cache.cpp
nngridder.cpp
non_cache.cpp
obmolopener.cpp
parallel_mc.cpp
parallel_progress.cpp
parse_pdbqt.cpp
pdb.cpp
PDBQTUtilities.cpp
quasi_newton.cpp
quaternion.cpp
random.cpp
result_info.cpp
ssd.cpp
szv_grid.cpp
terms.cpp
weighted_terms.cpp
conf_gpu.cu
gpucode.cu
gridmaker.cu
model.cu
non_cache_gpu.cu
precalculate_gpu.cu
tree_gpu.cu
)

include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${OPENBABEL2_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})

cuda_add_library(gninalib ${LIB_SRCS} OPTIONS ${NVCC_FLAGS_EXTRA} --shared)
set_target_properties(gninalib PROPERTIES OUTPUT_NAME gnina)
target_link_libraries(gninalib caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})
install(TARGETS gninalib DESTINATION lib)

#cuda_add_library(gninalib_static STATIC ${LIB_SRCS} OPTIONS ${NVCC_FLAGS_EXTRA})
#set_target_properties(gninalib_static PROPERTIES OUTPUT_NAME gnina)
#target_link_libraries(gninalib_static caffe ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})
