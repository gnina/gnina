
set(SERVER_SRCS
Logger.h
MinimizationQuery.cpp
MinimizationQuery.h
QueryManager.cpp
QueryManager.h
Reorienter.h
servercmds.h
server_common.h
server.cpp
../lib/CommandLine2/CommandLine.cpp
)

find_library(CUDA_DEV_RT NAMES cudadevrt PATHS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
cuda_add_executable(gninaserver ${SERVER_SRCS})
target_link_libraries(gninaserver caffe gninalib ${Boost_LIBRARIES}
    ${OPENBABEL2_LIBRARIES} ${CUDA_DEV_RT})

install(TARGETS gninaserver DESTINATION bin)
