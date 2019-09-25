
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

add_executable(gninaserver ${SERVER_SRCS})
target_link_libraries(gninaserver caffe gninalib ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})

install(TARGETS gninaserver DESTINATION bin)
