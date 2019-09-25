# FindLibMolGrid.cmake
# Try to find libmolgrid headers and libraries
# Defines:
#
#  LIBMOLGRID_FOUND - system has libmolgrid
#  LIBMOLGRID_INCLUDE_DIR - the libmolgrid include directory
#  LIBMOLGRID_LIBRAR& - Link to use libmolgrid
#


include(FindPackageHandleStandardArgs)

if(LIBMOLGRID_INCLUDE_DIR AND LIBMOLGRID_LIBRARY)
  # in cache already or user-specified
  find_package_handle_standard_args(libmolgrid  DEFAULT_MSG
                                  LIBMOLGRID_INCLUDE_DIR LIBMOLGRID_LIBRARY)
else()

  if(NOT LIBMOLGRID_INCLUDE_DIR)
    find_path(LIBMOLGRID_INCLUDE_DIR libmolgrid/libmolgrid.h
    PATHS
      $ENV{LIBMOLGRID_INCLUDE_DIR}
      ${CMAKE_INSTALL_PREFIX}/include
      /usr/include/
      /usr/local/include
      )
    if(LIBMOLGRID_INCLUDE_DIR)
       message(STATUS "Found libmolgrid include files at ${LIBMOLGRID_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT LIBMOLGRID_LIBRARY)
      find_library(LIBMOLGRID_LIBRARY NAMES molgrid
      PATHS
        $ENV{LIBMOLGRID_LIBRARY_DIR}
        ${CMAKE_INSTALL_PREFIX}/lib
        /usr/lib
        /usr/local/lib
        $ENV{LD_LIBRARY_PATH}
     )
    if(LIBMOLGRID_LIBRARY)
        message(STATUS "Found libmolgrid library at ${LIBMOLGRID_LIBRARY}")
    endif()

  endif()

  find_package_handle_standard_args(libmolgrid  DEFAULT_MSG
                                  LIBMOLGRID_INCLUDE_DIR LIBMOLGRID_LIBRARY)
  mark_as_advanced(LIBMOLGRID_INCLUDE_DIR LIBMOLGRID_LIBRARY LIBMOLGRID_LIBRARY_DIR)
endif()
