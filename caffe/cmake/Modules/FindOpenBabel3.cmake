# FindOpenBabel.cmake
# Try to find Open Babel headers and libraries
# Defines:
#
#  OPENBABEL3_FOUND - system has Open Babel
#  OPENBABEL3_INCLUDE_DIR - the Open Babel include directory
#  OPENBABEL3_LIBRARIES - Link these to use Open Babel
#  IF OPENBABEL_DIR is defined, will look there first

if(OPENBABEL3_INCLUDE_DIR AND OPENBABEL3_LIBRARIES)
  # in cache already or user-specified
  set(OPENBABEL3_FOUND TRUE)

else()

  if(NOT OPENBABEL3_INCLUDE_DIR)
      find_path(OPENBABEL3_INCLUDE_DIR openbabel/obconversion.h
        PATHS
        ${OPENBABEL_DIR}/include/openbabel3
        ${OPENBABEL_DIR}/include
        $ENV{OPENBABEL_INCLUDE_DIR}/openbabel3
        $ENV{OPENBABEL_INCLUDE_DIR}
        $ENV{OPENBABEL_INCLUDE_PATH}/openbabel3
        $ENV{OPENBABEL_INCLUDE_PATH}
        $ENV{OPENBABEL_DIR}/include/openbabel3
        $ENV{OPENBABEL_DIR}/include
        $ENV{OPENBABEL_PATH}/include/openbabel3
        $ENV{OPENBABEL_PATH}/include
        $ENV{OPENBABEL_BASE}/include/openbabel3
        $ENV{OPENBABEL_BASE}/include
        /usr/include/openbabel3
        /usr/include
        /usr/local/include/openbabel3
        /usr/local/include
        /usr/local/openbabel/include/openbabel3
        /usr/local/openbabel/include
        /usr/local/openbabel3/include/openbabel3
        /usr/local/openbabel3/include
        ~/include/openbabel3
        ~/include
      )
    if(OPENBABEL3_INCLUDE_DIR)
      message(STATUS "Found Open Babel include files at ${OPENBABEL3_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT OPENBABEL3_LIBRARIES)
  find_library(OPENBABEL3_LIBRARIES NAMES openbabel openbabel3
      HINTS
      ${OPENBABEL_DIR}/lib
      ${OPENBABEL_DIR}/windows-vc2008/build/src/Release
      ${OPENBABEL_DIR}/build/src/Release
      PATHS
      $ENV{OPENBABEL_DIR}/lib
      $ENV{OPENBABEL_DIR}/windows-vc2008/build/src/Release
      $ENV{OPENBABEL_PATH}/lib
      $ENV{OPENBABEL_BASE}/lib
      /usr/lib
      /usr/local/lib
      ~/lib
      $ENV{LD_LIBRARY_PATH}
    )
    if(OPENBABEL3_LIBRARIES)
      message(STATUS "Found Open Babel library at ${OPENBABEL3_LIBRARIES}")
    endif()
  endif()

  if(OPENBABEL3_INCLUDE_DIR AND OPENBABEL3_LIBRARIES)
    set(OPENBABEL3_FOUND TRUE)
message("Setting openbabel found ${OPENBABEL3_FOUND}")
  endif()

  mark_as_advanced(OPENBABEL3_INCLUDE_DIR OPENBABEL3_LIBRARIES)
endif()
