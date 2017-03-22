# FindOpenBabel.cmake
# Try to find Open Babel headers and libraries
# Defines:
#
#  OPENBABEL2_FOUND - system has Open Babel
#  OPENBABEL2_INCLUDE_DIR - the Open Babel include directory
#  OPENBABEL2_LIBRARIES - Link these to use Open Babel
#  IF OPENBABEL_DIR is defined, will look there first

if(OPENBABEL2_INCLUDE_DIR AND OPENBABEL2_LIBRARIES)
  # in cache already or user-specified
  set(OPENBABEL2_FOUND TRUE)

else()

  if(NOT OPENBABEL2_INCLUDE_DIR)
      find_path(OPENBABEL2_INCLUDE_DIR openbabel/obconversion.h
        PATHS
        ${OPENBABEL_DIR}/include/openbabel-2.0
        ${OPENBABEL_DIR}/include
        $ENV{OPENBABEL_INCLUDE_DIR}/openbabel-2.0
        $ENV{OPENBABEL_INCLUDE_DIR}
        $ENV{OPENBABEL_INCLUDE_PATH}/openbabel-2.0
        $ENV{OPENBABEL_INCLUDE_PATH}
        $ENV{OPENBABEL_DIR}/include/openbabel-2.0
        $ENV{OPENBABEL_DIR}/include
        $ENV{OPENBABEL_PATH}/include/openbabel-2.0
        $ENV{OPENBABEL_PATH}/include
        $ENV{OPENBABEL_BASE}/include/openbabel-2.0
        $ENV{OPENBABEL_BASE}/include
        /usr/include/openbabel-2.0
        /usr/include
        /usr/local/include/openbabel-2.0
        /usr/local/include
        /usr/local/openbabel/include/openbabel-2.0
        /usr/local/openbabel/include
        /usr/local/openbabel-2.0/include/openbabel-2.0
        /usr/local/openbabel-2.0/include
        ~/include/openbabel-2.0
        ~/include
      )
    if(OPENBABEL2_INCLUDE_DIR)
      message(STATUS "Found Open Babel include files at ${OPENBABEL2_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT OPENBABEL2_LIBRARIES)
  find_library(OPENBABEL2_LIBRARIES NAMES openbabel openbabel-2
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
    if(OPENBABEL2_LIBRARIES)
      message(STATUS "Found Open Babel library at ${OPENBABEL2_LIBRARIES}")
    endif()
  endif()

  if(OPENBABEL2_INCLUDE_DIR AND OPENBABEL2_LIBRARIES)
    set(OPENBABEL2_FOUND TRUE)
message("Setting openbabel found ${OPENBABEL2_FOUND}")
  endif()

  mark_as_advanced(OPENBABEL2_INCLUDE_DIR OPENBABEL2_LIBRARIES)
endif()

