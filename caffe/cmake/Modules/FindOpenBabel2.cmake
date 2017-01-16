# - Try to find OpenBabel2
# Once done this will define
#
#  OPENBABEL2_FOUND - system has OpenBabel2
#  OPENBABEL2_INCLUDE_DIR - the OpenBabel2 include directory
#  OPENBABEL2_LIBRARIES - Link these to use OpenBabel2
#
# Copyright (c) 2006, 2007 Carsten Niehaus, <cniehaus@gmx.de>
# Copyright (C) 2008 Marcus D. Hanwell <marcus@cryos.org>
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if (OPENBABEL2_INCLUDE_DIR AND OPENBABEL2_LIBRARIES AND OPENBABEL2_VERSION_MET)
  # in cache already
  set(OPENBABEL2_FOUND TRUE)

else (OPENBABEL2_INCLUDE_DIR AND OPENBABEL2_LIBRARIES AND OPENBABEL2_VERSION_MET)
  if(NOT WIN32)

    # Use the newer PkgConfig stuff
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(OPENBABEL2 openbabel-2.0>=2.2.0)

    # Maintain backwards compatibility with previous version of module
    if(OPENBABEL2_FOUND STREQUAL "1")
      set(OPENBABEL2_VERSION_MET TRUE)
      set(OPENBABEL2_INCLUDE_DIR ${OPENBABEL2_INCLUDE_DIRS})
    endif(OPENBABEL2_FOUND STREQUAL "1")

  else(NOT WIN32)
    set(OPENBABEL2_VERSION_MET TRUE)
  endif(NOT WIN32)

  if(OPENBABEL2_VERSION_MET)

    if(WIN32)
      if(NOT OPENBABEL2_INCLUDE_DIR)
        find_path(OPENBABEL2_INCLUDE_DIR openbabel-2.0/openbabel/obconversion.h
          PATHS
          ${_obIncDir}
          ${GNUWIN32_DIR}/include
          $ENV{OPENBABEL2_INCLUDE_DIR}
        )
        if(OPENBABEL2_INCLUDE_DIR)
          set(OPENBABEL2_INCLUDE_DIR ${OPENBABEL2_INCLUDE_DIR}/openbabel-2.0)
        endif(OPENBABEL2_INCLUDE_DIR)
      endif(NOT OPENBABEL2_INCLUDE_DIR)
    endif(WIN32)

    find_library(OPENBABEL2_LIBRARIES NAMES openbabel openbabel-2
      PATHS
      ${_obLinkDir}
      ${GNUWIN32_DIR}/lib
      $ENV{OPENBABEL2_LIBRARIES}
    )
  endif(OPENBABEL2_VERSION_MET)

  if(OPENBABEL2_INCLUDE_DIR AND OPENBABEL2_LIBRARIES AND OPENBABEL2_VERSION_MET)
    set(OPENBABEL2_FOUND TRUE)
  endif(OPENBABEL2_INCLUDE_DIR AND OPENBABEL2_LIBRARIES AND OPENBABEL2_VERSION_MET)

  if (OPENBABEL2_FOUND)
    if (NOT OpenBabel2_FIND_QUIETLY)
      message(STATUS "Found OpenBabel 2.2 or later: ${OPENBABEL2_LIBRARIES}")
    endif (NOT OpenBabel2_FIND_QUIETLY)
  else (OPENBABEL2_FOUND)
    if (OpenBabel2_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find OpenBabel 2.2 or later ")
    endif (OpenBabel2_FIND_REQUIRED)
  endif (OPENBABEL2_FOUND)

  mark_as_advanced(OPENBABEL2_INCLUDE_DIR OPENBABEL2_LIBRARIES)

endif (OPENBABEL2_INCLUDE_DIR AND OPENBABEL2_LIBRARIES AND OPENBABEL2_VERSION_MET)

# Search for Open Babel2 executable
if(OPENBABEL2_EXECUTABLE)

  # in cache already
  set(OPENBABEL2_EXECUTABLE_FOUND TRUE)

else(OPENBABEL2_EXECUTABLE)
  find_program(OPENBABEL2_EXECUTABLE NAMES babel
    PATHS
    [HKEY_CURRENT_USER\\SOFTWARE\\OpenBabel\ 2.2.0]
    $ENV{OPENBABEL2_EXECUTABLE}
  )

  if(OPENBABEL2_EXECUTABLE)
    set(OPENBABEL2_EXECUTABLE_FOUND TRUE)
  endif(OPENBABEL2_EXECUTABLE)

  if(OPENBABEL2_EXECUTABLE_FOUND)
    message(STATUS "Found OpenBabel2 executable: ${OPENBABEL2_EXECUTABLE}")
  endif(OPENBABEL2_EXECUTABLE_FOUND)

endif(OPENBABEL2_EXECUTABLE)
