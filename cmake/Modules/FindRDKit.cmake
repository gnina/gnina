# FindRDKit.cmake
# Placed in the public domain by NextMove Software in 2013
# Try to find RDKit headers and libraries
# Defines:
#
#  RDKIT_FOUND - system has RDKit
#  RDKIT_INCLUDE_DIR - the RDKit include directory
#  RDKIT_LIBRARIES - Link these to use RDKit
#
# References:
#  
#  http://nextmovesoftware.com/blog/2013/02/04/looking-for-a-c-cheminformatics-toolkit/
#  https://github.com/timvdm/MolDB/blob/master/cmake/modules/FindRDKit.cmake

include(FindPackageHandleStandardArgs)

if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
  # in cache already or user-specified
  find_package_handle_standard_args(RDKit  DEFAULT_MSG
                                  RDKIT_INCLUDE_DIR RDKIT_LIBRARIES)
else()

  if(NOT RDKIT_INCLUDE_DIR)
    if(WIN32)
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
        ${RDBASE}\\Code
        $ENV{RDKIT_INCLUDE_DIR}
        $ENV{RDKIT_INCLUDE_PATH}
        $ENV{RDKIT_BASE}\\Code
        $ENV{RDBASE}\\Code
        C:\\RDKit\\include
        C:\\RDKit\\Code
      )
    else()
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
          ${RDBASE}/Code
          $ENV{RDKIT_INCLUDE_DIR}
          $ENV{RDKIT_INCLUDE_PATH}
          $ENV{RDKIT_BASE}/Code
          $ENV{RDBASE}/Code
          /usr/include/rdkit
          /usr/local/include/rdkit
          /usr/local/rdkit/include/Code
          /usr/local/rdkit/include
          /usr/local/rdkit/Code
          ~/rdkit/Code
      )
    endif()
    if(RDKIT_INCLUDE_DIR)
       message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT RDKIT_LIBRARIES)
      find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers
      PATHS
        ${RDBASE}/lib
        $ENV{RDKIT_LIB_DIR}
        $ENV{RDKIT_LIB_PATH}
        $ENV{RDKIT_LIBRARIES}
        $ENV{RDKIT_BASE}/lib
        $ENV{RDBASE}/lib
        /usr/local/rdkit/lib
        ~/rdkit/lib
        ${RDKIT_LIBRARY_DIR}
        $ENV{LD_LIBRARY_PATH}

       #ignore default path, so search starts with above paths
       NO_DEFAULT_PATH
    )

    #run with default paths this time
    find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers)

    if(FILEPARSERS_LIB)
       GET_FILENAME_COMPONENT(RDKIT_LIBRARY_DIR ${FILEPARSERS_LIB} PATH)
       message(STATUS "Found RDKit libraries at ${RDKIT_LIBRARY_DIR}")

      # Note that the order of the following libraries is significant!!
      find_library(SMILESPARSE_LIB NAMES SmilesParse RDKitSmilesParse
                                   HINTS ${RDKIT_LIBRARY_DIR})
      find_library(DEPICTOR_LIB NAMES Depictor RDKitDepictor
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(GRAPHMOL_LIB NAMES GraphMol RDKitGraphMol
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGEOMETRYLIB_LIB NAMES RDGeometryLib RDKitRDGeometryLib
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGENERAL_LIB NAMES RDGeneral RDKitRDGeneral
                                 HINTS ${RDKIT_LIBRARY_DIR})

      #jhochuli - additional libraries for gninavis
      find_library(SUBSTRUCTMATCH_LIB NAMES SubstructMatch RDKitSubstructMatch
                                 HINTS ${RDKIT_LIBRARY_DIR})
      find_library(SUBGRAPHS_LIB NAMES Subgraphs RDKitSubgraphs
                                 HINTS ${RDKIT_LIBRARY_DIR})
      find_library(DATASTRUCTS_LIB NAMES DataStructs RDKitDataStructs
                                 HINTS ${RDKIT_LIBRARY_DIR})

      set (RDKIT_LIBRARIES ${FILEPARSERS_LIB} ${SMILESPARSE_LIB}
              ${SUBSTRUCTMATCH_LIB} ${GRAPHMOL_LIB} ${RDGEOMETRYLIB_LIB} ${RDGENERAL_LIB}
              ${SUBGRAPHS_LIB} ${DATASTRUCTS_LIB} ${DEPICTOR_LIB}
              )
    endif()
    if(RDKIT_LIBRARIES)
            message(STATUS "Found RDKit library files at ${RDKIT_LIBRARIES}")
    endif()
  endif()

  find_package_handle_standard_args(RDKit  DEFAULT_MSG
                                  RDKIT_INCLUDE_DIR RDKIT_LIBRARIES)
  mark_as_advanced(RDKIT_INCLUDE_DIR RDKIT_LIBRARIES RDKIT_LIBRARY_DIR)
endif()
