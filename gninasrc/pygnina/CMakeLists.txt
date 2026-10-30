
if(BUILD_CONTAINED_PYTHON)
  # create a python library with as few runtime dependencies as possible
  set(Boost_USE_STATIC_LIBS ON)
endif()

# boost_system is header-only since Boost 1.69 and has no cmake config
# Use the same Python3 interpreter found at the top level so the requested
# boost_python component matches the interpreter used everywhere else in the
# build (mixing this with the deprecated PythonInterp/PythonLibs modules can
# pick a different Python and ask Boost for a component it doesn't have).
find_package( Boost COMPONENTS filesystem python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR} REQUIRED )
find_package(Python3 COMPONENTS NumPy QUIET)
if(Python3_NumPy_FOUND)
  set(PYTHON_NUMPY_INCLUDE_DIR ${Python3_NumPy_INCLUDE_DIRS})
else()
  # Fall back to asking the interpreter (works even without cmake NumPy module)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE PYTHON_NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
#numpy header has unused function
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wno-unused-function")
endif()


set( PYGNINA_SOURCES
 python_streambuf.h
 bindings.cpp
)

set( PYGNINA_PY
 __init__.py
)



file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pygnina/")

# link

# create a python library 
if(NOT APPLE)
  add_link_options("LINKER:--exclude-libs,ALL")  # avoid exporting symbols (GNU ld only)
endif()
add_library(pygnina SHARED ${PYGNINA_SOURCES})
target_link_libraries(pygnina ${Boost_LIBRARIES} ${CUDA_LIBRARIES} gninalib_static )


# Suppress prefix "lib" because Python does not allow this prefix
set_target_properties(pygnina PROPERTIES
									PREFIX ""
									LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pygnina/")

if(APPLE)
  # On macOS the linker requires all symbols to be resolved at link time.
  # Python extension modules intentionally leave Python/Boost.Python symbols
  # unresolved (they are provided by the host interpreter at load time).
  # -undefined dynamic_lookup defers resolution to runtime, matching Linux behaviour.
  # SUFFIX ".so" is required too: CPython's import machinery doesn't recognize
  # the default macOS shared-library extension (.dylib) as a loadable module.
  set_target_properties(pygnina PROPERTIES
                                  LINK_FLAGS "-undefined dynamic_lookup"
                                  SUFFIX ".so")
endif()
									
                  
# copy module code
foreach(file ${PYGNINA_PY})
 configure_file(${file}  ${CMAKE_CURRENT_BINARY_DIR}/pygnina/)
endforeach()



# Determine site-packages destination using the selected Python interpreter
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('platlib', vars={'platbase': '${CMAKE_INSTALL_PREFIX}', 'base': '${CMAKE_INSTALL_PREFIX}'}))"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

install(TARGETS pygnina
  LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/pygnina"
)
install(FILES ${PYGNINA_PY}
  DESTINATION "${PYTHON_SITE_PACKAGES}/pygnina"
)
