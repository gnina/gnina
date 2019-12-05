
#Prep ourselves for compiling boost
find_package(Boost COMPONENTS unit_test_framework system REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})
include_directories(../gninasrc/lib)
include_directories(../gninavis)
include_directories(${LIBMOLGRID_INCLUDE})

# find python
find_package(PythonInterp)
if(NOT PYTHONINTERP_FOUND)
    message(FATAL_ERROR "ERROR: Python interpreter not found. Cannot continue - sorry.")
else()
    #check for pytest
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pytest --version
      OUTPUT_VARIABLE PYTEST_output
      ERROR_VARIABLE  PYTEST_error
      RESULT_VARIABLE PYTEST_result)
    if(${PYTEST_result})
      message(SEND_ERROR "Pytest package not available: ${PYTEST_error}")
    endif()
endif()

add_subdirectory(caffe)
add_subdirectory(gninagrid)
add_subdirectory(gnina)
