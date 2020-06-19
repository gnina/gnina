add_test(NAME pygninacaffe COMMAND ${PYTHON_EXECUTABLE} -m pytest -rsv WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(pygninacaffe PROPERTIES ENVIRONMENT
                    "PYTHONPATH=${CMAKE_SOURCE_DIR}/caffe/python:$ENV{PYTHONPATH}")
