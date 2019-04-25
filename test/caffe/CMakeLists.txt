


add_test(NAME pygninacaffe COMMAND ${PYTHON_EXECUTABLE} -m pytest -m "gnina specific caffe tests" -rsv WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
