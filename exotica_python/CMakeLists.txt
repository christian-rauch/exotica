cmake_minimum_required(VERSION 3.0.2)
project(exotica_python)

find_package(PkgConfig REQUIRED)
pkg_check_modules(MSGPACK QUIET msgpack)

if(MSGPACK_FOUND)
  add_definitions(-DMSGPACK_FOUND)
endif()

# Python2 uses char* and Python3 const char* for arguments so explicit casting
# of string literals would require define flags. Alternatively, deactivate
# warnings for now
add_definitions(-Wno-write-strings)

find_package(catkin REQUIRED COMPONENTS
  exotica_core
  geometry_msgs
  moveit_msgs
  pybind11_catkin
  shape_msgs
)

catkin_package(
  LIBRARIES
  CATKIN_DEPENDS exotica_core pybind11_catkin shape_msgs moveit_msgs geometry_msgs
)

pybind11_add_module(_pyexotica MODULE src/pyexotica.cpp)
target_link_libraries(_pyexotica PRIVATE ${PYTHON_LIBRARIES} ${catkin_LIBRARIES})
target_include_directories(_pyexotica PUBLIC ${catkin_INCLUDE_DIRS} ${pybind11_catkin_INCLUDE_DIRS}/pybind11_catkin) # the latter is a work-around until 2.4.3 is released
add_dependencies(_pyexotica ${catkin_EXPORTED_TARGETS})
set_target_properties(_pyexotica PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CATKIN_GLOBAL_PYTHON_DESTINATION}/pyexotica)
set(PYTHON_LIB_DIR ${CATKIN_DEVEL_PREFIX}/${CATKIN_GLOBAL_PYTHON_DESTINATION}/pyexotica)
add_custom_command(TARGET _pyexotica
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_LIB_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_pyexotica> ${PYTHON_LIB_DIR}/_pyexotica.so
    WORKING_DIRECTORY ${CATKIN_DEVEL_PREFIX}
COMMENT "Copying library files to python directory")

catkin_python_setup()

install(TARGETS _pyexotica LIBRARY DESTINATION ${CATKIN_GLOBAL_PYTHON_DESTINATION}/pyexotica)

if(CATKIN_ENABLE_TESTING)
  catkin_add_nosetests(test/test_no_unknown_initializer_types.py)
  catkin_add_nosetests(test/test_box_qp.py)
endif()
