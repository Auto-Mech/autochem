project(x2z)
cmake_minimum_required(VERSION 3.0)
find_package(PythonInterp REQUIRED)
if(NOT PYTHON_INSTALL_PATH)
    message(STATUS "PYTHON_INSTALL_PATH is not set")
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c
        "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
        OUTPUT_VARIABLE PYTHON_INSTALL_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Setting PYTHON_INSTALL_PATH: ${PYTHON_INSTALL_PATH}")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
find_package(pybind11 REQUIRED)
add_library(libx2z
    ${PROJECT_SOURCE_DIR}/src/libx2z/atom.cc
    ${PROJECT_SOURCE_DIR}/src/libx2z/chem.cc
    ${PROJECT_SOURCE_DIR}/src/libx2z/d3.cc
    ${PROJECT_SOURCE_DIR}/src/libx2z/linpack.cc
    ${PROJECT_SOURCE_DIR}/src/libx2z/math.cc
    ${PROJECT_SOURCE_DIR}/src/libx2z/units.cc)
add_executable(x2z ${PROJECT_SOURCE_DIR}/src/x2z.cc)
pybind11_add_module(pyx2z SHARED ${PROJECT_SOURCE_DIR}/src/pyx2z.cc)
target_link_libraries(x2z libx2z)
target_link_libraries(pyx2z PRIVATE libx2z)
install(TARGETS x2z DESTINATION bin)
install(TARGETS pyx2z LIBRARY DESTINATION ${PYTHON_INSTALL_PATH})
