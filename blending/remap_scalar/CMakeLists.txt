cmake_minimum_required(VERSION 3.0)
project(remap_scalar)

# Find the path to libiomp5.so
find_library(LIBIOMP5_PATH NAMES libiomp5.so)
set(FCODE ${CMAKE_CURRENT_SOURCE_DIR}/remap_scalar.f90)

# Compile the Fortran code using f2py
execute_process(
  COMMAND $ENV{BLENDINGPYTHON} -m numpy.f2py -DF2PY_REPORT_ATEXIT --fcompiler=intelem --build-dir ${CMAKE_CURRENT_BINARY_DIR} -c --f90flags='-fopenmp' -L${LIBIOMP5_PATH} -liomp5 -m remap_scalar ${FCODE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    RESULT_VARIABLE F2PY_RESULT
)

#install remap_scalar.so
file(GLOB_RECURSE SOFILE "${CMAKE_CURRENT_BINARY_DIR}/*remap_scalar*.so")
install(FILES ${SOFILE} DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME remap_scalar.so)
