cmake_minimum_required(VERSION 3.0)
project(raymond)

# Path to the Fortran source
set(FCODE "${CMAKE_CURRENT_SOURCE_DIR}/raymond.f")

# Python executable
if(DEFINED ENV{BLENDINGPYTHON})
  set(BPY "$ENV{BLENDINGPYTHON}")
else()
  find_program(BPY python3)
endif()

# Compiler flags
set(CFLAGS "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION -diag-disable=2650 -diag-disable=1478")
set(FFLAGS "-diag-disable=2650 -diag-disable=1478")

# Build command (shell wrapper)
set(F2PY_CMD "
CFLAGS='${CFLAGS}' FFLAGS='${FFLAGS}' FCFLAGS='${FFLAGS}' \
${BPY} -m numpy.f2py -c -m raymond ${FCODE} \
--fcompiler=intelem \
--build-dir \"${CMAKE_CURRENT_BINARY_DIR}\" \
-DF2PY_REPORT_ATEXIT
")

# Execute the build
execute_process(
    COMMAND /bin/sh -c "${F2PY_CMD}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    RESULT_VARIABLE F2PY_RESULT
    OUTPUT_VARIABLE F2PY_OUT
    ERROR_VARIABLE F2PY_ERR
)

# Print output and errors
message(STATUS "f2py stdout:\n${F2PY_OUT}")
message(STATUS "f2py stderr:\n${F2PY_ERR}")

# Fail CMake if f2py failed
if(NOT F2PY_RESULT EQUAL 0)
    message(FATAL_ERROR "f2py failed (exit ${F2PY_RESULT})")
endif()

# Install the resulting .so
file(GLOB_RECURSE SOFILE "${CMAKE_CURRENT_BINARY_DIR}/*raymond*.so")
install(FILES ${SOFILE} DESTINATION ${CMAKE_INSTALL_LIBDIR} RENAME raymond.so)

