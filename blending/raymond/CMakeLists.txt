cmake_minimum_required(VERSION 3.0)
project(raymond)

set(FCODE ${CMAKE_CURRENT_SOURCE_DIR}/raymond.f)

# Compile the Fortran code using f2py
execute_process(
    COMMAND $ENV{BLENDINGPYTHON} -m numpy.f2py  -DF2PY_REPORT_ATEXIT --fcompiler=intelem --build-dir ${CMAKE_CURRENT_BINARY_DIR} -c -m raymond ${FCODE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    RESULT_VARIABLE F2PY_RESULT
)

#install raymond.so
file(GLOB_RECURSE SOFILE "${CMAKE_CURRENT_BINARY_DIR}/*raymond*.so")
install(FILES ${SOFILE} DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME raymond.so)
