cmake_minimum_required(VERSION 3.0)
project(raymond)

enable_language(Fortran)
set(CMAKE_Fortran_COMPILER ifort)

# Set compiler flags
set(FFLAGS "-fopenmp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")

# Add the source file
add_library(raymond MODULE raymond.f)

# Link any necessary libraries
target_link_libraries(raymond PRIVATE iomp5)

# Set the output name of the module
set_target_properties(raymond PROPERTIES OUTPUT_NAME remap)

# Specify the build directory for f2py
set(F2PY_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Compile the Fortran code using f2py
execute_process(
    COMMAND f2py3.8 -DF2PY_REPORT_ATEXIT --build-dir ${F2PY_BUILD_DIR} -c -m raymond raymond.f
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE F2PY_RESULT
)

