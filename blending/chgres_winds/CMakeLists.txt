cmake_minimum_required(VERSION 3.0)
project(chgres_winds)

enable_language(Fortran)
set(CMAKE_Fortran_COMPILER ifort)

# Set compiler flags
set(FFLAGS "-fopenmp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
set(CFLAGS "-std=c99")

# Add the source file
add_library(chgres_winds MODULE chgres_winds.f90)

# Set the output name of the module
set_target_properties(chgres_winds PROPERTIES OUTPUT_NAME chgres_winds)

# Find the path to libiomp5.so
find_library(LIBIOMP5_PATH NAMES libiomp5.so)

# Specify the build directory for f2py
set(F2PY_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Compile the Fortran code using f2py
execute_process(
    COMMAND f2py -DF2PY_REPORT_ATEXIT --build-dir ${F2PY_BUILD_DIR} -c --f90flags='-fopenmp' -L${LIBIOMP5_PATH} -liomp5 -m chgres_winds chgres_winds.f90
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE F2PY_RESULT
)
