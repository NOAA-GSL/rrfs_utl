cmake_minimum_required(VERSION 3.0)
project(remap_dwinds)

enable_language(Fortran)
set(CMAKE_Fortran_COMPILER ifort)

# Set compiler flags
set(FFLAGS "-fopenmp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")

# Add the source file
add_library(remap_dwinds MODULE remap_dwinds.f90)

# Set the output name of the module
set_target_properties(remap_dwinds PROPERTIES OUTPUT_NAME remap_dwinds)

# Specify the build directory for f2py
set(F2PY_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Compile the Fortran code using f2py
execute_process(
    COMMAND f2py3.8 -DF2PY_REPORT_ATEXIT --build-dir ${F2PY_BUILD_DIR} -c --f90flags='-fopenmp' -L/pe/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin/libiomp5.so -liomp5 -m remap_dwinds remap_dwinds.f90
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE F2PY_RESULT
)

if(F2PY_RESULT EQUAL 0)
    message(STATUS "Fortran code compiled successfully with f2py")
else()
    message(FATAL_ERROR "Failed to compile Fortran code with f2py")
endif()

# Clean target
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove -f *.so *.mod *.o
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

