cmake_minimum_required(VERSION 3.0)
project(raymond)

enable_language(Fortran)
set(CMAKE_Fortran_COMPILER ifort)

# Set compiler flags
set(FFLAGS "-fopenmp")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
set(CFLAGS "-std=c99")

# Add the source file
add_library(raymond MODULE raymond.f)

# Set the output name of the module
set_target_properties(raymond PROPERTIES OUTPUT_NAME remap)

# Link any necessary libraries
target_link_libraries(raymond PRIVATE iomp5)

# Specify the build directory for f2py
set(F2PY_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Get the Python executable path - to generalize the f2py3.x command
execute_process(
    COMMAND python3 -c "import sys; print(sys.executable)"
    OUTPUT_VARIABLE PYTHON_EXECUTABLE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Compile the Fortran code using f2py
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py  -DF2PY_REPORT_ATEXIT --build-dir ${F2PY_BUILD_DIR} -c -m raymond raymond.f
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE F2PY_RESULT
)

# Generate a list of files matching the pattern
file(GLOB RAYMOND_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*cpython*.so")

# Install each file individually
foreach(FILE ${RAYMOND_FILES})
    get_filename_component(FILENAME ${FILE} NAME)
    install(
        FILES ${FILE}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RENAME ${FILENAME}
    )
endforeach()
