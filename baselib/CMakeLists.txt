file(GLOB src_all kinds.f90  module_map_utils.f90  module_ncio.f90  module_time.f90)
set_source_files_properties(${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

add_library(kinds STATIC kinds.f90)
add_library(${PROJECT_NAME}::kinds ALIAS kinds)

add_library(ncio STATIC module_ncio.f90)
add_library(${PROJECT_NAME}::ncio ALIAS ncio)

add_library (timelib STATIC module_time.f90)
add_library(${PROJECT_NAME}::timelib ALIAS timelib)

add_library (map_util STATIC kinds.f90 module_map_utils.f90)
add_library(${PROJECT_NAME}::map_util ALIAS map_util)

add_library (esggrid_util STATIC module_esggrid_util.f90)
add_library(${PROJECT_NAME}::esggrid_util ALIAS esggrid_util)
target_include_directories(esggrid_util PUBLIC ./regional_esg_grid.fd)
