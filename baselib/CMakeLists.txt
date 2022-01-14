cmake_minimum_required(VERSION 2.8)

file(GLOB src_all kinds.f90  module_map_utils.f90  module_ncio.f90  module_time.f90)
set_source_files_properties(${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${NETCDFINC} )

ADD_LIBRARY (kinds STATIC kinds.f90)
ADD_LIBRARY (ncio STATIC module_ncio.f90)
ADD_LIBRARY (timelib STATIC module_time.f90)
ADD_LIBRARY (map_util STATIC kinds.f90 module_map_utils.f90)
#ADD_LIBRARY (gsi_map_util STATIC kinds.f90 gridmod_gsi_map_utils.f90)
add_subdirectory(regional_esg_grid.fd)
