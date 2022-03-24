cmake_minimum_required(VERSION 2.8)

file(GLOB src_all module_update_bc.f90 module_io_fv3lam_bdy.f90 module_io_fv3lam_bk.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${NETCDFINC} )

set(EXE "update_bc.exe")
add_executable(${EXE} update_bc.f90 ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} ${LIBBASE} ${NETCDFLIB} ${LIBW3NCO} ${LIBPNG} ${LIBJASPER} ${LIBZ} ${MPILIB} ${LIBBACIO} ${LIBIOMP5})
