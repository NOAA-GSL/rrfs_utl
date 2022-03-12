cmake_minimum_required(VERSION 2.8)

file(GLOB src_all unfill_fv3lam_grid2.f90 gsl_adjust_soiltq_mod.f90 module_bkio_fv3lam.f90 genqsat.f90 constants.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${NETCDFINC} )

set(EXE "adjust_soiltq.exe")
add_executable(${EXE} adjust_soiltq.f90 ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} ${LIBBASE} ${LIBIP} ${NETCDFLIB} ${GSILIB} ${LIBG2} ${LIBG2TMPL} ${LIBW3NCO} ${LIBPNG} ${LIBJASPER} ${LIBZ} ${MPILIB} ${LIBBACIO} ${LIBSP} ${LIBIOMP5})
