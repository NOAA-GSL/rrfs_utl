file(GLOB src_all  grib2_read_snow_mod.f90 module_imssnow.f90 update_snowice_fv3lam.f90 process_imssnow_fv3lam.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${IPINC} ${G2INC} ${NETCDFINC} )

set(EXE "process_imssnow_fv3lam.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} ${LIBBASE} ${LIBIP} NetCDF::NetCDF_Fortran ${GSILIB} ${LIBG2} ${LIBG2TMPL} ${LIBW3NCO} ${LIBPNG} ${LIBJASPER} ${LIBZ} ${LIBBACIO} ${LIBSP} ${MPI_Fortran_LIBRARIES} )

