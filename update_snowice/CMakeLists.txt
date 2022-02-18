file(GLOB src_all  grib2_read_snow_mod.f90 module_imssnow.f90 update_snowice_fv3lam.f90 process_imssnow_fv3lam.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} )

set(EXE "process_imssnow_fv3lam.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} ${LIBBASE}
                              NetCDF::NetCDF_Fortran
                              ${GSILIB}
                              g2::g2_4
                              g2tmpl::g2tmpl
                              ${LIBPNG}
                              ${LIBJASPER}
                              ${LIBZ}
                              bacio::bacio_4
                              ip::ip_4
                              sp::sp_4
                              ${MPI_Fortran_LIBRARIES} )

