file(GLOB src_all  netCDFsub_geo.f90 surface_lakesurge_netcdf_mass.f90 )
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} )

set(EXE "lakesurgery.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} wrf_io::wrf_io
                              NetCDF::NetCDF_Fortran
                              ${GSILIB}
                              ${MPI_Fortran_LIBRARIES} )
