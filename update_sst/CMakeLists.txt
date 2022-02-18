### file(GLOB src_all grib2_read_mod.f90 netCDFsub_geo.f90  process_SST.f90  sstGlobal2RR.f90  update_SST_netcdf_mass.f90 update_SST_netcdf_fv3.f90)
file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} )

set(EXE "process_updatesst")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} wrf_io::wrf_io
                              NetCDF::NetCDF_Fortran
                              ${GSILIB}
                              g2::g2_4
                              g2tmpl::g2tmpl
                              w3nco::w3nco_4
                              bacio::bacio_4
                              ${MPI_Fortran_LIBRARIES} )

