### file(GLOB src_all grib2_read_mod.f90 netCDFsub_geo.f90  process_SST.f90  sstGlobal2RR.f90  update_SST_netcdf_mass.f90 update_SST_netcdf_fv3.f90)
file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} )

set(EXE "process_updatesst")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
target_link_libraries( ${EXE} PUBLIC wrf_io::wrf_io
                              PUBLIC NetCDF::NetCDF_Fortran
                              PUBLIC ${GSILIB}
                              PUBLIC g2::g2_4
                              PUBLIC g2tmpl::g2tmpl
                              PUBLIC w3nco::w3nco_4
                              PUBLIC bacio::bacio_4
                              PUBLIC MPI::MPI_Fortran )

if(OpenMP_Fortran_FOUND)
  target_link_libraries( ${EXE} PRIVATE OpenMP::OpenMP_Fortran)
endif()
