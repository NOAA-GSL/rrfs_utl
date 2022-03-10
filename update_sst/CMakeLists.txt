list(APPEND src_update_sst
  grib2_read_mod.f90
  module_gsi_rfv3io_sst.f90
  netCDFsub_geo.f90
  process_SST.f90
  sstGlobal2RR.f90
  update_SST_netcdf_fv3.f90
  update_SST_netcdf_mass.f90)

add_executable(process_updatesst.exe ${src_update_sst})
target_link_libraries(process_updatesst.exe PRIVATE wrf_io::wrf_io
                                            PRIVATE NetCDF::NetCDF_Fortran
                                            PRIVATE gsi::gsi
                                            PRIVATE g2::g2_4
                                            PRIVATE g2tmpl::g2tmpl
                                            PRIVATE w3nco::w3nco_4
                                            PRIVATE bacio::bacio_4
                                            PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(process_updatesst.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS process_updatesst.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
