list(APPEND src_nssl_mosaic
  module_read_NSSL_mosaic.f90
  netCDFsub.f90
  netCDFsub_geo.f90
  process_NSSL_mosaic.f90
  read_ncep_binary.f90
  read_nssl_binary.f90
  write_bufr_ref.f90
  read_grib2_mod.f90)

#list(APPEND src_nssl_mosaic_enkf
#  DART.f90
#  netCDFsub.f90
#  netCDFsub_geo.f90
#  process_NSSL_mosaic_to_dart_and_netcdf.f90
#  read_grib2_mod.f90
#  read_ncep_binary.f90
#  read_nssl_binary.f90
#  write_bufr_ref.f90
#  write_netcdf_ref.f90)

add_executable(process_NSSL_mosaic.exe ${src_nssl_mosaic})
target_link_libraries(process_NSSL_mosaic.exe PRIVATE ${PROJECT_NAME}::ncio
                                              PRIVATE NetCDF::NetCDF_Fortran
                                              PRIVATE gsi::gsi
                                              PRIVATE g2::g2_4
                                              PRIVATE g2tmpl::g2tmpl
                                              PRIVATE bufr::bufr_4
                                              PRIVATE bacio::bacio_4
                                              PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(process_NSSL_mosaic.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS process_NSSL_mosaic.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

#add_executable(process_NSSL_mosaic_enkf.exe ${src_nssl_mosaic_enkf})
#target_link_libraries(process_NSSL_mosaic_enkf.exe PRIVATE NetCDF::NetCDF_Fortran
#                                                   PRIVATE g2::g2_4
#                                                   PRIVATE g2tmpl::g2tmpl
#                                                   PRIVATE bufr::bufr_4
#                                                   PRIVATE bacio::bacio_4
#                                                   PRIVATE MPI::MPI_Fortran)
#
#if(OpenMP_Fortran_FOUND)
#  target_link_libraries(process_NSSL_mosaic_enkf.exe PRIVATE OpenMP::OpenMP_Fortran)
#endif()
#
#install(TARGETS process_NSSL_mosaic_enkf.exe
#        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
