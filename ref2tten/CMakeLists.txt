list(APPEND src_ref2tten
  build_missing_REFcone.f90
  convert_lghtn2ref.f90
  convert_stcst2ref.f90
  gsdcloudanalysis_ref2tten.f90
  module_gsi_rfv3io_tten.f90
  module_wrf_error.f90
  pbl_height.f90
  radar_ref2tten.f90
  read_Lightning_cld.f90
  smooth.f90
  vinterp_radar_ref.f90)

add_executable(ref2tten.exe ${src_ref2tten})
target_link_libraries(ref2tten.exe PRIVATE NetCDF::NetCDF_Fortran
                                   PRIVATE gsi::gsi
                                   PRIVATE MPI::MPI_Fortran)

install(
  TARGETS ref2tten.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

