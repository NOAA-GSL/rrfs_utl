list(APPEND src_update_snowice
  grib2_read_snow_mod.f90
  module_imssnow.f90
  update_snowice_fv3lam.f90
  process_imssnow_fv3lam.f90)

add_executable(process_imssnow_fv3lam.exe ${src_update_snowice})
target_link_libraries(process_imssnow_fv3lam.exe PRIVATE ${PROJECT_NAME}::ncio
                                                 PRIVATE NetCDF::NetCDF_Fortran
                                                 PRIVATE g2::g2_4
                                                 PRIVATE g2tmpl::g2tmpl
                                                 PRIVATE bacio::bacio_4
                                                 PRIVATE ip::ip_4
                                                 PRIVATE sp::sp_4
                                                 PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(process_imssnow_fv3lam.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS process_imssnow_fv3lam.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
