list(APPEND src_update_gvf
  nc_readwrt_mod.f90
  grib2_read_mod.f90)

add_executable(update_GVF.exe ${src_update_gvf} update_GVF.f90)
target_link_libraries(update_GVF.exe PRIVATE ${PROJECT_NAME}::ncio
                                     PRIVATE NetCDF::NetCDF_Fortran
                                     PRIVATE g2::g2_4
                                     PRIVATE g2tmpl::g2tmpl
                                     PRIVATE w3nco::w3nco_4
                                     PRIVATE bacio::bacio_4
                                     PRIVATE MPI::MPI_Fortran)

add_executable(gen_annual_maxmin_GVF.exe ${src_update_gvf} gen_annual_maxmin_GVF.f90)
target_link_libraries(gen_annual_maxmin_GVF.exe PRIVATE NetCDF::NetCDF_Fortran
                                                PRIVATE g2::g2_4
                                                PRIVATE g2tmpl::g2tmpl
                                                PRIVATE w3nco::w3nco_4
                                                PRIVATE bacio::bacio_4
                                                PRIVATE MPI::MPI_Fortran)

install(
  TARGETS update_GVF.exe gen_annual_maxmin_GVF.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
