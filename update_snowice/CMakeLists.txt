file(GLOB src_all  grib2_read_snow_mod.f90 module_imssnow.f90 update_snowice_fv3lam.f90 process_imssnow_fv3lam.f90)

include_directories( ${GSIINC} )

set(EXE "process_imssnow_fv3lam.exe")
add_executable(${EXE} ${src_all})
target_link_libraries( ${EXE} PUBLIC ${PROJECT_NAME}::ncio
                              PUBLIC NetCDF::NetCDF_Fortran
                              PUBLIC ${GSILIB}
                              PUBLIC g2::g2_4
                              PUBLIC g2tmpl::g2tmpl
                              PUBLIC bacio::bacio_4
                              PUBLIC ip::ip_4
                              PUBLIC sp::sp_4
                              PUBLIC MPI::MPI_Fortran )

if(OpenMP_Fortran_FOUND)
  target_link_libraries( ${EXE} PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
