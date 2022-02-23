file(GLOB src_all nc_readwrt_mod.f90 grib2_read_mod.f90)
#set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} )

set(EXE "update_GVF.exe")
add_executable(${EXE} ${src_all} update_GVF.f90  )
target_link_libraries( ${EXE} ${PROJECT_NAME}::ncio
                              NetCDF::NetCDF_Fortran
                              ${GSILIB}
                              g2::g2_4
                              g2tmpl::g2tmpl
                              w3nco::w3nco_4
                              bacio::bacio_4
                              MPI::MPI_Fortran )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})


set(EXE "gen_annual_maxmin_GVF.exe")
add_executable(${EXE} ${src_all} gen_annual_maxmin_GVF.f90)
target_link_libraries( ${EXE} NetCDF::NetCDF_Fortran
                              ${GSILIB}
                              g2::g2_4
                              g2tmpl::g2tmpl
                              w3nco::w3nco_4
                              bacio::bacio_4
                              MPI::MPI_Fortran )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
